#!/usr/bin/env python3
"""
make_tet_mesh.py  -  closed triangle surface (.glb)  ->  tetrahedral mesh for the tensor-mass simulator.

Method (no external meshing library needed, only numpy):
  1. read the surface from the .glb (chosen mesh parts, welded)
  2. rotate it into the simulator frame (+x towards the camera/probe, +z up) and scale it
  3. lay a regular grid of cell size h over it, split every cube into 6 tetrahedra (Kuhn split,
     the same split as makeBlockMesh() in tensor_mass.h, so neighbouring cubes conform)
  4. keep the tetrahedra whose centroid is inside the surface (ray-casting, majority of 3 rays)
  5. keep only the largest face-connected component
  6. clamp the vertices of the back layer (smallest x)  -> boundary condition
  7. write a C++ header with the arrays (so nothing has to be loaded from disk at run time)

Usage: make_tet_mesh.py <organ> <h> <out.h>      organ: liver | ear | bone
"""
import json, struct, sys
from collections import Counter, defaultdict
import numpy as np

import os
SRC = os.environ.get("ORGAN_GLB_DIR", ".") + "/"   # folder holding the .glb files

def read_glb(path):
    b = open(path, "rb").read()
    total = struct.unpack("<I", b[8:12])[0]; off = 12; js = bn = None
    while off < total:
        cl, ct = struct.unpack("<II", b[off:off + 8]); d = b[off + 8:off + 8 + cl]; off += 8 + cl
        if ct == 0x4E4F534A: js = json.loads(d)
        elif ct == 0x004E4942: bn = d
    return js, bn

def accessor(js, bn, i):
    a = js["accessors"][i]; bv = js["bufferViews"][a["bufferView"]]
    dt = {5126: "f4", 5123: "u2", 5125: "u4"}[a["componentType"]]
    nc = {"SCALAR": 1, "VEC3": 3}[a["type"]]
    o = bv.get("byteOffset", 0) + a.get("byteOffset", 0)
    arr = np.frombuffer(bn, dtype=dt, count=a["count"] * nc, offset=o)
    return arr.reshape(-1, nc) if nc > 1 else arr

def surface(path, names=None):
    js, bn = read_glb(path)
    Vs, Fs, o = [], [], 0
    for m in js["meshes"]:
        if names is not None and m["name"] not in names: continue
        for pr in m["primitives"]:
            V = accessor(js, bn, pr["attributes"]["POSITION"]).astype(float)
            F = accessor(js, bn, pr["indices"]).reshape(-1, 3).astype(int)
            Vs.append(V); Fs.append(F + o); o += len(V)
    V = np.vstack(Vs); F = np.vstack(Fs)
    q = np.round(V / 1e-5).astype(np.int64)
    u, first, inv = np.unique(q, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel(); Vw = V[first]; F = inv[F]
    F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])]
    return Vw, F

def edge_use(F):
    e = np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    return Counter(Counter(map(tuple, e)).values())

def signed_volume(V, F):
    return np.einsum("ij,ij->i", V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6.0

class RayInside:
    """point-in-closed-surface by ray casting along +x (triangles binned on the y-z plane)"""
    def __init__(self, V, F, cell):
        self.V, self.F, self.cell = V, F, cell
        T = V[F]                                   # (n,3,3)
        self.lo = T[:, :, 1:].min(axis=1); self.hi = T[:, :, 1:].max(axis=1)
        self.y0 = V[:, 1].min() - 1e-3; self.z0 = V[:, 2].min() - 1e-3
        self.bins = defaultdict(list)
        iy0 = ((self.lo[:, 0] - self.y0) / cell).astype(int); iy1 = ((self.hi[:, 0] - self.y0) / cell).astype(int)
        iz0 = ((self.lo[:, 1] - self.z0) / cell).astype(int); iz1 = ((self.hi[:, 1] - self.z0) / cell).astype(int)
        for t in range(len(F)):
            for a in range(iy0[t], iy1[t] + 1):
                for b in range(iz0[t], iz1[t] + 1): self.bins[(a, b)].append(t)
        self.bins = {k: np.array(v) for k, v in self.bins.items()}

    def crossings(self, p):
        key = (int((p[1] - self.y0) / self.cell), int((p[2] - self.z0) / self.cell))
        c = self.bins.get(key)
        if c is None: return 0
        A = self.V[self.F[c, 0]]; B = self.V[self.F[c, 1]]; C = self.V[self.F[c, 2]]
        # barycentric coordinates of (py,pz) in the projected triangle
        d = (B[:, 1] - A[:, 1]) * (C[:, 2] - A[:, 2]) - (C[:, 1] - A[:, 1]) * (B[:, 2] - A[:, 2])
        ok = np.abs(d) > 1e-18
        py = p[1] - A[:, 1]; pz = p[2] - A[:, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            u = (py * (C[:, 2] - A[:, 2]) - (C[:, 1] - A[:, 1]) * pz) / d
            v = ((B[:, 1] - A[:, 1]) * pz - py * (B[:, 2] - A[:, 2])) / d
        hit = ok & (u >= 0) & (v >= 0) & (u + v < 1)
        x = A[:, 0] + u * (B[:, 0] - A[:, 0]) + v * (C[:, 0] - A[:, 0])
        return int(np.count_nonzero(hit & (x > p[0])))

    def inside(self, p, rng):
        votes = 0
        for _ in range(3):
            q = p + np.array([0.0, rng.uniform(-1, 1), rng.uniform(-1, 1)]) * self.cell * 1e-3
            votes += self.crossings(q) & 1
        return votes >= 2

PRESETS = {
    # liver: capsule + bare-area patch together form the closed organ surface; the file is in metres.
    # glTF frame: +z anterior, +y superior, +x patient-left  ->  sim (x,y,z) = (gz, gx, gy):
    # the anterior surface faces the probe, z is up; the posterior layer is clamped.
    "liver": dict(file="8b4c6600-liver.glb",
                  parts={"VH_F_capsule_of_the_liver", "VH_F_bare_area_of_liver"},
                  axes=(2, 0, 1), signs=(1, 1, 1), longest=None, clamp_frac=0.15),
    # ear: one closed mesh, arbitrary units. glTF y = height, z = front-back, x = depth.
    # sim (x,y,z) = (gx, -gz, gy); scaled so the height is a typical adult auricle (~6.5 cm).
    # temporal bone (NIH 3DPX-014739, base temporal bone): CT-derived, in millimetres, one mesh.
    # +x is the lateral cortex (mastoid, faces the surgeon), +z up: axes pass straight through, mm -> m.
    # The medial slab (smallest x) is clamped as the skull it is attached to.
    "bone": dict(file="35bda85f-basetemporalbone.glb", parts=None,
                 axes=(0, 1, 2), signs=(1, 1, 1), longest=None, scale=0.001, clamp_layer=0.2),
    "ear": dict(file="7c5d9c51-ear.glb", parts=None,
                axes=(0, 2, 1), signs=(1, -1, 1), longest=("z", 0.065), clamp_frac=0.25),
}

# Optional source-file overrides make the generator usable with locally named GLBs.
for _organ, _env in (
    ("liver", "ORGAN_LIVER_GLB"),
    ("ear", "ORGAN_EAR_GLB"),
    ("bone", "ORGAN_BONE_GLB"),
):
    PRESETS[_organ]["file"] = os.environ.get(_env, PRESETS[_organ]["file"])

def build(organ, h):
    cfg = PRESETS[organ]
    V, F = surface(SRC + cfg["file"], cfg["parts"])
    eu = edge_use(F)
    assert set(eu) <= {2, 4} and 1 not in eu, "surface has open edges: %s" % eu   # 4 = a few touching (non-manifold) edges, harmless here
    # rotate into the sim frame
    V = V[:, cfg["axes"]] * np.array(cfg["signs"], float)
    if np.array(cfg["signs"]).prod() < 0: F = F[:, [0, 2, 1]]        # keep triangle orientation consistent
    if "scale" in cfg: V = V * cfg["scale"]
    if cfg["longest"] is not None:
        axis = {"x": 0, "y": 1, "z": 2}[cfg["longest"][0]]
        s = cfg["longest"][1] / (V[:, axis].max() - V[:, axis].min())
        V = V * s
    c = 0.5 * (V.min(0) + V.max(0)); V = V - c
    surf_vol = signed_volume(V, F)
    if surf_vol < 0: F = F[:, [0, 2, 1]]; surf_vol = -surf_vol

    lo, hi = V.min(0), V.max(0)
    n = np.ceil((hi - lo) / h).astype(int) + 1
    # centre the grid on the organ
    origin = 0.5 * (lo + hi) - 0.5 * n * h
    idx = lambda i, j, k: (k * (n[1] + 1) + j) * (n[0] + 1) + i
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    tets, nv = [], (n[0] + 1) * (n[1] + 1) * (n[2] + 1)
    for k in range(n[2]):
        for j in range(n[1]):
            for i in range(n[0]):
                for a, b, c3 in perms:
                    p0 = np.array([i, j, k]); e = np.eye(3, dtype=int)
                    p1 = p0 + e[a]; p2 = p1 + e[b]; p3 = p2 + e[c3]
                    tets.append([idx(*p0), idx(*p1), idx(*p2), idx(*p3)])
    tets = np.array(tets)
    gx, gy, gz = np.meshgrid(np.arange(n[0] + 1), np.arange(n[1] + 1), np.arange(n[2] + 1), indexing="ij")
    grid = np.zeros((nv, 3))
    ii = (gz * (n[1] + 1) + gy) * (n[0] + 1) + gx
    grid[ii.ravel()] = origin + h * np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)

    ray = RayInside(V, F, h * 0.5); rng = np.random.default_rng(1)
    cent = grid[tets].mean(axis=1)
    keep = np.array([ray.inside(p, rng) for p in cent])
    tets = tets[keep]

    # largest face-connected component
    faces = defaultdict(list)
    for t, T in enumerate(tets):
        for o in range(4):
            faces[tuple(sorted(T[[k for k in range(4) if k != o]]))].append(t)
    parent = list(range(len(tets)))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for ts in faces.values():
        if len(ts) == 2: parent[find(ts[0])] = find(ts[1])
    roots = Counter(find(t) for t in range(len(tets)))
    big = roots.most_common(1)[0][0]
    tets = tets[[t for t in range(len(tets)) if find(t) == big]]

    used = np.unique(tets); remap = -np.ones(nv, int); remap[used] = np.arange(len(used))
    P = grid[used]; T = remap[tets]
    # Element-volume and shape diagnostics.  The Kuhn grid split should
    # produce non-degenerate, similarly shaped tetrahedra; fail loudly if not.
    tet_edges = np.stack([
        P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]],
        P[T[:, 3]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 1]],
        P[T[:, 3]] - P[T[:, 1]], P[T[:, 3]] - P[T[:, 2]]
    ], axis=1)
    det = np.linalg.det(np.stack([
        P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]],
        P[T[:, 3]] - P[T[:, 0]]
    ], axis=1))
    tet_volumes = np.abs(det) / 6.0
    if np.any(tet_volumes <= max(h ** 3 * 1e-12, 1e-18)):
        raise ValueError("%s mesh contains degenerate/near-zero-volume tetrahedra" % organ)
    edge_sq_sum = np.sum(tet_edges * tet_edges, axis=(1, 2))
    # Dimensionless mean-ratio quality: 1 for a regular tetrahedron, 0 for degenerate.
    quality = 12.0 * (3.0 * tet_volumes) ** (2.0 / 3.0) / edge_sq_sum
    vol = float(tet_volumes.sum())
    # Boundary condition: clamp the back support region (smallest x).  The presets
    # use clamp_frac for liver/ear and clamp_layer for bone; both are fractions
    # of the organ's x-depth.  Fall back to the surface-normal rule only for
    # presets that do not specify either setting.
    bfaces = defaultdict(list)
    for t, tt in enumerate(T):
        for o in range(4):
            f = [tt[k] for k in range(4) if k != o]
            bfaces[tuple(sorted(f))].append((f, tt[o]))
    fixed_set = set()
    if "clamp_frac" in cfg or "clamp_layer" in cfg:
        fraction = float(cfg.get("clamp_frac", cfg.get("clamp_layer")))
        if not (0.0 < fraction < 1.0):
            raise ValueError("%s clamp fraction must be between 0 and 1" % organ)
        xmin, xmax = P[:, 0].min(), P[:, 0].max()
        fixed_set = set(np.where(
            P[:, 0] <= xmin + fraction * (xmax - xmin) + 1e-12
        )[0].tolist())
    else:
        for key, lst in bfaces.items():
            if len(lst) != 1:
                continue
            f, opp = lst[0]
            nrm = np.cross(P[f[1]] - P[f[0]], P[f[2]] - P[f[0]])
            norm_len = np.linalg.norm(nrm)
            if norm_len < 1e-15:
                continue
            if np.dot(nrm, P[opp] - P[f[0]]) > 0:
                nrm = -nrm
            nrm = nrm / norm_len
            if nrm[0] < -0.5:
                fixed_set.update(f)
    fixed = np.array(sorted(fixed_set))
    return dict(P=P, T=T, fixed=fixed, vol=vol, surf_vol=surf_vol,
                quality_min=float(quality.min()), quality_mean=float(quality.mean()),
                quality_p05=float(np.percentile(quality, 5)),
                surfV=V, surfF=F, n=n, h=h)

def write_header(m, organ, path):
    P, T, fx = m["P"], m["T"], m["fixed"]
    L = []
    L.append("// generated by make_tet_mesh.py - do not edit.  organ: %s, cell size %.4g m" % (organ, m["h"]))
    L.append("// %d vertices, %d tetrahedra, %d clamped vertices, mesh volume %.1f cm^3 (surface volume %.1f cm^3)" %
             (len(P), len(T), len(fx), m["vol"] * 1e6, m["surf_vol"] * 1e6))
    L.append("#pragma once")
    L.append("namespace organ_mesh {")
    def arr(name, typ, data, fmt, per):
        L.append("static const %s %s[] = {" % (typ, name))
        flat = data.ravel()
        for s in range(0, len(flat), per):
            L.append("  " + ",".join(fmt % v for v in flat[s:s + per]) + ",")
        L.append("};")
    arr("%s_verts" % organ, "double", P, "%.6g", 9)
    arr("%s_tets" % organ, "int", T, "%d", 12)
    arr("%s_fixed" % organ, "int", fx, "%d", 16)
    L.append("static const int %s_numVerts = %d;" % (organ, len(P)))
    L.append("static const int %s_numTets = %d;" % (organ, len(T)))
    L.append("static const int %s_numFixed = %d;" % (organ, len(fx)))
    L.append("}")
    open(path, "w").write("\n".join(L) + "\n")

if __name__ == "__main__":
    organ, h, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    m = build(organ, h)
    write_header(m, organ, out)
    P = m["P"]
    ratio = m["vol"] / m["surf_vol"] if m["surf_vol"] > 0 else float("nan")
    print("%s: h=%.4g  verts=%d tets=%d fixed=%d  mesh vol=%.1f cm3 vs surface vol=%.1f cm3 (ratio=%.3f)  extents(cm)=%s" %
          (organ, h, len(P), len(m["T"]), len(m["fixed"]), m["vol"] * 1e6, m["surf_vol"] * 1e6,
           ratio, ((P.max(0) - P.min(0)) * 100).round(1)))
    print("  tet quality (mean-ratio): min=%.3f p05=%.3f mean=%.3f; clamped fraction=%.1f%%" %
          (m["quality_min"], m["quality_p05"], m["quality_mean"],
           100.0 * len(m["fixed"]) / max(1, len(P))))
    if not (0.95 <= ratio <= 1.05):
        print("  WARNING: mesh/surface volume mismatch exceeds 5%%; inspect source surface, units, and inside test.")
    np.save(out + ".npy", dict(P=P, T=m["T"], fixed=m["fixed"], V=m["surfV"], F=m["surfF"]), allow_pickle=True)
