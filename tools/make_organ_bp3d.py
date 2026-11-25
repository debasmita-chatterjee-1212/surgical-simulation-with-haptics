"""Tetrahedral meshes of BodyParts3D organs (CC BY 4.0) with make_tet_mesh.py's mesher (mtm_base.py).
usage: python make_organ_bp3d.py  kidney|pancreas  cell_mm  out.h
BodyParts3D axes: x = patient left, y = posterior, z = superior, mm.  Simulation frame: +x faces the probe
(anterior), z up:  sim = (-y, x, z).  The posterior layer (smallest x) is clamped (organ attached behind it)."""
import sys, collections
import numpy as np
from scipy import ndimage as ndi
from skimage import measure
import mtm_base as mtm

FILES = {
  "kidney":   "/tmp/organs/kidney/BP1494_FMA3_0_partof_FMA7205_Left_kidney/FJ456_BP1494_FMA7205_Left kidney.obj",
  "pancreas": "/tmp/organs/pancreas/BP1548_FMA3_0_partof_FMA7198_Pancreas/FJ271_BP22_FMA7198_Pancreas.obj",
}

def read_obj(path):
    V = []; F = []
    for l in open(path):
        if l.startswith("v "): V.append(list(map(float, l.split()[1:4])))
        elif l.startswith("f "): F.append([int(t.split("/")[0]) - 1 for t in l.split()[1:]])
    return np.array(V), np.array(F)

def cap_holes(V, F):
    D = collections.Counter()
    for a, b, c in F:
        for u, v in ((a, b), (b, c), (c, a)): D[(u, v)] += 1
    bd = [(u, v) for (u, v) in D if (v, u) not in D]
    nxt = collections.defaultdict(list)
    for u, v in bd: nxt[u].append(v)
    seen = set(); loops = []
    for u, v in bd:
        if (u, v) in seen: continue
        loop = [u]; seen.add((u, v)); nx = v
        while nx != u:
            loop.append(nx)
            cand = [w for w in nxt[nx] if (nx, w) not in seen]
            if not cand: raise RuntimeError("open boundary does not close into a loop")
            seen.add((nx, cand[0])); nx = cand[0]
        loops.append(loop)
    newV = [v for v in V]; newF = [list(f) for f in F]
    for loop in loops:
        c = V[loop].mean(0); ci = len(newV); newV.append(c)
        for i in range(len(loop)):
            a, b = loop[i], loop[(i + 1) % len(loop)]
            newF.append([b, a, ci])
    print("closed %d holes (%d boundary edges)" % (len(loops), len(bd)))
    return np.array(newV), np.array(newF)

def voxel_close(V, F, h=0.5, sigma=0.9):
    """fragmented / open surface -> closed smooth surface: rasterise, seal, flood-fill, largest piece, marching cubes"""
    lo = V.min(0) - 4 * h; hi = V.max(0) + 4 * h
    n = np.ceil((hi - lo) / h).astype(int)
    surf = np.zeros(n, bool)
    for a, b, c in F:
        p0, p1, p2 = V[a], V[b], V[c]
        L = max(np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p1), np.linalg.norm(p0 - p2))
        k = int(np.ceil(L / (h * 0.4))) + 1
        u, v = np.meshgrid(np.linspace(0, 1, k), np.linspace(0, 1, k)); m = u + v <= 1
        P = p0 + np.outer(u[m], p1 - p0) + np.outer(v[m], p2 - p0)
        I = np.floor((P - lo) / h).astype(int); surf[I[:, 0], I[:, 1], I[:, 2]] = True
    sealed = ndi.binary_closing(surf, iterations=2)
    lab = ndi.label(~sealed)[0]
    solid = lab != lab[0, 0, 0]
    solid = ndi.binary_erosion(solid, iterations=1)           # undo the one-voxel skin added by rasterising + closing
    lab2, nl = ndi.label(solid)
    sizes = ndi.sum(solid, lab2, range(1, nl + 1))
    solid = lab2 == (1 + int(np.argmax(sizes)))
    print("voxel solid %.1f cm3 in %d piece(s), kept the largest (%.1f cm3)" % (ndi.sum(lab != lab[0,0,0]) * h**3 / 1000, nl, solid.sum() * h**3 / 1000))
    vol = ndi.gaussian_filter(np.pad(solid.astype(np.float32), 2), sigma)
    verts, faces, _, _ = measure.marching_cubes(vol, 0.5, spacing=(h, h, h))
    verts = verts - 2 * h + lo + 0.5 * h * 0                       # back to file coordinates [mm]
    return verts, faces[:, [0, 2, 1]] if False else faces

def load(organ):
    V, F = read_obj(FILES[organ])
    eu = collections.Counter()
    for a, b, c in F:
        for u, v in ((a, b), (b, c), (c, a)): eu[(min(u, v), max(u, v))] += 1
    open_edges = sum(1 for n in eu.values() if n == 1)
    if organ == "kidney" and open_edges < 50:
        V, F = cap_holes(V, F)
    else:
        V, F = voxel_close(V, F)
    return V * 0.001, F

def surface(path, names=None): return load(path)
mtm.surface = surface; mtm.SRC = ""
for o in FILES:
    mtm.PRESETS[o] = dict(file=o, parts=None, axes=(1, 0, 2), signs=(-1, 1, 1), longest=None)
mtm.PRESETS["kidney"]["clamp_frac"] = 0.15        # kidney: posterior layer clamped (retroperitoneal)
# pancreas: no clamp_frac -> the whole posterior-facing surface (normal.x < -0.5) is clamped, because the pancreas
# lies against the posterior abdominal wall along its whole length, not only at one end

if __name__ == "__main__":
    organ, h, out = sys.argv[1], float(sys.argv[2]) * 1e-3, sys.argv[3]
    m = mtm.build(organ, h)
    mtm.write_header(m, organ, out)
    P = m["P"]
    print("%s: h=%.2f mm verts=%d tets=%d fixed=%d mesh vol=%.1f cm3 vs surface %.1f cm3 (ratio %.3f) extents(mm)=%s" %
          (organ, h * 1e3, len(P), len(m["T"]), len(m["fixed"]), m["vol"] * 1e6, m["surf_vol"] * 1e6, m["vol"] / m["surf_vol"], ((P.max(0) - P.min(0)) * 1e3).round(1)))
    print("tet quality min=%.3f p05=%.3f mean=%.3f" % (m["quality_min"], m["quality_p05"], m["quality_mean"]))
    np.save(out + ".npy", dict(P=P, T=m["T"], fixed=m["fixed"], V=m["surfV"], F=m["surfF"]), allow_pickle=True)
