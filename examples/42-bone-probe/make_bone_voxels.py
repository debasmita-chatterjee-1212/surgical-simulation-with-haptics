#!/usr/bin/env python3
"""
make_bone_voxels.py  -  turn a bone surface (OBJ, millimetres) into the two files the bone-probe demo reads:
    <prefix>_voxels.bin    solid/air grid used for the contact (uint8, 255 = bone)
    <prefix>_surface.bin   smooth triangle surface for display (metres)

usage:  python make_bone_voxels.py  input.obj  [voxel_mm=0.5]  [prefix=temporal_bone]

Needs numpy, scipy, scikit-image.

How the surface is closed: the BodyParts3D meshes are open and in many pieces, so the triangles are
rasterised into the grid, the air outside is flood-filled from the grid corner, and everything that is not
outside counts as bone. Gaps smaller than a voxel are sealed this way.
Axes: the grid is stored as (x, y, z) = (file y, file z, file x), i.e. the file's +x side (the lateral side of a
LEFT temporal bone in BodyParts3D) ends up on top (+z), which is the side the probe comes from.
"""
import sys, struct
import numpy as np
from scipy import ndimage as ndi
from skimage import measure

src = sys.argv[1]
h = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5          # mm
prefix = sys.argv[3] if len(sys.argv) > 3 else "temporal_bone"

V, F = [], []
for line in open(src):
    if line.startswith("v "):
        V.append([float(t) for t in line.split()[1:4]])
    elif line.startswith("f "):
        F.append([int(t.split("/")[0]) - 1 for t in line.split()[1:]])
V = np.array(V); F = np.array(F)
lo = V.min(0) - 3 * h; hi = V.max(0) + 3 * h
n = np.ceil((hi - lo) / h).astype(int)
surf = np.zeros(n, bool)
for a, b, c in F:
    p0, p1, p2 = V[a], V[b], V[c]
    L = max(np.linalg.norm(p1 - p0), np.linalg.norm(p2 - p1), np.linalg.norm(p0 - p2))
    k = int(np.ceil(L / (h * 0.4))) + 1
    u, v = np.meshgrid(np.linspace(0, 1, k), np.linspace(0, 1, k)); m = u + v <= 1
    P = p0 + np.outer(u[m], p1 - p0) + np.outer(v[m], p2 - p0)
    I = np.floor((P - lo) / h).astype(int)
    surf[I[:, 0], I[:, 1], I[:, 2]] = True
lab = ndi.label(~surf)[0]
solid = lab != lab[0, 0, 0]
solid = ndi.binary_erosion(solid)            # removes the one-voxel skin added by rasterising the surface
print("solid volume %.2f cm3 (file header says about the same)" % (solid.sum() * h ** 3 / 1000))

g = np.transpose(solid, (1, 2, 0))           # (x, y, z) = (file y, file z, file x)
nx, ny, nz = g.shape
with open(prefix + "_voxels.bin", "wb") as fh:
    fh.write(struct.pack("<iiid", nx, ny, nz, h))
    fh.write((np.ascontiguousarray(g.transpose(2, 1, 0)).astype(np.uint8) * 255).tobytes())   # x fastest

# display surface: marching cubes on a lightly blurred copy, level 0.5
pad = 3
vol = ndi.gaussian_filter(np.pad(g.astype(np.float32), pad), 0.8)
verts, faces, _, _ = measure.marching_cubes(vol, 0.5, spacing=(h, h, h))
verts = verts - pad * h + 0.5 * h                         # back to the grid frame (voxel centres at (i+0.5)h)
verts = verts * 1e-3                                      # metres
with open(prefix + "_surface.bin", "wb") as fh:
    fh.write(struct.pack("<ii", len(verts), len(faces)))
    fh.write(verts.astype("<f4").tobytes())
    fh.write(faces.astype("<u4").tobytes())
print("grid %d x %d x %d, voxel %.2f mm, surface %d vertices %d triangles" % (nx, ny, nz, h, len(verts), len(faces)))
print("extent [mm] x %.1f  y %.1f  z %.1f" % (nx * h, ny * h, nz * h))
