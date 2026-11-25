42-bone-probe : rigid temporal bone, touch probe (force + deformation readout)

FILES (all five go in  E:\chai3d-master\chai3d-master\examples\GLFW\42-bone-probe\ ):
  42-bone-probe.cpp, voxel_bone.h, temporal_bone_voxels.bin, temporal_bone_surface.bin, make_bone_voxels.py

BUILD (same recipe as 41-gel-probe):
  1. Put the folder in examples\GLFW\ and touch the files:
       Get-ChildItem E:\chai3d-master\chai3d-master\examples\GLFW\42-bone-probe | ForEach-Object { $_.LastWriteTime = Get-Date }
  2. Add  42-bone-probe  to the foreach list on line 44 of examples\GLFW\CMakeLists.txt, then in build\ run:  cmake ..
  3. cmake --build . --config Release --target 42-bone-probe
  4. Run  E:\chai3d-master\chai3d-master\bin\win-x64\42-bone-probe.exe
     (the exe finds the .bin files next to the source; or pass the folder as the first argument)

OUTPUT: bone_probe_log.csv in the folder you run from (time, force, spring force, hertz indentation, ...).

NEW BONE LATER (e.g. OpenEar):  python make_bone_voxels.py  bone.obj  0.5  temporal_bone
Bone: BodyParts3D FJ580 (c) The Database Center for Life Science, CC Attribution 4.0 International.
