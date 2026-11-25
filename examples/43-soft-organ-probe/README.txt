43-soft-organ-probe: soft organ probing - liver, ear cartilage, left kidney, pancreas.
Keys: [1] liver  [2] ear  [3] kidney  [4] pancreas.  Mouse: move probe.  Arrows/scroll: depth.  [ ] softer/stiffer.  [r] reset  [e] edges  [f] full screen  [q] quit.
On screen: force on the probe (N, with a bar up to the force limit), deformation (maximum tissue displacement, mm), organ and contact state.

INSTALL: copy every file in this folder (not tools\) into
   examples\GLFW\43-soft-organ-probe\   (overwrite), then
   Get-ChildItem <that folder> | ForEach-Object { $_.LastWriteTime = Get-Date }
   cmake --build . --config Release --target 43-soft-organ-probe
tools\ holds the scripts that made the meshes; keep it OUTSIDE the build folder.

SOURCES
 liver     BodyParts3D FJ192   E = 15 kPa (Ozcan et al. 2011, Med Image Anal 15:45-52)
 kidney    BodyParts3D FJ456   E = 19.6 kPa (Snedeker et al. 2005, J Biomech 38:993-1001, as quoted in PMC7366619; check in the primary paper). nu 0.45, density 1050 ASSUMED.
 pancreas  BodyParts3D FJ271   E = 1.4 kPa, nu 0.46 (Pagnanelli et al. 2024, Biomedicines 12(1):91). Density 1050 ASSUMED.
 ear       NIH 3D model        E = 1.66 MPa (Griffin et al. 2016, Ann Biomed Eng 44(12):3460-3467)
 BodyParts3D: (c) The Database Center for Life Science, CC Attribution 4.0 International.

SETTINGS THAT ARE MINE (not from papers): hand-proxy spring 800 N/m, force limits (liver, ear 3 N; kidney 1.2 N; pancreas 0.1 N),
 maximum hand depth below the intact surface (kidney 11 mm, pancreas 10 mm). The limits keep the linear tensor-mass model valid:
 beyond ~50 % local strain it collapses the elements under a small tip.
