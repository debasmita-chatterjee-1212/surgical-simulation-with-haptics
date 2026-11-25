44-soft-organ-scalpel: soft organs (liver, ear cartilage, left kidney, pancreas) with a probe and a scalpel.
Keys: [1] liver [2] ear [3] kidney [4] pancreas | mouse: move | arrows / scroll: depth | [t] probe / scalpel | hold SPACE while pressing in: cut
      [c] clear cuts | [ ] softer / stiffer | [r] reset | [e] edges | [f] full screen | [q] quit
On screen: force (N, with bar), deformation (max tissue displacement), organ and contact state, and for the scalpel the cut state, depth and total cut length.

USE: move the scalpel onto the organ, press in a little (the blade starts cutting at 0.15 N), hold SPACE and press in further: the cut depth follows how far you
push (at most 4 mm: shallow cuts) while you drag across the surface. Release SPACE or lift off to end the cut. Crossing an existing cut feels weaker.

INSTALL (PowerShell): create examples\GLFW\44-soft-organ-scalpel, copy every .h and the .cpp into it (not tools\), add 44-soft-organ-scalpel to the example list
in examples\GLFW\CMakeLists.txt (after 43-soft-organ-probe), run cmake .. and build target 44-soft-organ-scalpel.

WHAT THE SCALPEL MODELS (scalpel.h)
 - Surface-only cuts: the tetrahedral mesh is not changed. A cut is a list of points anchored to surface triangles (triangle + barycentric weights),
   so it moves and stretches with the tissue; it is drawn as a dark gape with a black slit that widens with the depth.
 - Cutting force F = J * d (fracture mechanics, steady cutting): J = work of fracture per unit crack area, d = depth of the cut.
   The tissue is loaded with this force, limited to the organ's own force limit.
 - Not modelled: the groove does not open under tension, no bleeding, no friction on the blade, no cutting below 4 mm.

VALUES AND SOURCES
 J = 200 J/m^2 (= N/m)   LIVER: scalpel cutting of porcine liver, average fracture resistance 187-225 J/m^2 (Chanthasopeephan, Desai & Lau, Drexel,
                         "Determining fracture characteristics in scalpel cutting of soft tissue" - CHECK journal and year in the paper itself).
                         Bovine liver, needle insertion: 164 +/- 6 J/m^2 (Gokgol, Basdogan & Canadinc 2012, Med Eng Phys 34:882-891).
                         KIDNEY, PANCREAS, EAR: ASSUMED (liver value). Ear cartilage is certainly tougher: needs its own measured value.
 ASSUMED (mine, not from papers): start-of-cut force 0.15 N, maximum cut depth 4 mm, softness factor 0.6 on an existing cut, depth smoothing 0.08 s.
 Organ force limits (liver, ear 3 N; kidney 1.2 N; pancreas 0.1 N) cap the force the tissue is loaded with, so the pancreas shows at most 0.1 N.
