validation.cpp - checks the tensor-mass solver and the proxy contact against closed-form solutions. No graphics, no CHAI3D.

KEEP THIS FOLDER OUTSIDE THE CHAI3D build tree (for example E:\thesis-sim\validation\).

BUILD AND RUN (Developer PowerShell):
   cd E:\thesis-sim\validation
   cl /O2 /EHsc /std:c++17 validation.cpp /Fe:validation.exe
   .\validation.exe          (all tests, a few minutes)     or   .\validation.exe 3   (one test: 1, 2, 3 or 4)
Output: a table on screen and validation_results.csv.

TESTS (E = 15 kPa, liver value, Ozcan et al. 2011; nu = 0.30 and 0.45; several cell sizes)
 1 uniaxial compression of a clamped bar       delta = F L / (A E)
 2 cantilever, tip load                        Timoshenko beam: F L^3/(3 E I) + F L/(kappa G A), kappa = 5/6
 3 half-space loaded by the proxy load patch   w0 = 4 F (1 - nu^2) / (pi E rho)  (pressure (1 - r/rho)^2; derivation in the code)
 4 full proxy contact loop on the same block   reports force, proxy dent and tissue displacement, flags unsteady holds

example_output_from_my_run.txt is the result on my machine; your numbers should match.
