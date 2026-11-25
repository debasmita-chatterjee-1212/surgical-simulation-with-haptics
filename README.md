# Haptic surgical simulation in CHAI3D

Touch and cut virtual tissue with a mouse or touchpad or the Phantom Omni haptic device while force and deformation are shown on screen.
Built on [CHAI3D](https://www.chai3d.org/) for Windows (Visual Studio, CMake, GLFW examples).

| Example | What it does |
| --- | --- |
| `examples/42-bone-probe` | Rigid left temporal bone (BodyParts3D, 0.5 mm voxels); probe with continuous force and Hertz indentation readout |
| `examples/43-soft-organ-probe` | Liver, ear cartilage, left kidney, pancreas as deformable tetrahedral meshes; probe with force and maximum tissue displacement |
| `examples/44-soft-organ-scalpel` | The same four organs with a probe and a scalpel (shallow surface cuts up to 4 mm, cutting force F = J d) |
| `validation/` | Headless checks of the solver and contact model against closed-form solutions, plus the script that draws the error figures |
| `tools/` | Python scripts that build the organ meshes from BodyParts3D files |

## Install and run (Windows)
1. Get CHAI3D and build it once (Visual Studio, CMake).
2. From a PowerShell in this repository: `.\scripts\install.ps1 -Chai3dRoot "<path to chai3d>"` (copies the examples and registers them in `examples\GLFW\CMakeLists.txt`).
3. `cd <chai3d>\build; cmake ..; cmake --build . --config Release --target 44-soft-organ-scalpel`
4. Run `<chai3d>\bin\win-x64\44-soft-organ-scalpel.exe`.

Keys (43 and 44): `1-4` organ, mouse moves the instrument, arrows or scroll set depth, `[` `]` softer or stiffer tissue, `r` reset, `e` edges, `f` full screen, `q` quit.
In 44: `t` switches probe and scalpel, hold `SPACE` while pressing in to cut, `c` clears cuts.

## Validation
`validation/validation.cpp` (build: `cl /O2 /EHsc /std:c++17 validation.cpp`) runs four tests; results are in `validation/example_output_from_my_run.txt` and `validation/validation_results.csv`.
Summary at the finest mesh: bar in compression -1.2 % (nu 0.30) and -3.6 % (nu 0.45); cantilever -7.8 % and -18.2 %; half-space under the proxy load patch -11.9 % and -18.0 %.
The model is stiffer than theory on coarse meshes (locking of linear tetrahedra) and converges with refinement. The force readout of the contact loop is qualitative.
Figures: `python validation/figures/make_plots.py validation/figures/validation_results.csv` (needs Python and matplotlib).

## Material values
See [docs/MATERIALS.md](docs/MATERIALS.md). Each value is marked cited, checked or assumed; densities, the ear Poisson ratio and the bone Poisson ratio and density are assumed.

## Data and licences
- Liver, kidney, pancreas and temporal bone geometry: [BodyParts3D](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/desc.html), (c) The Database Center for Life Science, CC Attribution 4.0 International. The mesh headers and voxel files in this repository are derived from it.
- Ear geometry: NIH 3D model; check its licence before making this repository public.
- CHAI3D is not included and has its own licence.
- This repository has no licence file yet: add one before sharing the code.
