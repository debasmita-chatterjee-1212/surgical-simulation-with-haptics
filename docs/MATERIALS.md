# Material parameters

"Cited" = a published source gives the value. "Checked" = the paper was opened when this table was written (November 2026).
Values still without a source: all densities, the ear Poisson ratio, the bone Poisson ratio and density.

| Organ | Parameter | Value | Source and status |
| --- | --- | --- | --- |
| Liver | E | 15 kPa | Ozcan et al. 2011, Med Image Anal 15:45-52 (cited; not re-opened) |
| Liver | nu | 0.45 | Low end of the 0.45-0.495 range used for soft tissue (Fougeron et al. 2023). Numerical choice: linear tetrahedra lock as nu approaches 0.5 |
| Liver | density | 1060 kg/m^3 | Assumed within realistic values |
| Kidney | E | 19.6 kPa | Human cortex, low-strain compression, Snedeker et al. 2005 (J Biomech 38:993-1001), as quoted by Hattab et al. 2020 (Sci Rep 10:11849). The Snedeker abstract gives only failure data: confirm in the full text. Alternative: whole human kidney 5.6 kPa at 1 %/s (Johnson et al. 2021, Biomechanics 1(2):264-280) |
| Kidney | nu, density | 0.45, 1050 kg/m^3 | nu as liver; density assumed within realistic values |
| Pancreas | E, nu | 1.4 kPa, 0.46 | Pagnanelli et al. 2024, Biomedicines 12(1):91 (cited; not re-opened) |
| Pancreas | density | 1050 kg/m^3 | Assumed within realistic values |
| Ear cartilage | E | 1.66 MPa | Griffin et al. 2016, Ann Biomed Eng 44(12):3460-3467: human auricular cartilage, indentation, mean 1.66 +/- 0.63 MPa. Checked |
| Ear cartilage | nu, density | 0.40, 1100 kg/m^3 | Assumed; Griffin et al. report |
| Temporal bone | E | 10 GPa | Inside the 7.5-15.5 GPa measured in parietal and frontal bone (Motherway et al. 2009, J Biomech 42(13)). No temporal-bone value found |
| Temporal bone | nu, density | 0.3, 1900 kg/m^3 | Common modelling values, not measured on temporal bone. Only E and nu enter the bone readout |
| Scalpel, liver | J | 200 J/m^2 | Porcine liver, scalpel: 187-225 J/m^2 (Chanthasopeephan, Desai and Lau, Drexel; confirm journal and year). Bovine liver, needle: 164 +/- 6 J/m^2 (Gokgol et al. 2012, Med Eng Phys 34:882-891) |
| Scalpel, other organs | J | 200 J/m^2 | Assumed (liver value) within realistic values. No measured scalpel toughness for ear cartilage found |

Settings chosen by the author from experimental observations: coupling spring 800 N/m; force limits (liver and ear 3 N, kidney 1.2 N, pancreas 0.1 N);
hand-depth limits (kidney 11 mm, pancreas 10 mm); cut start force 0.15 N; maximum cut depth 4 mm; softness factor 0.6 on an existing cut.
