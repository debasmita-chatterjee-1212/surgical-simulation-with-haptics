//==============================================================================
// organs.h  -  organ presets (mesh + material + probe) for the surgical simulator.
//
// Material sources (cite these in the thesis):
//   liver : E = 15 kPa  - human liver, ex vivo, storage modulus ~10-20 kPa without fibrosis
//           Ozcan MU, Ocal S, Basdogan C, Dogusoy G, Tokat Y (2011) "Characterization of
//           frequency-dependent material properties of human liver and its pathologies using
//           an impact hammer", Medical Image Analysis 15:45-52.
//   ear   : E = 1.66 MPa - human auricular cartilage, indentation, mean of 14 sites
//           Griffin MF, Premakumar Y, Seifalian AM, Szarko M, Butler PEM (2016) "Biomechanical
//           characterisation of the human auricular cartilages; implications for tissue
//           engineering", Annals of Biomedical Engineering 44(12):3460-3467.
//   ASSUMED (not from these papers): Poisson ratios (0.45 liver, 0.40 cartilage), densities
//           (1060 and 1100 kg/m^3). Replace by values you can cite.
// Meshes: liver = BodyParts3D FJ192 (c) The Database Center for Life Science, CC BY 4.0;
//         ear   = NIH 3D Print Exchange model; both converted by make_tet_mesh.py.
// kidney, pancreas = BodyParts3D FJ456 / FJ271 (c) The Database Center for Life Science, CC BY 4.0.
// The soft gel block and the temporal bone are separate programs (41-gel-probe, 42-bone-probe).
//==============================================================================
#pragma once

#include "surgical_sim.h"
#include "mesh_liver.h"
#include "mesh_ear.h"
#include "mesh_kidney.h"
#include "mesh_pancreas.h"

namespace tmass {

enum OrganId { ORGAN_BLOCK = 0, ORGAN_LIVER = 1, ORGAN_EAR = 2, ORGAN_BONE = 3, ORGAN_KIDNEY = 4, ORGAN_PANCREAS = 5, ORGAN_COUNT = 6 };

inline const char* organName(int id)
{
    switch (id)
    {
    case ORGAN_LIVER: return "liver (E = 15 kPa)";
    case ORGAN_EAR:   return "ear cartilage (E = 1.66 MPa)";
    case ORGAN_BONE:  return "temporal bone (E = 10 GPa, PLACEHOLDER value)";
    case ORGAN_KIDNEY:   return "left kidney (E = 19.6 kPa)";
    case ORGAN_PANCREAS: return "pancreas (E = 1.4 kPa)";
    default:          return "soft gel block (E = 7 kPa)";
    }
}

inline SimConfig organPreset(int id)
{
    SimConfig c;                                   // default = soft block
    c.probeShaftLength = 0.03; c.maxIndentation = 0.010;   // gel: shaft must be several times longer than the indentation
    if (id == ORGAN_LIVER)
    {
        c.meshVerts = organ_mesh::liver_verts; c.meshTets = organ_mesh::liver_tets; c.meshFixed = organ_mesh::liver_fixed;
        c.meshNumVerts = organ_mesh::liver_numVerts; c.meshNumTets = organ_mesh::liver_numTets;
        c.meshNumFixed = organ_mesh::liver_numFixed;
        c.youngModulus = 15.0e3; c.poisson = 0.45; c.density = 1060.0;
        c.dampingPerMass = 30.0;
        c.probeRadius = 0.00125; c.drillRadius = 0.014; // small contact tip; cutting radius remains independent
        c.contactStiffness = 400.0; c.contactMaxForce = 3.0;
        c.contactStiffnessPerArea = 20.0 * c.youngModulus / c.probeRadius;   // penalty penetration ~ few % of the indentation
        c.probeMaxSpeed = 0.6;
        c.massScale = 4.0;
        c.probeShaftLength = 0.05; c.maxIndentation = 0.020;
        c.proxyContact = true; c.couplingStiffness = 800.0;   // proxy (god-object) contact, see surgical_sim.h
        c.cuttingEnergyPerVolume = 1.0e4; // prototype calibration [J/m^3]
    }
    else if (id == ORGAN_EAR)
    {
        c.meshVerts = organ_mesh::ear_verts; c.meshTets = organ_mesh::ear_tets; c.meshFixed = organ_mesh::ear_fixed;
        c.meshNumVerts = organ_mesh::ear_numVerts; c.meshNumTets = organ_mesh::ear_numTets;
        c.meshNumFixed = organ_mesh::ear_numFixed;
        c.youngModulus = 1.66e6; c.poisson = 0.40; c.density = 1100.0;
        c.dampingPerMass = 80.0;
        c.probeRadius = 0.00125; c.drillRadius = 0.005; // small contact tip; cutting radius remains independent
        c.contactStiffness = 4000.0; c.contactMaxForce = 3.0;
        c.contactStiffnessPerArea = 15.0 * c.youngModulus / c.probeRadius;
        c.probeMaxSpeed = 0.2;
        c.massScale = 2000.0;
        c.probeShaftLength = 0.03; c.maxIndentation = 0.005;
        c.proxyContact = true; c.couplingStiffness = 800.0;
        c.cuttingEnergyPerVolume = 2.0e6; // prototype calibration [J/m^3]
    }
    else if (id == ORGAN_KIDNEY)
    {
        c.meshVerts = organ_mesh::kidney_verts; c.meshTets = organ_mesh::kidney_tets; c.meshFixed = organ_mesh::kidney_fixed;
        c.meshNumVerts = organ_mesh::kidney_numVerts; c.meshNumTets = organ_mesh::kidney_numTets;
        c.meshNumFixed = organ_mesh::kidney_numFixed;
        c.youngModulus = 19.6e3; c.poisson = 0.45; c.density = 1050.0;   // E: human renal cortex, low strain (Snedeker et al. 2005, as quoted in PMC7366619); nu, density ASSUMED
        c.dampingPerMass = 30.0;
        c.probeRadius = 0.00125; c.drillRadius = 0.008;
        c.contactStiffness = 400.0;
        c.contactMaxForce = 1.2;   // [N] force limit: above this the linear model collapses the elements under the tip (volume ratio < 0.3)
        c.probeMaxSpeed = 0.6;
        c.massScale = 4.0;
        c.proxyContact = true; c.couplingStiffness = 800.0; c.loadRadiusEdges = 2.4; c.proxyMaxHandDepth = 0.011;
        c.cuttingEnergyPerVolume = 1.0e4; // prototype calibration [J/m^3], not validated
    }
    else if (id == ORGAN_PANCREAS)
    {
        c.meshVerts = organ_mesh::pancreas_verts; c.meshTets = organ_mesh::pancreas_tets; c.meshFixed = organ_mesh::pancreas_fixed;
        c.meshNumVerts = organ_mesh::pancreas_numVerts; c.meshNumTets = organ_mesh::pancreas_numTets;
        c.meshNumFixed = organ_mesh::pancreas_numFixed;
        c.youngModulus = 1.4e3; c.poisson = 0.46; c.density = 1050.0;    // E, nu: Pagnanelli et al. 2024 (indentation of human pancreas; nu = 0.46 was their assumption); density ASSUMED
        c.dampingPerMass = 20.0;
        c.probeRadius = 0.00125; c.drillRadius = 0.006;
        c.contactStiffness = 400.0;
        c.contactMaxForce = 0.1;   // [N] force limit: E = 1.4 kPa is ~14x softer than the kidney, so the same dimple needs ~14x less force
        c.probeMaxSpeed = 0.6;
        c.massScale = 1.0;
        c.proxyContact = true; c.couplingStiffness = 800.0; c.loadRadiusEdges = 2.4; c.proxyMaxHandDepth = 0.010;
        c.cuttingEnergyPerVolume = 5.0e3; // prototype calibration [J/m^3], not validated
    }
    return c;
}

// how the window should frame and drive each organ
struct OrganView
{
    double keySpeed;       // [m/s] UP/DOWN arrow speed (probe depth)
    double scrollStep;     // [m] per scroll tick
    double cameraDist;     // [m]
    double workspaceHalfY, workspaceHalfZ;   // [m] range the mouse covers
};

inline OrganView organView(int id)
{
    OrganView v;
    if (id == ORGAN_LIVER)    { v.keySpeed = 0.06;  v.scrollStep = 0.004;  v.cameraDist = 0.55;  v.workspaceHalfY = 0.20;  v.workspaceHalfZ = 0.15; }
    else if (id == ORGAN_EAR) { v.keySpeed = 0.012; v.scrollStep = 0.0008; v.cameraDist = 0.14;  v.workspaceHalfY = 0.05;  v.workspaceHalfZ = 0.04; }
    else if (id == ORGAN_KIDNEY)   { v.keySpeed = 0.05;  v.scrollStep = 0.003;  v.cameraDist = 0.30;  v.workspaceHalfY = 0.09;  v.workspaceHalfZ = 0.09; }
    else if (id == ORGAN_PANCREAS) { v.keySpeed = 0.05;  v.scrollStep = 0.003;  v.cameraDist = 0.32;  v.workspaceHalfY = 0.10;  v.workspaceHalfZ = 0.07; }
    else if (id == ORGAN_BONE){ v.keySpeed = 0.015; v.scrollStep = 0.001;  v.cameraDist = 0.30;  v.workspaceHalfY = 0.09;  v.workspaceHalfZ = 0.07; }
    else                      { v.keySpeed = 0.08;  v.scrollStep = 0.004;  v.cameraDist = 0.30;  v.workspaceHalfY = 0.13;  v.workspaceHalfZ = 0.10; }
    return v;
}

} // namespace tmass
