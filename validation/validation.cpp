// validation.cpp - headless verification of the tensor-mass solver and of the proxy contact against
// closed-form (analytical) solutions. No graphics, no CHAI3D needed.
//
//  Test 1  uniaxial compression of a clamped bar        analytical: delta = F L / (A E)
//  Test 2  cantilever beam, tip load                    analytical: Timoshenko beam, delta = F L^3/(3 E I) + F L/(kappa G A)
//  Test 3  half-space loaded by the proxy's load patch   analytical: centre displacement of a half-space (derivation at test3)
//  Test 4  full proxy contact loop on the same block    force-indentation consistency (no closed form for the loop itself)
// Each test is repeated on several mesh sizes to show convergence.
//
// Build (Developer PowerShell):   cl /O2 /EHsc /std:c++17 validation.cpp /Fe:validation.exe
// Run:                            .\validation.exe        (writes validation_results.csv next to it)
#include "surgical_sim.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <utility>
#include <string>
using namespace tmass;

static FILE* csv = nullptr;
static const double PI = 3.14159265358979323846;

// ---- relax a loaded body to static equilibrium (damped explicit dynamics) ----------------------------------------
static void relax(TensorMassDynamics& dyn, const std::vector<int>& probeNodes, int axis)
{
    auto meanDisp = [&]() {
        double s = 0.0;
        for (int v : probeNodes) { Vec3 u = dyn.positions()[(size_t)v] - dyn.model().restPositions()[(size_t)v]; s += (axis == 0 ? u.x : axis == 1 ? u.y : u.z); }
        return s / (double)probeNodes.size();
    };
    double prev = 0.0;
    for (int it = 0; it < 4000; ++it)
    {
        for (int k = 0; k < 100; ++k) dyn.step();
        const double m = meanDisp();
        if (it > 5 && std::fabs(m - prev) < 1e-6 * std::fabs(m)) break;
        prev = m;
    }
}

// consistent nodal load weights of the x = max face of a regular nx*ny*nz grid (interior 1, edge 1/2, corner 1/4)
static void endFace(int nx, int ny, int nz, std::vector<int>& nodes, std::vector<double>& w)
{
    double sum = 0.0;
    for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j)
    {
        nodes.push_back((k * (ny + 1) + j) * (nx + 1) + nx);
        const double wy = (j == 0 || j == ny) ? 0.5 : 1.0, wz = (k == 0 || k == nz) ? 0.5 : 1.0;
        w.push_back(wy * wz); sum += wy * wz;
    }
    for (double& x : w) x /= sum;
}

static TensorMassDynamics makeBar(int nx, int ny, int nz, double h, double E, double nu, double damping)
{
    std::vector<Vec3> V; std::vector<TensorMassModel::Tet> T;
    makeBlockMesh(nx, ny, nz, h, V, T);
    TensorMassModel m; m.init(V, T, lameLambda(E, nu), lameMu(E, nu));
    TensorMassDynamics d; d.init(m, 1060.0, damping, 1e-4);
    for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) d.fixVertex((k * (ny + 1) + j) * (nx + 1) + 0);   // clamp x = 0
    d.setTimeStep(0.4 * d.estimateCriticalTimeStep());
    return d;
}

// ---------------------------------------------------------------------------------------------------------------------
static void test1()
{
    printf("\nTEST 1  uniaxial compression of a clamped bar (L = 60 mm, 20 x 20 mm, E = 15 kPa)\n");
    printf("  nu    cell h [mm]   tets     FEM [mm]   analytic [mm]   error [%%]\n");
    const double E = 15e3, L = 0.060, W = 0.020, F = 0.10;
    for (double nu : {0.30, 0.45})
        for (double h : {0.010, 0.005, 0.004, 0.0025})
        {
            const int nx = (int)std::lround(L / h), ny = (int)std::lround(W / h);
            TensorMassDynamics d = makeBar(nx, ny, ny, h, E, nu, 60.0);
            std::vector<int> nodes; std::vector<double> w; endFace(nx, ny, ny, nodes, w);
            for (size_t i = 0; i < nodes.size(); ++i) d.setExternalForce(nodes[i], Vec3(-F * w[i], 0, 0));
            relax(d, nodes, 0);
            double u = 0.0; for (int v : nodes) u += (d.positions()[(size_t)v].x - d.model().restPositions()[(size_t)v].x);
            u = -u / (double)nodes.size();
            const double an = F * L / (W * W * E);
            printf("  %.2f  %8.1f  %8zu  %9.4f  %12.4f  %10.2f\n", nu, h * 1e3, d.model().numTets(), u * 1e3, an * 1e3, 100.0 * (u - an) / an);
            if (csv) fprintf(csv, "uniaxial,nu=%.2f,%g,%zu,%g,%g,%g\n", nu, h * 1e3, d.model().numTets(), u * 1e3, an * 1e3, 100.0 * (u - an) / an);
        }
}

static void test2()
{
    printf("\nTEST 2  cantilever beam, tip load (L = 100 mm, 20 x 20 mm, E = 15 kPa), Timoshenko beam theory\n");
    printf("  nu    cell h [mm]   tets     FEM [mm]   analytic [mm]   error [%%]\n");
    const double E = 15e3, L = 0.100, W = 0.020, F = 1.2e-3;
    for (double nu : {0.30, 0.45})
        for (double h : {0.010, 0.005, 0.004, 0.0025})
        {
            const int nx = (int)std::lround(L / h), ny = (int)std::lround(W / h);
            TensorMassDynamics d = makeBar(nx, ny, ny, h, E, nu, 25.0);
            std::vector<int> nodes; std::vector<double> w; endFace(nx, ny, ny, nodes, w);
            for (size_t i = 0; i < nodes.size(); ++i) d.setExternalForce(nodes[i], Vec3(0, -F * w[i], 0));
            relax(d, nodes, 1);
            double u = 0.0; for (int v : nodes) u += (d.positions()[(size_t)v].y - d.model().restPositions()[(size_t)v].y);
            u = -u / (double)nodes.size();
            const double I = W * W * W * W / 12.0, A = W * W, G = E / (2.0 * (1.0 + nu));
            const double an = F * L * L * L / (3.0 * E * I) + F * L / ((5.0 / 6.0) * G * A);
            printf("  %.2f  %8.1f  %8zu  %9.4f  %12.4f  %10.2f\n", nu, h * 1e3, d.model().numTets(), u * 1e3, an * 1e3, 100.0 * (u - an) / an);
            if (csv) fprintf(csv, "cantilever,nu=%.2f,%g,%zu,%g,%g,%g\n", nu, h * 1e3, d.model().numTets(), u * 1e3, an * 1e3, 100.0 * (u - an) / an);
        }
}

// Test 3. The proxy spreads its force F over the surface nodes within rho of the contact point with weights
// p(r) ~ (1 - r/rho)^2 (mass-proportional, i.e. area-proportional on a regular surface). For a pressure p(r) on the
// surface of an elastic half-space the displacement under the centre is (Boussinesq/Love)
//        w0 = (1 - nu^2)/(pi E) * Int p(r)/r dA = 2 (1 - nu^2)/E * Int_0^rho p(r) dr.
// With p(r) = p0 (1 - r/rho)^2:  F = pi p0 rho^2 / 6  and  Int p dr = p0 rho / 3, hence
//        w0 = 4 F (1 - nu^2) / (pi E rho).
static void test4()
{
    printf("\nTEST 4  proxy contact on the same block: force-indentation consistency (E = 15 kPa)\n");
    printf("  (the contact force is the hand-proxy spring, so this tests the whole contact loop, not only the solid)\n");
    printf("  nu    cell h [mm]  tets    rho [mm]  hand depth   F [N]   FEM dent [mm]  analytic [mm]  err [%%]  tissue max disp [mm]  err [%%]  state\n");
    const double E = 15e3;
    for (double nu : {0.45})
    for (double h : {0.008, 0.005})
    {
        SimConfig c;
        c.nx = (int)std::lround(0.040 / h); c.ny = c.nz = (int)std::lround(0.080 / h); c.h = h;
        c.youngModulus = E; c.poisson = nu; c.density = 1060.0; c.massScale = 4.0; c.dampingPerMass = 30.0;
        c.proxyContact = true; c.couplingStiffness = 800.0; c.contactMaxForce = 5.0; c.proxyMaxHandDepth = 0.0;
        c.probeRadius = 0.00125; c.probeMaxSpeed = 0.6; c.forceLogPath = ""; c.numThreads = 0;
        SurgicalSim sim; sim.init(c);
        const double top = sim.boundsMax().x, frame = 1.0 / 60.0;
        double hx = top + 0.004;
        // approach until contact; the zero of the dent is the rest height of the proxy on the intact top face
        double x0 = 0.0; bool touched = false;
        for (int i = 0; i < 600 && !touched; ++i) { hx -= 0.00005; sim.update(frame, Vec3(hx, 0, 0), false); if (sim.inContact()) touched = true; }
        x0 = top + c.probeRadius;     // the undeformed proxy rests one radius above the flat top face
        if (!touched) { printf("  h=%.1f mm: no contact found\n", h * 1e3); continue; }
        const double rho = sim.loadRadius();
        const double x0hand = x0;   // hand depth is measured below the rest position of the proxy
        for (double depthTarget : {0.001, 0.002})
        {
            const double hxEnd = x0hand - depthTarget;
            while (hx > hxEnd) { hx = (std::max)(hxEnd, hx - 0.0001); sim.update(frame, Vec3(hx, 0, 0), false); }
            double fs = 0.0, xs = 0.0, ds = 0.0, xmin = 1e9, xmax = -1e9; int n = 0;
            for (int i = 0; i < 240; ++i) { Vec3 f = sim.update(frame, Vec3(hx, 0, 0), false); if (i >= 180) { fs += f.length(); xs += sim.probe().x; ds += sim.maxDisplacement(); xmin = (std::min)(xmin, sim.probe().x); xmax = (std::max)(xmax, sim.probe().x); ++n; } }
            const double F = fs / n, dent = x0 - xs / n, pp = xmax - xmin, tissue = ds / n;
            const bool steady = pp < 0.05 * dent + 5e-5;
            const double an = 4.0 * F * (1.0 - nu * nu) / (PI * E * rho);
            printf("  %.2f  %7.1f  %8zu  %7.2f  %8.1f mm  %6.3f  %12.3f  %13.3f  %9.1f  %12.3f  %9.1f  %s\n", nu, h * 1e3, sim.activeTets(), rho * 1e3, depthTarget * 1e3, F, dent * 1e3, an * 1e3, 100.0 * (dent - an) / an, tissue * 1e3, 100.0 * (tissue - an) / an, steady ? "static" : "UNSTEADY (not a valid comparison)");
            if (csv) fprintf(csv, "proxy_indentation,nu=%.2f;h=%g,%zu,rho=%g;hand=%g;F=%g;tissue_max_disp=%g,%g,%g,%g\n", nu, h * 1e3, sim.activeTets(), rho * 1e3, depthTarget * 1e3, F, tissue * 1e3, dent * 1e3, an * 1e3, 100.0 * (dent - an) / an);
        }
        fflush(stdout);
    }
}


// Test 3. The load patch used by the proxy, applied directly as nodal forces: the tissue response without the contact loop.
static void test3()
{
    printf("\nTEST 3  half-space indented by the proxy's load patch applied directly (F = 0.3 N, patch radius rho = 10 mm)\n");
    printf("  analytic: w0 = 4 F (1 - nu^2) / (pi E rho)   [pressure p ~ (1 - r/rho)^2 on an elastic half-space, centre displacement]\n");
    printf("  nu    cell h [mm]   tets     FEM [mm]   analytic [mm]   error [%%]\n");
    const double E = 15e3, F = 0.3, rho = 0.010;
    for (double nu : {0.30, 0.45})
        for (double h : {0.008, 0.006, 0.005, 0.004})
        {
            const int nx = (int)std::lround(0.040 / h), ny = 2 * (int)std::lround(0.040 / h);
            std::vector<Vec3> V; std::vector<TensorMassModel::Tet> T; makeBlockMesh(nx, ny, ny, h, V, T);
            for (Vec3& v : V) v += Vec3(-0.5 * nx * h, -0.5 * ny * h, -0.5 * ny * h);
            TensorMassModel m; m.init(V, T, lameLambda(E, nu), lameMu(E, nu));
            TensorMassDynamics d; d.init(m, 1060.0, 30.0, 1e-4);
            for (int k = 0; k <= ny; ++k) for (int j = 0; j <= ny; ++j) d.fixVertex((k * (ny + 1) + j) * (nx + 1) + 0);
            d.setTimeStep(0.4 * d.estimateCriticalTimeStep());
            std::vector<std::pair<int, double> > w; double sw = 0.0; int centre = -1; double best = 1.0;
            for (int k = 0; k <= ny; ++k) for (int j = 0; j <= ny; ++j)
            {
                const int id = (k * (ny + 1) + j) * (nx + 1) + nx; const Vec3 p = m.restPositions()[(size_t)id];
                const double r = std::sqrt(p.y * p.y + p.z * p.z);
                if (r < best) { best = r; centre = id; }
                if (r < rho) { const double f = 1.0 - r / rho; w.push_back(std::make_pair(id, f * f)); sw += f * f; }
            }
            for (size_t i = 0; i < w.size(); ++i) d.setExternalForce(w[i].first, Vec3(-F * w[i].second / sw, 0, 0));
            std::vector<int> cn(1, centre); relax(d, cn, 0);
            const double u = -(d.positions()[(size_t)centre].x - d.model().restPositions()[(size_t)centre].x);
            const double an = 4.0 * F * (1.0 - nu * nu) / (PI * E * rho);
            printf("  %.2f  %8.1f  %8zu  %9.4f  %12.4f  %10.2f\n", nu, h * 1e3, d.model().numTets(), u * 1e3, an * 1e3, 100.0 * (u - an) / an);
            if (csv) fprintf(csv, "patch_indentation,nu=%.2f,%g,%zu,%g,%g,%g\n", nu, h * 1e3, d.model().numTets(), u * 1e3, an * 1e3, 100.0 * (u - an) / an);
        }
}

int main(int argc, char** argv)
{
    csv = fopen("validation_results.csv", "w");
    if (csv) fprintf(csv, "test,case,h_mm_or_info,tets,fem_or_info,analytic_or_info,error_percent\n");
    printf("Validation of the tensor-mass solver and the proxy contact against analytical solutions\n");
    const std::string only = (argc > 1) ? argv[1] : "";
    if (only.empty() || only == "1") test1();
    if (only.empty() || only == "2") test2();
    if (only.empty() || only == "3") test3();
    if (only.empty() || only == "4") test4();
    printf("\nHOW TO READ THE RESULTS\n  Test 1: solver reproduces the modulus E (error is the end-clamp effect and, for nu = 0.45, volumetric locking of linear tetrahedra).\n  Test 2: constant-strain tetrahedra are too stiff in bending on coarse meshes (shear locking); the error falls as the mesh is refined.\n  Test 3: tissue response to the proxy load patch approaches the half-space solution as the mesh is refined (a finite 40 mm block adds ~10 %% stiffness).\n  Test 4: the contact loop is not static at larger depths: the proxy tip does not follow the deepest point of the dimple (contact on/off cycling).\n          Use the tissue max displacement column for the tissue response and report the contact loop as a qualitative interaction model.\n\nResults also written to validation_results.csv\n");
    if (csv) fclose(csv);
    return 0;
}
