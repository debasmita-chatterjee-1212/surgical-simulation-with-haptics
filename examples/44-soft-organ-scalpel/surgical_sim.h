//==============================================================================
// surgical_sim.h
//
// Step 3 of the surgical simulator: everything the CHAI3D window needs, with no
// CHAI3D dependency (so it is tested on its own):
//   - a soft-tissue block (tensor-mass model, one face clamped)
//   - a spherical probe: palpation by penalty contact, or drilling (cutting)
//   - boundary-surface extraction for rendering (outward-oriented triangles)
//   - removal of fragments that a cut has separated from the clamped base
//   - real-time stepping with a fixed, provably stable time step
//
// The probe position comes from outside (mouse now, Phantom Omni later) and the
// force on the probe is returned, which is exactly what a haptic device needs.
//==============================================================================
#pragma once

#include "tensor_mass_dynamics.h"
#include <algorithm>
#include <array>
#include <map>
#include <queue>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <thread>
#include <fstream>
#include <string>

namespace tmass {

struct SurfaceTri { int a, b, c; };

struct TriangleClosestPoint
{
    Vec3 point;
    double w0, w1, w2;
};

inline TriangleClosestPoint closestPointOnTriangle(const Vec3& p,
                                                    const Vec3& a,
                                                    const Vec3& b,
                                                    const Vec3& c)
{
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return TriangleClosestPoint{a, 1.0, 0.0, 0.0};

    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return TriangleClosestPoint{b, 0.0, 1.0, 0.0};

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        return TriangleClosestPoint{a + ab * v, 1.0 - v, v, 0.0};
    }

    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return TriangleClosestPoint{c, 0.0, 0.0, 1.0};

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        return TriangleClosestPoint{a + ac * w, 1.0 - w, 0.0, w};
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return TriangleClosestPoint{b + (c - b) * w, 0.0, 1.0 - w, w};
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    return TriangleClosestPoint{a + ab * v + ac * w, 1.0 - v - w, v, w};
}


// Boundary of the active tetrahedra as outward-pointing triangles (the outward
// side of a cavity left by a cut points into the cavity).
inline std::vector<SurfaceTri> extractBoundaryTriangles(const TensorMassModel& m)
{
    struct Entry { int count; SurfaceTri tri; };
    std::map<std::array<int, 3>, Entry> faces;
    const std::vector<Vec3>& rest = m.restPositions();
    for (size_t t = 0; t < m.numTets(); t++)
    {
        if (!m.isActive((int)t)) continue;
        const TensorMassModel::Tet& T = m.tets()[t];
        for (int o = 0; o < 4; o++)
        {
            int f[3], n = 0;
            for (int k = 0; k < 4; k++) if (k != o) f[n++] = T[k];
            const Vec3 nrm = cross(rest[f[1]] - rest[f[0]], rest[f[2]] - rest[f[0]]);
            if (dot(nrm, rest[T[o]] - rest[f[0]]) > 0.0) std::swap(f[1], f[2]);   // face the normal away from the tet
            std::array<int, 3> key = {{f[0], f[1], f[2]}};
            std::sort(key.begin(), key.end());
            auto it = faces.find(key);
            if (it == faces.end()) { Entry e; e.count = 1; e.tri = SurfaceTri{f[0], f[1], f[2]}; faces[key] = e; }
            else it->second.count++;
        }
    }
    std::vector<SurfaceTri> out;
    for (std::map<std::array<int, 3>, Entry>::const_iterator it = faces.begin(); it != faces.end(); ++it)
        if (it->second.count == 1) out.push_back(it->second.tri);
    return out;
}

struct SimConfig
{
    int nx = 8, ny = 8, nz = 6;          // cells; block is centred on the origin
    double h = 0.01;                     // cell size [m]  -> 8 x 8 x 6 cm
    double youngModulus = 7.0e3;         // [Pa] soft tissue / gel (lower = softer, deforms more)
    double poisson = 0.35;
    double density = 1000.0;             // [kg/m^3]
    double dampingPerMass = 25.0;        // [1/s]
    double probeRadius = 0.00125;        // [m] contact-proxy radius (2.5 mm tip diameter)
    double drillRadius = 0.010;          // [m] tetrahedra whose centroid is inside are cut
    double contactStiffness = 400.0;     // [N/m] penalty stiffness of the probe
    double contactMaxForce = 3.0;        // [N] per vertex
    // > 0: area-weighted contact on the surface triangles (3 sample points per triangle, force
    // spread to the triangle's vertices). The penalty stiffness is then per unit area [N/m^3] and
    // the result does not depend on how coarse the surface mesh is. 0 = per-vertex contact.
    double contactStiffnessPerArea = 8.0e6; // calibrated so coarse block triangles still produce measurable contact
    // Progressive cutting: energy density [J/m^3] and available cutting power [W].
    // These are prototype calibration parameters, not validated tissue constants.
    double cuttingEnergyPerVolume = 2.0e3;
    double cuttingPower = 0.05;
    double probeMaxSpeed = 0.6;          // [m/s] the probe is rate-limited (no teleporting)
    double safety = 0.4;                 // dt = safety * critical time step
    double maxDt = 2.0e-3;               // [s]
    int maxStepsPerUpdate = 400;
    int numThreads = 0;                  // physics threads; 0 = automatic (cores - 1, at most 6), 1 = serial

    // Optional: an arbitrary tetrahedral mesh instead of the block (see mesh_*.h, make_tet_mesh.py).
    // verts: 3 doubles per vertex [m]; tets: 4 ints per tetrahedron; fixed: clamped vertex indices.
    const double* meshVerts = nullptr;
    const int* meshTets = nullptr;
    const int* meshFixed = nullptr;
    int meshNumVerts = 0, meshNumTets = 0, meshNumFixed = 0;

    // Mass scaling for stiff tissue: density is multiplied by this factor, which raises the stable
    // time step by sqrt(massScale). Quasi-static behaviour (palpation force/displacement) is
    // unchanged; only the fast dynamics are slowed down. 1 = physical density.
    double massScale = 1.0;

    // Rigid mode (thesis Ch. 6 hybrid, uncut stiff bone): the bone's elastic response to haptic-level
    // forces is far below a micron (see test_organs: static solve), so the uncut part is not
    // integrated at all. The probe still feels a force-limited penalty contact (contactStiffnessPerArea,
    // total clipped to contactMaxForce), and cutting still removes tetrahedra. Real-time at any stiffness.
    bool rigid = false;

    // Capsule contact: the probe is a sphere (probeRadius) at the tip plus a cylindrical shaft of the
    // same radius that extends probeShaftLength along +x (away from the tissue). A shaft that touches
    // the tissue surface keeps it from slipping round a small tip when the probe is pressed deep.
    double probeShaftLength = 0.012;     // [m]
    // How far BELOW the intact outer surface the probe axis may be pushed when palpating [m].
    // 0 = old behaviour (the tip centre cannot go below the surface).
    double maxIndentation = 0.0;
    // Rigid (bone) drilling: the burr cannot enter uncut bone. Its reaction force is the contact
    // force against the bone, limited to contactMaxForce, shown and logged live. PLACEHOLDER until a
    // measured drilling-force model (cutting force vs. speed, depth, burr size) is available.
    // Contact radius of the burr as a fraction of the cutting radius: the cut cavity is slightly wider than
    // the contact sphere, so the burr does not feel the walls of its own hole as resistance.
    double burrContactFraction = 0.8;
    // Proxy (god-object) contact for soft organs. The instrument is a small sphere that can never enter the
    // tissue: it is a proxy that stays on the surface while the hand (mouse / device) may go deeper. The
    // hand is tied to the proxy by a virtual spring (couplingStiffness); the spring force, along the
    // surface normal and limited to contactMaxForce, is both the force the user feels and the load that
    // is applied to the tissue where the sphere touches it. There is no penalty stiffness and no
    // tunnelling, so the force does not drop to zero when the tip is pressed deep, and the contact no longer
    // limits the time step. See: Zilles & Salisbury 1995 (god-object); Ruspini et al. 1997 (proxy).
    bool   proxyContact = false;
    double couplingStiffness = 800.0;    // [N/m] hand-proxy spring
    // The proxy load is spread over the surface nodes within max(loadRadiusTips * radius, loadRadiusEdges * mean
    // surface edge) of the contact point (weights (1 - d/rho)^2, mass-proportional). A point load on one cell makes
    // an unphysical needle-like pit in a mesh whose cells are larger than the tip; the real contact patch of a
    // small tip on soft tissue (Hertz) is several mm wide.
    double loadRadiusTips = 3.0, loadRadiusEdges = 1.2;
    // Largest depth [m] the hand may go below the intact surface (0 = no limit). A linear-elastic model cannot
    // represent the dimple of a small tip beyond ~50 % strain, so the hand stops where the model stays valid.
    double proxyMaxHandDepth = 0.0;
    int    smoothingPasses = 8;          // Taubin passes (lambda|mu) on the exposed surface; 0 = raw staircase
    double contactTolerance = 0.6;      // contact counts while the surface is within (1 + tol) * radius of the proxy
    // Force log (CSV, appended): empty string = off.
    std::string forceLogPath = "force_log.csv";
};

class SurgicalSim
{
public:
    void init(const SimConfig& cfg)
    {
        if (!(cfg.youngModulus > 0.0) || !std::isfinite(cfg.youngModulus) ||
            !(cfg.poisson > -1.0 && cfg.poisson < 0.5) || !std::isfinite(cfg.poisson) ||
            !(cfg.density > 0.0) || !std::isfinite(cfg.density) ||
            !(cfg.massScale > 0.0) || !std::isfinite(cfg.massScale) ||
            cfg.dampingPerMass < 0.0 || !std::isfinite(cfg.dampingPerMass) ||
            !(cfg.safety > 0.0 && cfg.safety <= 1.0) ||
            !(cfg.maxDt > 0.0) || !std::isfinite(cfg.maxDt) || cfg.maxStepsPerUpdate <= 0 ||
            !(cfg.probeRadius > 0.0) || !std::isfinite(cfg.probeRadius) ||
            !(cfg.drillRadius > 0.0) || !std::isfinite(cfg.drillRadius) ||
            cfg.contactStiffness < 0.0 || !std::isfinite(cfg.contactStiffness) ||
            cfg.contactStiffnessPerArea < 0.0 || !std::isfinite(cfg.contactStiffnessPerArea) ||
            cfg.contactMaxForce < 0.0 || !std::isfinite(cfg.contactMaxForce) ||
            !(cfg.cuttingEnergyPerVolume > 0.0) || !std::isfinite(cfg.cuttingEnergyPerVolume) ||
            !(cfg.cuttingPower >= 0.0) || !std::isfinite(cfg.cuttingPower) ||
            cfg.probeMaxSpeed < 0.0 || !std::isfinite(cfg.probeMaxSpeed))
            throw std::invalid_argument("SurgicalSim: invalid material, mass, contact, or timestep configuration");
        m_cfg = cfg;
        std::vector<Vec3> V;
        std::vector<TensorMassModel::Tet> T;
        const bool customMesh = (cfg.meshVerts != nullptr && cfg.meshTets != nullptr);
        if (customMesh)
        {
            if (cfg.meshNumVerts <= 0 || cfg.meshNumTets <= 0 || cfg.meshNumFixed < 0 ||
                (cfg.meshNumFixed > 0 && cfg.meshFixed == nullptr))
                throw std::invalid_argument("SurgicalSim: custom mesh counts/pointers are inconsistent");
            for (int i = 0; i < cfg.meshNumVerts; i++)
                V.push_back(Vec3(cfg.meshVerts[3 * i], cfg.meshVerts[3 * i + 1], cfg.meshVerts[3 * i + 2]));
            for (int t = 0; t < cfg.meshNumTets; t++)
            {
                TensorMassModel::Tet tt = {{cfg.meshTets[4 * t], cfg.meshTets[4 * t + 1],
                                            cfg.meshTets[4 * t + 2], cfg.meshTets[4 * t + 3]}};
                T.push_back(tt);
            }
        }
        else
        {
            makeBlockMesh(cfg.nx, cfg.ny, cfg.nz, cfg.h, V, T);
            const Vec3 shift(-0.5 * cfg.nx * cfg.h, -0.5 * cfg.ny * cfg.h, -0.5 * cfg.nz * cfg.h);
            for (size_t i = 0; i < V.size(); i++) V[i] += shift;
        }

        m_boundsMin = Vec3(1e300, 1e300, 1e300);
        m_boundsMax = Vec3(-1e300, -1e300, -1e300);
        for (size_t i = 0; i < V.size(); i++)
        {
            m_boundsMin = Vec3((std::min)(m_boundsMin.x, V[i].x), (std::min)(m_boundsMin.y, V[i].y), (std::min)(m_boundsMin.z, V[i].z));
            m_boundsMax = Vec3((std::max)(m_boundsMax.x, V[i].x), (std::max)(m_boundsMax.y, V[i].y), (std::max)(m_boundsMax.z, V[i].z));
        }

        TensorMassModel model;
        model.init(V, T, lameLambda(cfg.youngModulus, cfg.poisson), lameMu(cfg.youngModulus, cfg.poisson));
        m_dyn.init(model, cfg.density * cfg.massScale, cfg.dampingPerMass, 1e-4);
        {
            int th = cfg.numThreads;
            if (th <= 0) th = (std::max)(1, (std::min)(6, (int)std::thread::hardware_concurrency() - 1));
            m_dyn.setThreads(th);
        }
        m_cutsSinceFullEstimate = 0;

        if (customMesh)
        {
            for (int i = 0; i < cfg.meshNumFixed; i++) m_dyn.fixVertex(cfg.meshFixed[i]);
        }
        else
        {
            // clamp the back face (x = minimum)
            for (int k = 0; k <= cfg.nz; k++)
                for (int j = 0; j <= cfg.ny; j++) m_dyn.fixVertex((k * (cfg.ny + 1) + j) * (cfg.nx + 1) + 0);
        }

        rebuildSurface();
        // the intact outer surface (rest positions): used to limit how deep palpation may push
        m_intactTris.clear();
        for (size_t i = 0; i < m_surface.size(); i++)
        {
            m_intactTris.push_back(m_dyn.model().restPositions()[m_surface[i].a]);
            m_intactTris.push_back(m_dyn.model().restPositions()[m_surface[i].b]);
            m_intactTris.push_back(m_dyn.model().restPositions()[m_surface[i].c]);
        }

        recomputeTimeStep();

        m_deepestDrillX = 1e9;
        m_deepestClearedX = 1e9;
        m_simTime = 0.0;
        m_logRows = 0;
        openForceLog();
        m_acc = 0.0;
        m_probeInit = false;
        m_surfaceChanged = true;
        m_lastReaction = Vec3();
        m_numContacts = 0;
        m_numTets0 = m_dyn.model().numActiveTets();
        m_cutWork.assign(m_dyn.model().numTets(), 0.0);
        m_lastMaxDisplacement = 0.0;
    }

    void reset() { init(m_cfg); }

    // Advance the simulation by wallDt seconds of real time.
    //   target   : desired probe centre [m]
    //   drilling : true -> cut tissue inside the drill sphere; false -> palpate (contact force)
    // Returns the force the tissue exerts on the probe [N] (send this to the haptic device).
    Vec3 update(double wallDt, const Vec3& target, bool drilling)
    {
        if (!(wallDt > 0.0) || !std::isfinite(wallDt)) wallDt = 0.0;
        wallDt = (std::min)(wallDt, 0.05);

        // Palpation (no cutting) cannot push the probe centre deeper than the intact
        // outer surface beneath it: beyond that the penalty contact would let the tissue
        // slip around the probe (it "pops out"). Where the user has drilled, the probe may
        // follow into the pocket.
        const bool proxy = m_cfg.proxyContact && !m_cfg.rigid;
        Vec3 goal = target;
        if (!drilling && proxy && m_cfg.proxyMaxHandDepth > 0.0)
        {
            const double xl = intactSurfaceX(goal) - m_cfg.proxyMaxHandDepth;
            if (goal.x < xl) goal.x = xl;
        }
        if (!drilling && !proxy)
        {
            const double xLimit = (std::min)(intactSurfaceX(goal) - m_cfg.maxIndentation, m_deepestDrillX);
            if (goal.x < xLimit) goal.x = xLimit;
        }
        else if (m_cfg.rigid)
        {
            // rigid bone: the burr can only advance into material it has already cleared
            const double xLimit = (std::min)(intactSurfaceX(goal), m_deepestClearedX);
            if (goal.x < xLimit) goal.x = xLimit;
        }

        m_hand = goal;
        if (!m_probeInit) { m_probe = goal; m_probeInit = true; }
        else if (!(proxy && !drilling))        // the proxy moves itself, step by step, in proxyStep()
        {
            Vec3 d = goal - m_probe;
            const double len = d.length();
            const double maxMove = m_cfg.probeMaxSpeed * wallDt;
            if (len > maxMove && len > 0.0) d = d * (maxMove / len);
            m_probe += d;
        }

        if (drilling) { m_deepestDrillX = (std::min)(m_deepestDrillX, m_probe.x); drill(wallDt); }

        if (m_cfg.rigid)
        {
            m_dyn.clearExternalForces();
            Vec3 f = drilling ? applyContact(m_cfg.burrContactFraction * m_cfg.drillRadius, m_cfg.probeShaftLength) : applyContact(m_cfg.probeRadius, m_cfg.probeShaftLength);
            const double len = f.length();
            if (len > m_cfg.contactMaxForce) f = f * (m_cfg.contactMaxForce / len);   // device-style force limit
            m_lastReaction = f;
            m_simTime += wallDt;
            logForce(drilling);
            return m_lastReaction;
        }

        m_acc += wallDt;
        Vec3 sum;
        int steps = 0;
        if (proxy && !drilling) refreshNear(true);
        while (m_acc >= m_dt && steps < m_cfg.maxStepsPerUpdate)
        {
            m_dyn.clearExternalForces();
            if (!drilling) sum += proxy ? proxyStep(m_dt) : applyContact(m_cfg.probeRadius, m_cfg.probeShaftLength);
            m_dyn.step();
            m_acc -= m_dt;
            steps++;
        }
        if (proxy && !drilling) proxyRecover();
        if (steps >= m_cfg.maxStepsPerUpdate) m_acc = 0.0;       // fall behind instead of spiralling
        if (steps > 0) m_lastReaction = sum * (1.0 / steps);
        else if (!drilling) m_lastReaction = Vec3();
        if (drilling) { m_lastReaction = Vec3(); m_numContacts = 0; }

        m_lastMaxDisplacement = m_dyn.maxDisplacement();
        m_simTime += wallDt;

        // safety net: if anything ever diverges, start over instead of showing garbage
        if (!std::isfinite(m_dyn.maxDisplacement()))
        {
            const bool wasInit = m_probeInit;
            const Vec3 p = m_probe;
            init(m_cfg);
            m_probeInit = wasInit;
            m_probe = p;
            m_resetCount++;
        }
        logForce(drilling);
        return m_lastReaction;
    }

    // true once after the visible surface changed (a cut happened or the sim was reset)
    bool consumeSurfaceChanged() { const bool c = m_surfaceChanged; m_surfaceChanged = false; return c; }

    const std::vector<Vec3>& positions() const { return m_dyn.positions(); }
    // positions for drawing: the simulated positions plus the surface-smoothing offset (identical when smoothing is off)
    const std::vector<Vec3>& displayPositions()
    {
        const std::vector<Vec3>& P = m_dyn.positions();
        m_displayPos = P;
        if (m_smoothOffset.size() == P.size())
            for (size_t q = 0; q < m_exposed.size(); ++q) { const size_t v = (size_t)m_exposed[q]; m_displayPos[v] += m_smoothOffset[v]; }
        return m_displayPos;
    }
    const std::vector<SurfaceTri>& surface() const { return m_surface; }
    const Vec3& probe() const { return m_probe; }
    double timeStep() const { return m_dt; }
    double meshCriticalDt() const { return m_dtMesh; }
    double contactCriticalDt() const { return m_dtContact; }
    size_t activeTets() const { return m_dyn.model().numActiveTets(); }
    size_t initialTets() const { return m_numTets0; }
    int numContacts() const { return m_numContacts; }
    int proxyRecoveries() const { return m_recoveries; }
    bool proxyIsInside() const { return pointInside(m_probe); }
    bool inContact() const { return m_inContact; }
    // strongest contact of the last proxy step: surface triangle index, barycentric weights of its a,b,c vertices, outward normal
    bool contactInfo(int& tri, double bary[3], Vec3& normal) const
    {
        if (!m_inContact || m_lastTri < 0) return false;
        tri = m_lastTri; bary[0] = m_lastBary[0]; bary[1] = m_lastBary[1]; bary[2] = m_lastBary[2]; normal = m_lastN;
        return true;
    }
    // upper limit of the proxy load / reaction force [N] (e.g. while cutting); negative = the organ's contactMaxForce
    void setLoadLimit(double limit) { m_loadLimit = limit; }
    double rawSpringForce() const { return m_rawSpringForce; }
    const Vec3& hand() const { return m_hand; }
    double maxDisplacement() const { return m_dyn.maxDisplacement(); }
    double lastMaxDisplacement() const { return m_lastMaxDisplacement; }
    double cuttingProgress(int tet) const { return (tet >= 0 && (size_t)tet < m_cutWork.size()) ? m_cutWork[(size_t)tet] : 0.0; }
    int resetCount() const { return m_resetCount; }
    const TensorMassDynamics& dynamics() const { return m_dyn; }
    const SimConfig& config() const { return m_cfg; }
    const Vec3& boundsMin() const { return m_boundsMin; }      // bounding box of the intact organ
    const Vec3& boundsMax() const { return m_boundsMax; }

private:
    // Recompute both mesh and contact stability limits. This is also called after
    // cutting because removing tetrahedra changes the lumped masses and stiffness.
    //
    // Removing a tetrahedron removes stiffness AND lumped mass, so the stable time step is
    // not guaranteed to stay the same. After a cut the mesh limit is therefore re-estimated,
    // but cheaply: a warm-started power iteration (100 steps instead of 400). Every
    // 20 cuts a full cold estimate replaces it. The contact limit (cheap) is always redone.
    void recomputeTimeStep(bool afterCut = false)
    {
        const double inf = std::numeric_limits<double>::infinity();
        if (m_cfg.rigid) m_dtMesh = inf;
        else if (afterCut && ++m_cutsSinceFullEstimate < 20) m_dtMesh = m_dyn.estimateCriticalTimeStep(100, true);
        else { m_cutsSinceFullEstimate = 0; m_dtMesh = m_dyn.estimateCriticalTimeStep(); }
        double dtContact = inf;
        if (m_cfg.proxyContact && !m_cfg.rigid)
        {
            // The virtual spring loads a patch of surface nodes (weights (1 - d/rho)^2, mass-proportional), so the
            // relevant mass is that of a patch, not of the lightest corner node. A patch holds at least ~3 nodes of
            // typical mass: use 3 x the 25th-percentile mass of the exposed nodes.
            std::vector<double> ms;
            for (size_t q = 0; q < m_exposed.size(); ++q)
            {
                const int v = m_exposed[q];
                const double mass = m_dyn.masses()[(size_t)v];
                if (!m_dyn.isFixed(v) && mass > 0.0) ms.push_back(mass);
            }
            double minMass = inf;
            if (!ms.empty()) { std::sort(ms.begin(), ms.end()); minMass = 3.0 * ms[ms.size() / 4]; }
            if (std::isfinite(minMass)) dtContact = 2.0 * std::sqrt(minMass / m_cfg.couplingStiffness);
        }
        else if (m_cfg.contactStiffnessPerArea > 0.0)
        {
            std::vector<double> tributary(m_dyn.model().numVertices(), 0.0);
            const std::vector<Vec3>& R = m_dyn.model().restPositions();
            for (size_t i = 0; i < m_surface.size(); ++i)
            {
                const SurfaceTri& tri = m_surface[i];
                const double area = 0.5 * cross(R[tri.b] - R[tri.a], R[tri.c] - R[tri.a]).length();
                tributary[(size_t)tri.a] += area / 3.0;
                tributary[(size_t)tri.b] += area / 3.0;
                tributary[(size_t)tri.c] += area / 3.0;
            }
            double worst = 0.0;
            for (size_t q = 0; q < m_exposed.size(); ++q)
            {
                const int v = m_exposed[q];
                const double mass = m_dyn.masses()[(size_t)v];
                if (m_dyn.isFixed(v) || !(mass > 0.0)) continue;
                worst = (std::max)(worst, m_cfg.contactStiffnessPerArea * tributary[(size_t)v] / mass);
            }
            if (worst > 0.0 && std::isfinite(worst)) dtContact = 2.0 / std::sqrt(worst);
        }
        else if (m_cfg.contactStiffness > 0.0)
        {
            double minMass = inf;
            for (size_t q = 0; q < m_exposed.size(); ++q)
            {
                const int v = m_exposed[q];
                const double mass = m_dyn.masses()[(size_t)v];
                if (!m_dyn.isFixed(v) && mass > 0.0) minMass = (std::min)(minMass, mass);
            }
            if (std::isfinite(minMass)) dtContact = 2.0 * std::sqrt(minMass / m_cfg.contactStiffness);
        }
        m_dtContact = dtContact;
        if (m_cfg.rigid) m_dt = 1.0e-3;
        else
        {
            const double stable = (std::min)(m_dtMesh, m_dtContact);
            m_dt = std::isfinite(stable) ? (std::min)(m_cfg.safety * stable, m_cfg.maxDt) : m_cfg.maxDt;
            if (!(m_dt > 0.0) || !std::isfinite(m_dt))
                throw std::runtime_error("SurgicalSim: could not determine a valid integration timestep");
        }
        m_dyn.setTimeStep(m_dt);
        m_acc = 0.0;
    }

    // x of the intact outer surface directly beneath the probe axis (ray cast along -x on the
    // initial surface triangles); -infinity when the probe is beside the organ (nothing to limit).
    double intactSurfaceX(const Vec3& p) const
    {
        double best = -1e300;
        for (size_t i = 0; i < m_intactTris.size(); i += 3)
        {
            const Vec3 &A = m_intactTris[i], &B = m_intactTris[i + 1], &C = m_intactTris[i + 2];
            const double d = (B.y - A.y) * (C.z - A.z) - (C.y - A.y) * (B.z - A.z);
            if (std::fabs(d) < 1e-18) continue;
            const double py = p.y - A.y, pz = p.z - A.z;
            const double u = (py * (C.z - A.z) - (C.y - A.y) * pz) / d;
            const double v = ((B.y - A.y) * pz - py * (B.z - A.z)) / d;
            if (u < 0.0 || v < 0.0 || u + v > 1.0) continue;
            best = (std::max)(best, A.x + u * (B.x - A.x) + v * (C.x - A.x));
        }
        return best;
    }

    // ---- proxy contact -------------------------------------------------------------------------
    // surface triangles close to the proxy (bounding boxes on the current, deformed positions)
    void refreshNear(bool force)
    {
        const double R = m_cfg.probeRadius;
        const double margin = 0.006 + 3.0 * R;
        if (!force && m_nearValid && (m_probe - m_nearCenter).length() < 0.25 * margin) return;
        m_near.clear();
        const std::vector<Vec3>& P = m_dyn.positions();
        const double lo = margin;
        for (size_t i = 0; i < m_surface.size(); ++i)
        {
            const SurfaceTri& t = m_surface[i];
            const Vec3 &a = P[(size_t)t.a], &b = P[(size_t)t.b], &c = P[(size_t)t.c];
            if ((std::min)({a.x, b.x, c.x}) > m_probe.x + lo || (std::max)({a.x, b.x, c.x}) < m_probe.x - lo ||
                (std::min)({a.y, b.y, c.y}) > m_probe.y + lo || (std::max)({a.y, b.y, c.y}) < m_probe.y - lo ||
                (std::min)({a.z, b.z, c.z}) > m_probe.z + lo || (std::max)({a.z, b.z, c.z}) < m_probe.z - lo) continue;
            m_near.push_back((int)i);
        }
        m_nearVerts.clear();
        for (size_t q = 0; q < m_near.size(); ++q)
        {
            const SurfaceTri& t = m_surface[(size_t)m_near[q]];
            m_nearVerts.push_back(t.a); m_nearVerts.push_back(t.b); m_nearVerts.push_back(t.c);
        }
        std::sort(m_nearVerts.begin(), m_nearVerts.end());
        m_nearVerts.erase(std::unique(m_nearVerts.begin(), m_nearVerts.end()), m_nearVerts.end());
        m_nearCenter = m_probe;
        m_nearValid = true;
    }

    // One physics step of the proxy: move towards the hand (speed-limited), push out of the surface,
    // load the tissue with the spring force along the contact normal. Returns the force on the hand.
    Vec3 proxyStep(double dt)
    {
        const double R = m_cfg.probeRadius;
        // 1. move towards the hand
        Vec3 d = m_hand - m_probe;
        const double len = d.length();
        const double maxMove = m_cfg.probeMaxSpeed * dt;
        if (len > maxMove && len > 0.0) d = d * (maxMove / len);
        // god-object rule: while touching, the proxy never moves INTO the surface (it may slide or retreat); it
        // follows a receding surface only in the small steps it takes once the contact has been lost.
        for (int pass = 0; pass < 2; ++pass)
            for (size_t k = 0; k < m_prevNormals.size(); ++k)
            {
                const double c = dot(d, m_prevNormals[k]);
                if (c < 0.0) d -= m_prevNormals[k] * c;
            }
        m_probe += d;
        refreshNear(false);
        const std::vector<Vec3>& P = m_dyn.positions();
        const bool useOff = (m_smoothOffset.size() == P.size());
        auto sp = [&](int v) -> Vec3 { return useOff ? P[(size_t)v] + m_smoothOffset[(size_t)v] : P[(size_t)v]; };
        // 2. push the sphere out of the surface (never by more than R, so it cannot jump through a thin wall)
        for (int it = 0; it < 4; ++it)
        {
            bool moved = false;
            for (size_t q = 0; q < m_near.size(); ++q)
            {
                const SurfaceTri& t = m_surface[(size_t)m_near[q]];
                const Vec3 a = sp(t.a), b = sp(t.b), c = sp(t.c);
                const TriangleClosestPoint cp = closestPointOnTriangle(m_probe, a, b, c);
                const Vec3 delta = m_probe - cp.point;
                const double dist = delta.length();
                if (!(dist < R)) continue;
                Vec3 n;
                double push = R - dist;
                const Vec3 nf0 = cross(b - a, c - a); const double nfl = nf0.length();
                if (dist > 1e-10 && nfl > 1e-14 && dot(delta, nf0) < -0.5 * dist * nfl)
                {
                    // the centre is BEHIND this face (inside the organ): the way out is along the outward face normal
                    n = nf0 * (1.0 / nfl);
                    push = R + dist;
                }
                else if (dist > 1e-10) n = delta * (1.0 / dist);
                else
                {
                    const Vec3 nn = cross(b - a, c - a); const double nl = nn.length();
                    if (!(nl > 1e-14)) continue;
                    n = nn * (1.0 / nl);                       // outward
                }
                m_probe += n * (std::min)(push, 2.0 * R);
                moved = true;
            }
            if (!moved) break;
        }
        // 3. contacts: every triangle within (1 + tol) R
        m_contactsTmp.clear();
        double wsum = 0.0; Vec3 nsum;
        const double reach = R * (1.0 + m_cfg.contactTolerance);
        for (size_t q = 0; q < m_near.size(); ++q)
        {
            const SurfaceTri& t = m_surface[(size_t)m_near[q]];
            const Vec3 a = sp(t.a), b = sp(t.b), c = sp(t.c);
            const TriangleClosestPoint cp = closestPointOnTriangle(m_probe, a, b, c);
            const Vec3 delta = m_probe - cp.point;
            const double dist = delta.length();
            if (!(dist < reach) || !(dist > 1e-10)) continue;
            ContactSample s; s.tri = m_near[q]; s.cp = cp; s.w = (reach - dist); s.n = delta * (1.0 / dist);
            m_contactsTmp.push_back(s);
            wsum += s.w; nsum += s.n * s.w;
        }
        m_numContacts = (int)m_contactsTmp.size();
        m_inContact = m_numContacts > 0;
        m_prevNormals.clear();
        for (size_t k = 0; k < m_contactsTmp.size() && k < 8; ++k) m_prevNormals.push_back(m_contactsTmp[k].n);
        if (!m_inContact || !(wsum > 0.0)) return Vec3();
        const double nl = nsum.length();
        if (!(nl > 1e-12)) return Vec3();
        const Vec3 nbar = nsum * (1.0 / nl);                  // outward
        // 4. spring force along the normal, clipped to the device limit
        const Vec3 gap = m_hand - m_probe;
        double pen = -dot(gap, nbar);                         // how far the hand is inside the tissue
        if (pen < 0.0) pen = 0.0;
        double mag = m_cfg.couplingStiffness * pen;
        m_rawSpringForce = mag;
        const double cap = (m_loadLimit >= 0.0) ? m_loadLimit : m_cfg.contactMaxForce;
        if (cap > 0.0 && mag > cap) mag = cap;
        {   // remember the strongest contact (triangle + barycentric position + normal): used by the scalpel
            size_t bi = 0;
            for (size_t k = 1; k < m_contactsTmp.size(); ++k) if (m_contactsTmp[k].w > m_contactsTmp[bi].w) bi = k;
            m_lastTri = m_contactsTmp[bi].tri;
            m_lastBary[0] = m_contactsTmp[bi].cp.w0; m_lastBary[1] = m_contactsTmp[bi].cp.w1; m_lastBary[2] = m_contactsTmp[bi].cp.w2;
            m_lastN = nbar;
        }
        const Vec3 onTissue = nbar * (-mag);
        // contact point: weighted mean of the closest points
        Vec3 q;
        for (size_t k = 0; k < m_contactsTmp.size(); ++k) q += m_contactsTmp[k].cp.point * (m_contactsTmp[k].w / wsum);
        const double rho = (std::max)(m_cfg.loadRadiusTips * R, m_cfg.loadRadiusEdges * m_meanEdge);
        m_loadW.clear(); double usum = 0.0;
        for (size_t k = 0; k < m_nearVerts.size(); ++k)
        {
            const int v = m_nearVerts[k];
            if (m_dyn.isFixed(v)) continue;
            const double d = (sp(v) - q).length();
            if (!(d < rho)) continue;
            const double f = 1.0 - d / rho;
            const double u = f * f * m_dyn.masses()[(size_t)v];
            if (!(u > 0.0)) continue;
            m_loadW.push_back(std::make_pair(v, u)); usum += u;
        }
        if (!(usum > 0.0)) return Vec3();
        for (size_t k = 0; k < m_loadW.size(); ++k)
            m_dyn.addExternalForce(m_loadW[k].first, onTissue * (m_loadW[k].second / usum));
        return nbar * mag;                                    // reaction on the user: pushes the hand out
    }

    // Is point p inside the organ? Ray parity along (almost) +x over the current surface.
    bool pointInside(const Vec3& p) const
    {
        const std::vector<Vec3>& P = m_dyn.positions();
        const bool useOff = (m_smoothOffset.size() == P.size());   // same (smoothed) surface the contact uses
        auto sp = [&](int v) -> Vec3 { return useOff ? P[(size_t)v] + m_smoothOffset[(size_t)v] : P[(size_t)v]; };
        const Vec3 dir(1.0, 1.3e-3, 2.1e-3);
        int cnt = 0;
        for (size_t i = 0; i < m_surface.size(); ++i)
        {
            const SurfaceTri& t = m_surface[i];
            const Vec3 a = sp(t.a), b = sp(t.b), c = sp(t.c);
            if ((std::max)({a.x, b.x, c.x}) < p.x) continue;
            const Vec3 e1 = b - a, e2 = c - a;
            const Vec3 h = cross(dir, e2);
            const double det = dot(e1, h);
            if (std::fabs(det) < 1e-18) continue;
            const double inv = 1.0 / det;
            const Vec3 sv = p - a;
            const double u = dot(sv, h) * inv;
            if (u < 0.0 || u > 1.0) continue;
            const Vec3 qv = cross(sv, e1);
            const double v = dot(dir, qv) * inv;
            if (v < 0.0 || u + v > 1.0) continue;
            if (dot(e2, qv) * inv > 0.0) cnt++;
        }
        return (cnt & 1) != 0;
    }

    // Safety net, once per frame: the proxy is a rigid sphere and must never be inside the organ. If it is
    // (e.g. the tissue moved through it), put it back on the nearest surface point, outside.
    void proxyRecover()
    {
        if (!pointInside(m_probe)) { m_goodProbe = m_probe; m_goodValid = true; return; }
        const std::vector<Vec3>& P = m_dyn.positions();
        const bool useOff = (m_smoothOffset.size() == P.size());
        auto sp = [&](int v) -> Vec3 { return useOff ? P[(size_t)v] + m_smoothOffset[(size_t)v] : P[(size_t)v]; };
        double best = 1e30; Vec3 bp, bn;
        for (size_t i = 0; i < m_surface.size(); ++i)
        {
            const SurfaceTri& t = m_surface[i];
            const Vec3 a = sp(t.a), b = sp(t.b), c = sp(t.c);
            const TriangleClosestPoint cp = closestPointOnTriangle(m_probe, a, b, c);
            const double d = (m_probe - cp.point).length();
            if (d < best)
            {
                const Vec3 nf = cross(b - a, c - a); const double nl = nf.length();
                if (!(nl > 1e-14)) continue;
                best = d; bp = cp.point; bn = nf * (1.0 / nl);
            }
        }
        if (best > 1e29) return;
        m_probe = bp + bn * (1.05 * m_cfg.probeRadius);
        m_prevNormals.clear(); m_prevNormals.push_back(bn);
        m_nearValid = false;
        m_recoveries++;
    }

    // Closest-point triangle contact. Unlike fixed barycentric samples, this detects
    // contact anywhere on a triangle, including between widely spaced mesh vertices.
    Vec3 applyContact(double radius, double shaftLength)
    {
        Vec3 reaction;
        int contacts = 0;
        const std::vector<Vec3>& P = m_dyn.positions();
        for (size_t i = 0; i < m_surface.size(); ++i)
        {
            const SurfaceTri& tri = m_surface[i];
            const int v[3] = {tri.a, tri.b, tri.c};
            const Vec3& a = P[(size_t)v[0]];
            const Vec3& b = P[(size_t)v[1]];
            const Vec3& c = P[(size_t)v[2]];
            const double area = 0.5 * cross(b - a, c - a).length();
            if (!(area > 1e-16)) continue;

            // cheap reject: bounding box of the triangle against the bounding box of the capsule
            if ((std::min)({a.x, b.x, c.x}) > m_probe.x + shaftLength + radius || (std::max)({a.x, b.x, c.x}) < m_probe.x - radius ||
                (std::min)({a.y, b.y, c.y}) > m_probe.y + radius || (std::max)({a.y, b.y, c.y}) < m_probe.y - radius ||
                (std::min)({a.z, b.z, c.z}) > m_probe.z + radius || (std::max)({a.z, b.z, c.z}) < m_probe.z - radius) continue;

            // closest point of the triangle to the capsule axis: sample the tip and the shaft
            TriangleClosestPoint cp = closestPointOnTriangle(m_probe, a, b, c);
            Vec3 axisPoint = m_probe;
            double dist = (cp.point - m_probe).length();
            for (double s = radius; s <= shaftLength + 1e-12 && dist > 1e-10; s += radius)
            {
                const Vec3 q(m_probe.x + s, m_probe.y, m_probe.z);
                const TriangleClosestPoint cq = closestPointOnTriangle(q, a, b, c);
                const double dq = (cq.point - q).length();
                if (dq < dist) { dist = dq; cp = cq; axisPoint = q; }
            }
            const Vec3 delta = cp.point - axisPoint;
            if (!(dist < radius)) continue;

            // At the intact-surface clamp the probe centre can lie exactly on the
            // triangle. Use the inward-facing contact direction in that degenerate
            // case instead of discarding the contact because its distance is zero.
            Vec3 direction;
            if (dist > 1e-10) direction = delta * (1.0 / dist);
            else
            {
                Vec3 normal = cross(b - a, c - a);
                const double nl = normal.length();
                if (!(nl > 1e-14)) continue;
                direction = normal * (-1.0 / nl); // force applied to tissue, into its volume
            }

            // With area stiffness, integrate a penalty pressure over this triangle.
            // Otherwise preserve the legacy N/m spring scale as a total triangle spring.
            const double k = (m_cfg.contactStiffnessPerArea > 0.0)
                           ? m_cfg.contactStiffnessPerArea * area
                           : m_cfg.contactStiffness;
            double magnitude = k * (radius - dist);
            if (m_cfg.contactMaxForce > 0.0)
                magnitude = (std::min)(magnitude, m_cfg.contactMaxForce);
            const Vec3 force = direction * magnitude;

            const double w[3] = {cp.w0, cp.w1, cp.w2};
            for (int kVertex = 0; kVertex < 3; ++kVertex)
                if (!m_dyn.isFixed(v[kVertex]) && w[kVertex] > 0.0)
                    m_dyn.addExternalForce(v[kVertex], force * w[kVertex]);

            reaction -= force;
            ++contacts;
        }

        const double forceLength = reaction.length();
        if (m_cfg.contactMaxForce > 0.0 && forceLength > m_cfg.contactMaxForce)
            reaction = reaction * (m_cfg.contactMaxForce / forceLength);
        m_numContacts = contacts;
        return reaction;
    }

    // Progressive drilling. The drill allocates its available energy across nearby
    // tetrahedra and removes at most one tetrahedron per update. A tet is removed only
    // after its accumulated work reaches energyDensity * tetVolume.
    void drill(double wallDt)
    {
        if (!(wallDt > 0.0) || !(m_cfg.cuttingPower > 0.0)) return;
        const std::vector<Vec3>& P = m_dyn.positions();
        const TensorMassModel& model = m_dyn.model();
        std::vector<int> candidates;
        for (size_t t = 0; t < model.numTets(); ++t)
        {
            if (!model.isActive((int)t)) continue;
            const TensorMassModel::Tet& T = model.tets()[t];
            const Vec3 centroid = (P[(size_t)T[0]] + P[(size_t)T[1]] +
                                   P[(size_t)T[2]] + P[(size_t)T[3]]) * 0.25;
            if ((centroid - m_probe).length() <= m_cfg.drillRadius)
                candidates.push_back((int)t);
        }
        if (candidates.empty()) return;

        const double workPerCandidate = m_cfg.cuttingPower * wallDt /
                                        (double)candidates.size();
        int ready = -1;
        double greatestProgress = -1.0;
        for (size_t i = 0; i < candidates.size(); ++i)
        {
            const int t = candidates[i];
            const double volume = model.tetVolume(t);
            if (!(volume > 0.0) || !std::isfinite(volume)) continue;
            m_cutWork[(size_t)t] += workPerCandidate;
            const double threshold = m_cfg.cuttingEnergyPerVolume * volume;
            const double progress = m_cutWork[(size_t)t] / threshold;
            if (progress >= 1.0 && progress > greatestProgress)
            {
                ready = t;
                greatestProgress = progress;
            }
        }
        if (ready < 0) return;

        {
            const TensorMassModel::Tet& T = model.tets()[(size_t)ready];
            const double cx = 0.25 * (P[(size_t)T[0]].x + P[(size_t)T[1]].x + P[(size_t)T[2]].x + P[(size_t)T[3]].x);
            m_deepestClearedX = (std::min)(m_deepestClearedX, cx);
        }
        m_dyn.cutTetrahedron(ready);
        m_cutWork[(size_t)ready] = 0.0;
        removeFloatingFragments();
        rebuildSurface();
        recomputeTimeStep(true);
        m_surfaceChanged = true;
    }

    // material that is no longer connected to the clamped base would drift freely: remove it
    void removeFloatingFragments()
    {
        const TensorMassModel& m = m_dyn.model();
        const size_t nT = m.numTets();
        std::vector<std::vector<int> > vertTets(m.numVertices());
        for (size_t t = 0; t < nT; t++)
            if (m.isActive((int)t)) for (int k = 0; k < 4; k++) vertTets[m.tets()[t][k]].push_back((int)t);

        std::vector<char> seen(nT, 0);
        std::vector<int> doomed;
        for (size_t s = 0; s < nT; s++)
        {
            if (!m.isActive((int)s) || seen[s]) continue;
            std::vector<int> comp;
            bool anchored = false;
            std::queue<int> q;
            q.push((int)s);
            seen[s] = 1;
            while (!q.empty())
            {
                const int t = q.front(); q.pop();
                comp.push_back(t);
                for (int k = 0; k < 4; k++)
                {
                    const int v = m.tets()[t][k];
                    if (m_dyn.isFixed(v)) anchored = true;
                    for (size_t u = 0; u < vertTets[v].size(); u++)
                    {
                        const int t2 = vertTets[v][u];
                        if (!seen[t2]) { seen[t2] = 1; q.push(t2); }
                    }
                }
            }
            if (!anchored) doomed.insert(doomed.end(), comp.begin(), comp.end());
        }
        for (size_t i = 0; i < doomed.size(); i++) m_dyn.cutTetrahedron(doomed[i]);
    }

    void rebuildSurface()
    {
        m_surface = extractBoundaryTriangles(m_dyn.model());
        std::vector<char> on(m_dyn.model().numVertices(), 0);
        for (size_t i = 0; i < m_surface.size(); i++) { on[m_surface[i].a] = on[m_surface[i].b] = on[m_surface[i].c] = 1; }
        m_exposed.clear();
        for (size_t v = 0; v < on.size(); v++) if (on[v]) m_exposed.push_back((int)v);
        m_nearValid = false;
        buildSmoothOffsets();
    }

    // Surface smoothing. The tetrahedral surface of a grid-cut mesh is a staircase (+- half a cell), which is
    // rougher than the sphere that touches it and makes the contact jitter. Each exposed vertex gets a fixed
    // offset = (Taubin-smoothed rest position) - (rest position); contact and display use position + offset.
    // Taubin smoothing (lambda | mu, Taubin 1995) removes the staircase without shrinking the organ.
    void buildSmoothOffsets()
    {
        const size_t nv = m_dyn.model().numVertices();
        m_smoothOffset.assign(nv, Vec3());
        if (!m_cfg.proxyContact || m_cfg.rigid || m_surface.empty()) return;
        const std::vector<Vec3>& R0 = m_dyn.model().restPositions();
        std::vector<std::vector<int> > nb(nv);
        for (size_t i = 0; i < m_surface.size(); ++i)
        {
            const int v[3] = {m_surface[i].a, m_surface[i].b, m_surface[i].c};
            for (int k = 0; k < 3; ++k) { nb[(size_t)v[k]].push_back(v[(k + 1) % 3]); nb[(size_t)v[k]].push_back(v[(k + 2) % 3]); }
        }
        for (size_t q = 0; q < m_exposed.size(); ++q)
        {
            std::vector<int>& l = nb[(size_t)m_exposed[q]];
            std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end());
        }
        {   double sum = 0.0; size_t cnt = 0;
            for (size_t i = 0; i < m_surface.size(); ++i)
            {
                const int v[3] = {m_surface[i].a, m_surface[i].b, m_surface[i].c};
                for (int k = 0; k < 3; ++k) { sum += (R0[(size_t)v[k]] - R0[(size_t)v[(k + 1) % 3]]).length(); ++cnt; }
            }
            m_meanEdge = cnt ? sum / (double)cnt : 0.0;
        }
        std::vector<Vec3> X(R0.begin(), R0.end()), Y(R0.begin(), R0.end());
        const double lam = 0.5, mu = -0.53;
        for (int it = 0; it < 2 * m_cfg.smoothingPasses; ++it)
        {
            const double f = (it % 2 == 0) ? lam : mu;
            for (size_t q = 0; q < m_exposed.size(); ++q)
            {
                const int v = m_exposed[q];
                const std::vector<int>& l = nb[(size_t)v];
                if (l.empty() || m_dyn.isFixed(v)) { Y[(size_t)v] = X[(size_t)v]; continue; }
                Vec3 avg;
                for (size_t k = 0; k < l.size(); ++k) avg += X[(size_t)l[k]];
                avg = avg * (1.0 / (double)l.size());
                Y[(size_t)v] = X[(size_t)v] + (avg - X[(size_t)v]) * f;
            }
            X.swap(Y);
        }
        for (size_t q = 0; q < m_exposed.size(); ++q)
        {
            const size_t v = (size_t)m_exposed[q];
            m_smoothOffset[v] = X[v] - R0[v];
        }
    }

    SimConfig m_cfg;
    TensorMassDynamics m_dyn;
    std::vector<SurfaceTri> m_surface;
    std::vector<int> m_exposed;
    struct ContactSample { int tri; TriangleClosestPoint cp; double w; Vec3 n; };
    Vec3 m_probe, m_lastReaction, m_hand, m_nearCenter;
    std::vector<int> m_near;
    std::vector<ContactSample> m_contactsTmp;
    std::vector<Vec3> m_smoothOffset, m_displayPos;
    std::vector<int> m_nearVerts;
    int m_lastTri = -1; double m_lastBary[3] = {0, 0, 0}; Vec3 m_lastN; double m_loadLimit = -1.0;
    Vec3 m_goodProbe; bool m_goodValid = false, m_proxyRecheck = false; int m_recoveries = 0;
    std::vector<Vec3> m_prevNormals;       // outward normals of the contacts of the previous step
    std::vector<std::pair<int, double> > m_loadW;
    double m_meanEdge = 0.0;
    bool m_nearValid = false, m_inContact = false;
    double m_rawSpringForce = 0.0;
    bool m_probeInit = false, m_surfaceChanged = true;
    int m_cutsSinceFullEstimate = 0;
    double m_dt = 1e-4, m_dtMesh = 0.0, m_dtContact = 0.0, m_acc = 0.0;
    size_t m_numTets0 = 0;
    int m_numContacts = 0, m_resetCount = 0;
    std::vector<Vec3> m_intactTris;       // rest positions of the intact outer surface triangles (3 per triangle)
    Vec3 m_boundsMin, m_boundsMax;
    double m_deepestDrillX = 1e9;
    double m_deepestClearedX = 1e9;       // x of the deepest tetrahedron the burr has removed (rigid bone)
    double m_simTime = 0.0;
    long m_logRows = 0;
    std::ofstream m_log;

    void openForceLog()
    {
        if (m_log.is_open()) m_log.close();
        if (m_cfg.forceLogPath.empty()) return;
        const bool fresh = !std::ifstream(m_cfg.forceLogPath).good();
        m_log.open(m_cfg.forceLogPath, std::ios::app);
        if (m_log && fresh) m_log << "time_s,instrument,force_N,fx_N,fy_N,fz_N,tip_x_m,tip_y_m,tip_z_m,contacts,active_tets\n";
    }
    void logForce(bool drilling)
    {
        if (!m_log.is_open()) return;
        if ((m_logRows++ % 3) != 0) return;                          // about 20 rows per second
        const Vec3& f = m_lastReaction;
        m_log << m_simTime << ',' << (drilling ? "drill" : "probe") << ',' << f.length() << ',' << f.x << ',' << f.y << ',' << f.z << ','
              << m_probe.x << ',' << m_probe.y << ',' << m_probe.z << ',' << m_numContacts << ',' << m_dyn.model().numActiveTets() << '\n';
        if ((m_logRows % 60) == 0) m_log.flush();
    }
    std::vector<double> m_cutWork;
    double m_lastMaxDisplacement = 0.0;         // deepest probe x reached while drilling
};

} // namespace tmass
