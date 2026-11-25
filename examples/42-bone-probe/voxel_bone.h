//==============================================================================
// voxel_bone.h  -  voxel model of bone with a spherical burr (standalone, no CHAI3D).
//
// Bone is a regular grid of density values (0 = air, 1 = solid bone), the way the
// published temporal-bone drilling simulators represent it:
//   Agus M, Giachetti A, Gobbetti E, Zanetti G, Zorcolo A (2003) "Real-time haptic and visual
//     simulation of bone dissection", Presence 12(1):110-122
//   Morris D, Sewell C, Barbagli F, Salisbury K, Blevins NH, Girod S (2006) "Visuohaptic
//     simulation of bone surgery for training and evaluation", IEEE CG&A 26(6):48-57.
//
// The burr is a sphere that follows the target (mouse / haptic device) through a virtual
// spring. Inside bone the burr cannot move freely: it advances only as fast as it can cut.
//
// CUTTING FORCE LAW (calibrated to measured data, see below):
//   F = a * depth^b * speed^c        [N, depth in mm, speed in mm/s]
// with a = 0.917, b = 1.22, c = 0.32: a least-squares fit (done by us, log-linear, R^2 = 0.965,
// mean error 0.10 N) of the cortical-bone means in Table A2 of
//   Dillon NP, Kratchman LB, Dietrich MS, Labadie RF, Webster RJ, Withrow TJ (2013)
//   "An experimental evaluation of the force requirements for robotic mastoidectomy",
//   Otol Neurotol 34(7):e93-e102.
// Conditions of that data: 5 mm fluted burr (radius 2.5 mm), 80,000 rpm, 20 degree angle, depths
// 0.62-1.64 mm, speeds 1-8 mm/s, formalin-fixed human temporal bone, robot-guided. Using the law
// for a spherical burr, other angles and hand-held feed is an EXTRAPOLATION.
//
// What the model does each haptic sub-step (1 kHz):
//   free space  : the burr follows the target, force = 0.
//   in contact  : the burr is pulled towards the target by the spring (force F, limited to maxForce,
//                 the force a haptic device can render). It moves at the speed at which the cutting
//                 law gives F, and removes volume = (engaged cross-section) x (speed) x dt, as in Dillon.
//   Net effect  : the force the user feels is the cutting force, it grows with depth and speed, and it
//                 relaxes to zero when the user stops pushing.
//==============================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace vbone {

constexpr double kPi = 3.14159265358979323846;

struct V3
{
    double x = 0, y = 0, z = 0;
    V3() {}
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    V3 operator+(const V3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    V3 operator-(const V3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    V3 operator*(double s) const { return {x * s, y * s, z * s}; }
    double length() const { return std::sqrt(x * x + y * y + z * z); }
};
inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

struct BurrConfig
{
    double burrRadius        = 0.0025;   // [m] (5 mm burr, as in the force data)
    double couplingStiffness = 800.0;    // [N/m] spring between target and burr
    double maxForce          = 3.3;      // [N] force limit (haptic-device peak)
    double freeSpeedMax      = 0.30;     // [m/s] burr speed in free space
    double maxCutSpeed       = 0.008;    // [m/s] 8 mm/s = upper end of the measured speed range
    double minDepthMM        = 0.02;     // [mm] depth floor inside the force law
    // force law  F = a * depth^b * speed^c   (N, mm, mm/s)
    double lawA = 0.917, lawDepthExp = 1.22, lawSpeedExp = 0.32;
    double depthGain         = 1.0;      // calibration of the engaged depth (see test_voxel_bone)
    double forceScaleBone    = 1.0;      // multiplies the law (1 = cortical data; mastoid trabecular ~ 0.5-0.65)
    double solidThreshold    = 0.5;      // a voxel counts as bone above this density
    double substep           = 0.001;    // [s] physics step ( 1 kHz )
    std::string logPath;                 // CSV log of force / depth / speed (empty = off)
    // ---- probe-only mode: the bone is rigid and is NOT removed; the sphere is a touch probe ----
    bool   probeOnly         = false;
    double boneE             = 10.0e9;   // [Pa]  PLACEHOLDER - replace with a cited value
    double boneNu            = 0.30;     //       PLACEHOLDER - replace with a cited value
    double toolE             = 210.0e9;  // [Pa]  steel tip (handbook value)
    double toolNu            = 0.30;
    double globalCompliance  = 0.0;      // [m/N] whole-bone compliance from the static FE solve (0 = not yet available)
    double contactSmoothing  = 1.0;      // 1 = blur the contact copy of the density once ([1 2 1]^3) to hide the voxel staircase
};

class VoxelBone
{
public:
    // density: nx*ny*nz values in [0,1], x fastest. origin = corner of voxel (0,0,0), h = voxel size [m]
    void init(int nx, int ny, int nz, double h, const V3& origin, const std::vector<float>& density, const BurrConfig& cfg)
    {
        if (nx < 2 || ny < 2 || nz < 2 || !(h > 0.0)) throw std::invalid_argument("VoxelBone: bad grid");
        if (density.size() != (size_t)nx * ny * nz) throw std::invalid_argument("VoxelBone: density size mismatch");
        m_nx = nx; m_ny = ny; m_nz = nz; m_h = h; m_origin = origin; m_cfg = cfg;
        m_rho0 = density; m_rho = density;
        m_rhoC.clear();
        if (m_cfg.probeOnly) buildContactField();
        m_p = V3(); m_started = false;
        m_removed = 0.0; m_time = 0.0; m_changedBox = false;
        m_lastForce = V3(); m_lastDepthMM = 0; m_lastSpeed = 0; m_contact = false; m_logRows = 0;
        if (m_log.is_open()) m_log.close();
        if (!m_cfg.logPath.empty())
        {
            m_log.open(m_cfg.logPath);
            if (m_log) m_log << "time_s,force_N,fx,fy,fz,depth_mm,speed_mm_s,removed_mm3,burr_x,burr_y,burr_z,in_contact,raw_force_N,hertz_um,contact_deformation_um\n";
        }
    }
    void reset()
    {
        m_rho = m_rho0; m_removed = 0.0; m_started = false; m_lastForce = V3(); m_contact = false; m_changedBox = true;
    }

    // advance by wallDt seconds towards `target`; returns the force to render on the user [N]
    V3 update(double wallDt, const V3& target)
    {
        if (!m_started) { m_p = target; m_started = true; }
        int n = (int)std::ceil(wallDt / m_cfg.substep);
        n = (std::max)(1, (std::min)(n, 50));
        const double dt = wallDt / n;
        V3 fSum; int nContact = 0;
        for (int s = 0; s < n; ++s)
        {
            V3 f = substep(dt, target);
            fSum = fSum + f;
            if (m_contact) ++nContact;
        }
        m_lastForce = fSum * (1.0 / n);
        m_time += wallDt;
        logRow();
        return m_lastForce;
    }

    // ---- results / state -----------------------------------------------------------
    const V3& burr() const { return m_p; }
    const V3& lastForce() const { return m_lastForce; }
    double lastDepthMM() const { return m_lastDepthMM; }      // drill mode: engaged depth; probe mode: hand penetration beyond the surface
    double rawForceN() const { return m_rawForce; }           // spring force before the device limit
    double hertzIndentationUM() const { return m_hertz * 1e6; }
    double contactDeformationUM() const { return m_dispUM; }  // Hertz indentation + whole-bone compliance x force
    double lastSpeedMMs() const { return m_lastSpeed * 1e3; }
    bool inContact() const { return m_contact; }
    double removedVolumeMM3() const { return m_removed * 1e9; }
    int nx() const { return m_nx; } int ny() const { return m_ny; } int nz() const { return m_nz; }
    double voxelSize() const { return m_h; }
    const V3& origin() const { return m_origin; }
    const std::vector<float>& density() const { return m_rho; }
    double solidVolumeMM3() const { double v = 0; for (float r : m_rho) v += r; return v * m_h * m_h * m_h * 1e9; }
    BurrConfig& config() { return m_cfg; }

    // cutting-force law and its inverse (N, mm, mm/s)
    double lawForce(double depthMM, double speedMMs) const
    {
        depthMM = (std::max)(depthMM, m_cfg.minDepthMM);
        return m_cfg.forceScaleBone * m_cfg.lawA * std::pow(depthMM, m_cfg.lawDepthExp) * std::pow((std::max)(speedMMs, 0.0), m_cfg.lawSpeedExp);
    }
    double lawSpeedMMs(double forceN, double depthMM) const
    {
        depthMM = (std::max)(depthMM, m_cfg.minDepthMM);
        const double base = m_cfg.forceScaleBone * m_cfg.lawA * std::pow(depthMM, m_cfg.lawDepthExp);
        if (!(forceN > 0.0)) return 0.0;
        return std::pow(forceN / base, 1.0 / m_cfg.lawSpeedExp);
    }

    // area of the circular segment of the burr engaged to depth d (the cross-section of the cut)
    static double segmentArea(double R, double d)
    {
        d = (std::min)((std::max)(d, 0.0), 2.0 * R);
        if (d <= 0.0) return 0.0;
        const double a = R - d;
        return R * R * std::acos(a / R) - a * std::sqrt((std::max)(2.0 * R * d - d * d, 0.0));
    }

private:
    struct Contact { double overlapVolume = 0.0; double depth = 0.0; V3 normal; };

    // sphere/bone overlap -> equivalent depth of a spherical cap with the same volume
    Contact evaluate(const V3& p) const { return evaluateIn(p, m_rho, m_cfg.solidThreshold); }
    Contact evaluateIn(const V3& p, const std::vector<float>& field, double normalThreshold) const
    {
        Contact c;
        const double R = m_cfg.burrRadius, R2 = R * R, h = m_h, h3 = h * h * h;
        const int i0 = (std::max)(0, (int)std::floor((p.x - R - m_origin.x) / h)), i1 = (std::min)(m_nx - 1, (int)std::floor((p.x + R - m_origin.x) / h));
        const int j0 = (std::max)(0, (int)std::floor((p.y - R - m_origin.y) / h)), j1 = (std::min)(m_ny - 1, (int)std::floor((p.y + R - m_origin.y) / h));
        const int k0 = (std::max)(0, (int)std::floor((p.z - R - m_origin.z) / h)), k1 = (std::min)(m_nz - 1, (int)std::floor((p.z + R - m_origin.z) / h));
        V3 nsum;
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i)
                {
                    const float rho = field[idx(i, j, k)];
                    if (rho <= 0.0f) continue;
                    const V3 x(m_origin.x + (i + 0.5) * h, m_origin.y + (j + 0.5) * h, m_origin.z + (k + 0.5) * h);
                    const V3 d = p - x;
                    const double r2 = dot(d, d);
                    if (r2 >= R2) continue;
                    c.overlapVolume += rho * h3;
                    if (rho > normalThreshold && r2 > 1e-18)
                    {
                        const double r = std::sqrt(r2);
                        nsum = nsum + d * ((R - r) * rho / r);
                    }
                }
        const double nl = nsum.length();
        if (nl > 1e-18) c.normal = nsum * (1.0 / nl);
        // invert  V(d) = pi d^2 (3R - d) / 3
        const double V = c.overlapVolume;
        if (V > 0.0)
        {
            double lo = 0.0, hi = 2.0 * R;
            for (int it = 0; it < 40; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                const double vm = kPi * mid * mid * (3.0 * R - mid) / 3.0;
                (vm < V) ? lo = mid : hi = mid;
            }
            c.depth = 0.5 * (lo + hi);
        }
        return c;
    }

    V3 substep(double dt, const V3& target)
    {
        if (m_cfg.probeOnly) return probeSubstep(dt, target);
        const double R = m_cfg.burrRadius;
        const V3 gap = target - m_p;
        const double gl = gap.length();
        const V3 dir = (gl > 1e-12) ? gap * (1.0 / gl) : V3();
        // 1. try to move towards the target (limited speed)
        const double stepFree = (std::min)(gl, m_cfg.freeSpeedMax * dt);
        const V3 cand = m_p + dir * stepFree;
        const Contact cc = evaluate(cand);
        const double capMin = 1e-3 * (kPi * R * R * R * 4.0 / 3.0);     // 0.1 % of the sphere volume
        if (cc.overlapVolume < capMin)
        {
            m_p = cand;
            m_contact = false; m_lastDepthMM = 0.0; m_lastSpeed = 0.0;
            return V3();
        }
        // 2. in contact: the burr moves only at the speed the cutting law allows for the applied force
        m_contact = true;
        double F = m_cfg.couplingStiffness * gl;
        if (F > m_cfg.maxForce) F = m_cfg.maxForce;
        const Contact c0 = evaluate(m_p);
        const double depthMM = m_cfg.depthGain * (std::max)(c0.depth, cc.depth) * 1e3;
        double speed = lawSpeedMMs(F, depthMM) * 1e-3;                        // [m/s]
        speed = (std::min)(speed, m_cfg.maxCutSpeed);
        speed = (std::min)(speed, gl / dt);
        const V3 move = dir * (speed * dt);
        m_p = m_p + move;
        m_lastDepthMM = depthMM;
        m_lastSpeed = speed;
        // 3. remove material:  volume = engaged cross-section x speed x dt
        const double dRealM = (std::max)(c0.depth, cc.depth);
        const double dV = segmentArea(R, dRealM) * speed * dt;
        erode(m_p, dV);
        // force on the user: spring pulling back towards the burr, limited to maxForce
        return dir * (-F);
    }

    // ---- probe-only: rigid frictionless constraint. The sphere cannot enter the bone; the hand (target)
    //      can, through the virtual spring, and the user feels the spring force along the surface normal.
    void buildContactField()
    {
        m_rhoC = m_rho;
        if (m_cfg.contactSmoothing <= 0.0) return;
        const int nx = m_nx, ny = m_ny, nz = m_nz;
        std::vector<float> tmp(m_rhoC.size());
        const float w[3] = {0.25f, 0.5f, 0.25f};
        for (int axis = 0; axis < 3; ++axis)
        {
            for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
            {
                float a = 0.f;
                for (int t = -1; t <= 1; ++t)
                {
                    int ii = i + (axis == 0 ? t : 0), jj = j + (axis == 1 ? t : 0), kk = k + (axis == 2 ? t : 0);
                    ii = (std::max)(0, (std::min)(nx - 1, ii)); jj = (std::max)(0, (std::min)(ny - 1, jj)); kk = (std::max)(0, (std::min)(nz - 1, kk));
                    a += w[t + 1] * m_rhoC[idx(ii, jj, kk)];
                }
                tmp[idx(i, j, k)] = a;
            }
            m_rhoC.swap(tmp);
        }
        // remove the faint halo outside the surface: iso-level 0.5 stays where it was
        for (float& v : m_rhoC) { v = (v - 0.35f) / 0.3f; v = v < 0.f ? 0.f : (v > 1.f ? 1.f : v); }
    }

    double hertzIndentation(double F) const
    {
        if (!(F > 0.0)) return 0.0;
        const double Es = 1.0 / ((1.0 - m_cfg.boneNu * m_cfg.boneNu) / m_cfg.boneE + (1.0 - m_cfg.toolNu * m_cfg.toolNu) / m_cfg.toolE);
        return std::pow(3.0 * F / (4.0 * Es * std::sqrt(m_cfg.burrRadius)), 2.0 / 3.0);
    }

    V3 probeSubstep(double dt, const V3& target)
    {
        const double R = m_cfg.burrRadius;
        const V3 gap = target - m_p;
        const double gl = gap.length();
        const V3 dir = (gl > 1e-12) ? gap * (1.0 / gl) : V3();
        V3 cand = m_p + dir * (std::min)(gl, m_cfg.freeSpeedMax * dt);
        const double capMin = 1e-3 * (kPi * R * R * R * 4.0 / 3.0);
        Contact c = evaluateIn(cand, m_rhoC, 0.0);
        if (c.overlapVolume < capMin)
        {
            m_p = cand; m_contact = false;
            m_lastDepthMM = 0; m_lastSpeed = 0; m_rawForce = 0; m_hertz = 0; m_dispUM = 0;
            return V3();
        }
        m_contact = true;
        if (c.normal.length() > 0.5) m_lastN = c.normal;
        const V3 n = m_lastN;                                     // outward from the bone (last good one if degenerate)
        // push the sphere back out of the bone along the normal (a few iterations), never more than R at once
        V3 q = cand; Contact cq = c;
        for (int it = 0; it < 6 && cq.overlapVolume >= capMin; ++it)
        {
            const double push = (std::min)((std::max)(cq.depth, 0.25 * m_h), R);
            const V3 nn = (cq.normal.length() > 0.5) ? cq.normal : n;
            q = q + nn * push;
            cq = evaluateIn(q, m_rhoC, 0.0);
        }
        if (cq.overlapVolume < capMin) m_p = q;                   // else: stay where we were
        // spring force along the normal, from the hand's penetration beyond the sphere
        const V3 gap2 = target - m_p;
        const double pen = (std::max)(0.0, -dot(gap2, n));        // how far the hand is inside the bone [m]
        const double raw = m_cfg.couplingStiffness * pen;
        double F = raw; if (F > m_cfg.maxForce) F = m_cfg.maxForce;
        m_rawForce = raw;
        m_lastDepthMM = pen * 1e3;                                // = hand penetration (the sphere itself stays on the surface)
        m_lastSpeed = 0;
        m_hertz = hertzIndentation(F);
        m_dispUM = (m_hertz + m_cfg.globalCompliance * F) * 1e6;
        return n * F;
    }

    void erode(const V3& p, double dV)
    {
        if (!(dV > 0.0)) return;
        const double R = m_cfg.burrRadius, R2 = R * R, h = m_h, h3 = h * h * h;
        const int i0 = (std::max)(0, (int)std::floor((p.x - R - m_origin.x) / h)), i1 = (std::min)(m_nx - 1, (int)std::floor((p.x + R - m_origin.x) / h));
        const int j0 = (std::max)(0, (int)std::floor((p.y - R - m_origin.y) / h)), j1 = (std::min)(m_ny - 1, (int)std::floor((p.y + R - m_origin.y) / h));
        const int k0 = (std::max)(0, (int)std::floor((p.z - R - m_origin.z) / h)), k1 = (std::min)(m_nz - 1, (int)std::floor((p.z + R - m_origin.z) / h));
        // weights: voxels deeper inside the burr give way first
        double wsum = 0.0;
        for (int pass = 0; pass < 2; ++pass)
        {
            for (int k = k0; k <= k1; ++k)
                for (int j = j0; j <= j1; ++j)
                    for (int i = i0; i <= i1; ++i)
                    {
                        float& rho = m_rho[idx(i, j, k)];
                        if (rho <= 0.0f) continue;
                        const V3 x(m_origin.x + (i + 0.5) * m_h, m_origin.y + (j + 0.5) * m_h, m_origin.z + (k + 0.5) * m_h);
                        const V3 d = p - x; const double r2 = dot(d, d);
                        if (r2 >= R2) continue;
                        const double w = (R - std::sqrt(r2)) * rho;
                        if (pass == 0) wsum += w;
                        else
                        {
                            const double take = (std::min)((double)rho, dV * w / wsum / h3);
                            rho -= (float)take; if (rho < 1e-4f) rho = 0.0f;
                            m_removed += take * h3;
                        }
                    }
            if (!(wsum > 0.0)) return;
        }
    }

    void logRow()
    {
        if (!m_log.is_open()) return;
        if ((m_logRows++ % 5) != 0) return;
        m_log << m_time << ',' << m_lastForce.length() << ',' << m_lastForce.x << ',' << m_lastForce.y << ',' << m_lastForce.z << ','
              << m_lastDepthMM << ',' << m_lastSpeed * 1e3 << ',' << m_removed * 1e9 << ',' << m_p.x << ',' << m_p.y << ',' << m_p.z << ','
              << (m_contact ? 1 : 0) << ',' << m_rawForce << ',' << m_hertz * 1e6 << ',' << m_dispUM << '\n';
        if ((m_logRows % 100) == 0) m_log.flush();
    }

    size_t idx(int i, int j, int k) const { return ((size_t)k * m_ny + j) * m_nx + i; }

    int m_nx = 0, m_ny = 0, m_nz = 0;
    double m_h = 0.0;
    V3 m_origin;
    BurrConfig m_cfg;
    std::vector<float> m_rho0, m_rho, m_rhoC;
    V3 m_lastN = V3(0, 0, 1);
    double m_rawForce = 0.0, m_hertz = 0.0, m_dispUM = 0.0;
    V3 m_p, m_lastForce;
    bool m_started = false, m_contact = false, m_changedBox = false;
    double m_removed = 0.0, m_time = 0.0, m_lastDepthMM = 0.0, m_lastSpeed = 0.0;
    long m_logRows = 0;
    std::ofstream m_log;
};

} // namespace vbone
