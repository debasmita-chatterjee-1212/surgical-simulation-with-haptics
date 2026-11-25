// scalpel.h - shallow surface cuts with a scalpel (no change of the tetrahedral mesh).
//
// What is modelled
//  * The blade only cuts the surface layer: the cut is a groove that is stored as a list of points ANCHORED to the surface
//    triangles (triangle + barycentric weights), so the groove moves and stretches with the tissue.
//  * Cutting force (fracture mechanics, steady cutting): F = J * d   [N]
//        J = fracture toughness [J/m^2 = N/m], d = depth of the cut [m]  (work of fracture per unit crack area, Chanthasopeephan et al.).
//    The depth follows how far the hand is pushed below the proxy, limited to maxDepth ("shallow" cuts).
//  * The tissue is loaded with the cutting force (never more than the organ's own force limit), so the surface is dragged along.
//  * Where the blade crosses an existing cut the tissue is weaker: force limit x softness (a simple heuristic, not measured).
// Not modelled: the groove does not open under tension, no bleeding, no friction on the blade.
#pragma once
#include "surgical_sim.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace tmass {

struct ScalpelParams
{
    double toughness = 200.0;     // J [N/m]. LIVER: Chanthasopeephan et al. scalpel cutting of porcine liver, 187-225 J/m^2;
                                  // Gokgol et al. 2012 (bovine liver, needle) 164 J/m^2. For the other organs this is ASSUMED (liver value).
    double maxDepth = 0.004;      // [m] deepest cut: "shallow cuts" only
    double startForce = 0.15;     // [N] pressing force at which the blade starts to cut. ASSUMED.
    double sampleSpacing = 0.0005;// [m] distance between stored groove points
    double depthTau = 0.08;       // [s] smoothing of the cut depth
    double resumeTime = 0.5;      // [s] a stroke resumes if the blade is back within this time while SPACE stays down
    double lostTime = 0.15;       // [s] contact lost for this long ends a stroke
    double softness = 0.6;        // force limit factor on an existing cut. ASSUMED heuristic.
    double nearGroove = 0.0015;   // [m] distance that counts as "on a cut"
    size_t maxSamples = 6000;     // oldest strokes are dropped beyond this
};

struct CutSample { int tri; double b[3]; double depth; };
struct CutStroke { std::vector<CutSample> s; };

class Scalpel
{
public:
    ScalpelParams params;

    void clear() { m_strokes.clear(); m_cutting = false; m_depth = 0.0; m_lost = 0.0; m_softened = false; m_length = 0.0; m_total = 0; m_keyHeld = false; m_sinceStop = 1e9; m_haveLast = false; }

    // call BEFORE sim.update(): sets the force limit for this frame
    void beforeStep(SurgicalSim& sim) const
    {
        const double cap = sim.config().contactMaxForce;
        double limit = -1.0;
        if (m_cutting) limit = (std::min)(params.toughness * (std::max)(m_depth, 0.0005), cap > 0.0 ? cap : 1e9);
        if (m_softened) limit = (limit >= 0.0 ? limit : (cap > 0.0 ? cap : 1e9)) * params.softness;
        sim.setLoadLimit(limit);
    }

    // call AFTER sim.update()
    void afterStep(SurgicalSim& sim, bool cutKey, double dt)
    {
        int tri; double b[3]; Vec3 n;
        if (!sim.contactInfo(tri, b, n))
        {
            m_lost += dt;
            m_sinceStop += dt;
            if (m_cutting && m_lost > params.lostTime) { m_cutting = false; m_sinceStop = 0.0; }
            m_depth += (0.0 - m_depth) * (1.0 - std::exp(-dt / params.depthTau));
            m_softened = false;
            return;
        }
        m_lost = 0.0;
        const Vec3 gapV = sim.probe() - sim.hand();
        const double gap = (std::max)(0.0, dot(gapV, n));                 // how far the hand is below the proxy, along the normal
        const std::vector<Vec3>& D = sim.displayPositions();
        const std::vector<SurfaceTri>& S = sim.surface();
        const SurfaceTri& t = S[(size_t)tri];
        const Vec3 pc = D[(size_t)t.a] * b[0] + D[(size_t)t.b] * b[1] + D[(size_t)t.c] * b[2];

        if (!m_cutting)
        {
            m_sinceStop += dt;
            if (!cutKey) m_keyHeld = false;
            if (cutKey && sim.rawSpringForce() >= params.startForce)
            {
                // the same stroke goes on if the blade only left the tissue for a moment while SPACE stayed down
                const bool resume = m_keyHeld && m_sinceStop < params.resumeTime && !m_strokes.empty();
                m_cutting = true; m_keyHeld = true;
                if (!resume) { m_strokes.push_back(CutStroke()); m_depth = 0.0; m_haveLast = false; }
            }
        }
        else if (!cutKey || gap < 0.0001) { m_cutting = false; m_sinceStop = 0.0; if (!cutKey) m_keyHeld = false; }

        if (m_cutting)
        {
            const double target = (std::min)(gap, params.maxDepth);
            m_depth += (target - m_depth) * (1.0 - std::exp(-dt / params.depthTau));
            if (!m_haveLast || (pc - m_lastPos).length() >= params.sampleSpacing)
            {
                CutSample c; c.tri = tri; c.b[0] = b[0]; c.b[1] = b[1]; c.b[2] = b[2]; c.depth = m_depth;
                m_strokes.back().s.push_back(c); m_total++;
                if (m_haveLast) m_length += (pc - m_lastPos).length();
                m_lastPos = pc; m_haveLast = true;
                trim();
            }
        }
        else m_depth += (0.0 - m_depth) * (1.0 - std::exp(-dt / params.depthTau));

        // on an existing cut? (ignore the newest part of the stroke being drawn)
        m_softened = false;
        const double r2 = params.nearGroove * params.nearGroove;
        for (size_t si = 0; si < m_strokes.size() && !m_softened; ++si)
        {
            const std::vector<CutSample>& v = m_strokes[si].s;
            size_t end = v.size();
            if (m_cutting && si + 1 == m_strokes.size()) end = (end > 12) ? end - 12 : 0;
            for (size_t k = 0; k < end; ++k)
            {
                const SurfaceTri& u = S[(size_t)v[k].tri];
                const Vec3 q = D[(size_t)u.a] * v[k].b[0] + D[(size_t)u.b] * v[k].b[1] + D[(size_t)u.c] * v[k].b[2];
                const Vec3 d = q - pc;
                if (dot(d, d) < r2) { m_softened = true; break; }
            }
        }
    }

    bool cutting() const { return m_cutting; }
    bool onCut() const { return m_softened; }
    double depth() const { return m_depth; }                 // [m] current depth of the cut
    double length() const { return m_length; }                // [m] total length of cut made since the last clear
    size_t numSamples() const { return m_total; }
    const std::vector<CutStroke>& strokes() const { return m_strokes; }

private:
    void trim()
    {
        while (m_total > params.maxSamples && m_strokes.size() > 1)
        {
            m_total -= m_strokes.front().s.size();
            m_strokes.erase(m_strokes.begin());
        }
    }
    std::vector<CutStroke> m_strokes;
    bool m_cutting = false, m_softened = false, m_haveLast = false, m_keyHeld = false;
    double m_sinceStop = 1e9;
    double m_depth = 0.0, m_lost = 0.0, m_length = 0.0;
    size_t m_total = 0;
    Vec3 m_lastPos;
};

} // namespace tmass
