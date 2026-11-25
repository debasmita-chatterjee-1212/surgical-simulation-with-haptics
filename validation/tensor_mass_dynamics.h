//==============================================================================
// tensor_mass_dynamics.h
//
// Step 2 of the surgical simulator: explicit Newtonian integration of the
// tensor-mass model (thesis section 4.4.3), on top of tensor_mass.h.
//
// Equation of motion (thesis):      m_i P'' = gamma_i P' + F_i
// Explicit central-difference scheme (thesis):
//     (m/dt^2 - gamma/(2dt)) P(t+1) = F + (2m/dt^2) P(t) - (m/dt^2 + gamma/(2dt)) P(t-1)
//
// SIGN CONVENTION: with gamma as written in the thesis, damping needs gamma < 0.
// Here the damping coefficient is c = -gamma >= 0, stored per unit mass:
// c_i = dampingPerMass * m_i.   Then
//     A = m/dt^2 + c/(2dt),   B = m/dt^2 - c/(2dt)
//     P(t+1) = ( F + (2m/dt^2) P(t) - B P(t-1) ) / A
//
// Explicit integration is only conditionally stable: dt must stay below
// 2/omega_max, where omega_max^2 is the largest eigenvalue of M^-1 K.
//==============================================================================
#pragma once

#include "tensor_mass.h"
#include <limits>
#include <random>
#include <stdexcept>

namespace tmass {

class TensorMassDynamics
{
public:
    // density [kg/m^3], dampingPerMass [1/s], dt [s]
    void init(const TensorMassModel& model, double density, double dampingPerMass, double dt)
    {
        if (!(density > 0.0) || !std::isfinite(density)) throw std::invalid_argument("density must be finite and positive");
        if (dampingPerMass < 0.0 || !std::isfinite(dampingPerMass)) throw std::invalid_argument("damping must be finite and non-negative");
        validateTimeStep(dt);
        m_model = model;
        m_rho = density;
        m_damp = dampingPerMass;
        m_dt = dt;
        const size_t n = m_model.numVertices();
        m_P = m_Pprev = m_model.restPositions();
        m_F.assign(n, Vec3());
        m_Pnext.assign(n, Vec3());
        m_powerVec.clear();
        m_fixed.assign(n, 0);
        m_ext.assign(n, Vec3());
        recomputeMasses();
    }

    void setTimeStep(double dt) { validateTimeStep(dt); m_dt = dt; }
    void setDampingPerMass(double d)
    {
        if (d < 0.0 || !std::isfinite(d)) throw std::invalid_argument("damping must be finite and non-negative");
        m_damp = d;
    }

    // pin a vertex at its rest position (boundary condition)
    void fixVertex(int i)
    {
        checkVertex(i);
        m_fixed[(size_t)i] = 1;
        m_P[(size_t)i] = m_Pprev[(size_t)i] = m_model.restPositions()[(size_t)i];
    }
    void setExternalForce(int i, const Vec3& f) { checkVertex(i); m_ext[(size_t)i] = f; }
    void addExternalForce(int i, const Vec3& f) { checkVertex(i); m_ext[(size_t)i] += f; }
    void clearExternalForces() { std::fill(m_ext.begin(), m_ext.end(), Vec3()); }

    // set an initial displacement field (zero initial velocity)
    void setDisplacement(const std::vector<Vec3>& u)
    {
        for (size_t i = 0; i < u.size() && i < m_P.size(); i++)
        {
            if (m_fixed[i]) continue;
            m_P[i] = m_Pprev[i] = m_model.restPositions()[i] + u[i];
        }
    }

    // number of threads used by step() and the time-step estimate (1 = serial).
    // Small meshes always run serially: the hand-over costs more than it saves.
    void setThreads(int n)
    {
        m_threads = (std::max)(1, n);
        if (m_threads > 1) { if (!m_pool) m_pool = std::make_shared<WorkerPool>(m_threads); else m_pool->setThreads(m_threads); }
        else m_pool.reset();
    }
    int threads() const { return m_threads; }

    // one explicit time step. Every vertex gathers its own force and writes only its
    // own new position, so the loop can be split over threads without any locking.
    void step()
    {
        const size_t n = m_P.size();
        if (m_Pnext.size() != n) m_Pnext.assign(n, Vec3());
        const double inv_dt2 = 1.0 / (m_dt * m_dt);
        const double dt = m_dt, damp = m_damp;
        const std::vector<Vec3>& rest = m_model.restPositions();
        auto body = [&](size_t b, size_t e)
        {
            for (size_t i = b; i < e; ++i)
            {
                if (m_fixed[i])       { m_Pnext[i] = rest[i]; continue; }
                if (m_mass[i] <= 0.0) { m_Pnext[i] = m_P[i]; continue; }     // detached by cutting: frozen
                const double m = m_mass[i];
                const double c = damp * m;
                const double A = m * inv_dt2 + c / (2.0 * dt);
                const double B = m * inv_dt2 - c / (2.0 * dt);
                const Vec3 F = m_model.forceAt((int)i, m_P);
                m_Pnext[i] = ((F + m_ext[i]) + m_P[i] * (2.0 * m * inv_dt2) - m_Pprev[i] * B) * (1.0 / A);
            }
        };
        runParallel(n, body);
        m_Pprev.swap(m_P);          // Pprev <- P (old Pprev is now in m_P and gets overwritten next)
        m_P.swap(m_Pnext);          // P <- new positions, Pnext <- spare buffer
    }

    // cut: remove the tetrahedron's stiffness (thesis: 10 tensor updates) and its mass
    bool cutTetrahedron(int t)
    {
        if (!m_model.removeTetrahedron(t)) return false;
        recomputeMasses();
        return true;
    }

    // Largest angular frequency of M^-1 K over the free vertices (power iteration,
    // Rayleigh quotient => converges from below, so a small safety factor is applied).
    // warmStart = true continues from the vector found last time (plus a little noise so a
    // new local mode near a cut can still be picked up); this needs far fewer iterations.
    double estimateMaxAngularFrequency(int iterations = 400, bool warmStart = false) const
    {
        const size_t n = m_P.size();
        std::vector<int> freeVertices;
        for (size_t i = 0; i < n; ++i) if (isFree((int)i)) freeVertices.push_back((int)i);
        if (freeVertices.empty()) return 0.0;
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        std::vector<Vec3>& x = m_powerVec;
        const bool warm = warmStart && x.size() == n;
        if (!warm) x.assign(n, Vec3());
        for (int i : freeVertices)
        {
            const size_t k = (size_t)i;
            if (!warm) x[k] = Vec3(U(rng), U(rng), U(rng));
            else x[k] += Vec3(U(rng), U(rng), U(rng)) * 0.02 * (std::fabs(x[k].x) + std::fabs(x[k].y) + std::fabs(x[k].z) + 1e-12);
        }
        for (size_t i = 0; i < n; ++i) if (!isFree((int)i)) x[i] = Vec3();
        std::vector<Vec3> y(n);

        double lambda = 0.0;
        const int nf = (int)freeVertices.size();
        for (int it = 0; it < iterations; ++it)
        {
            auto body = [&](size_t b, size_t e)
            {
                for (size_t q = b; q < e; ++q)
                {
                    const size_t k = (size_t)freeVertices[q];
                    y[k] = m_model.stiffnessRow((int)k, x);                    // K x
                }
            };
            runParallel((size_t)nf, body);
            double xKx = 0.0, xMx = 0.0, yMy = 0.0;
            for (int i : freeVertices)
            {
                const size_t k = (size_t)i;
                xKx += dot(x[k], y[k]);
                xMx += m_mass[k] * dot(x[k], x[k]);
                y[k] = y[k] * (1.0 / m_mass[k]);                              // M^-1 K x
                yMy += m_mass[k] * dot(y[k], y[k]);
            }
            if (!(xMx > 0.0) || !std::isfinite(xMx)) return 0.0;
            lambda = xKx / xMx;
            if (!(yMy > 0.0) || !std::isfinite(yMy)) return 0.0;
            const double invY = 1.0 / std::sqrt(yMy);
            for (int i : freeVertices) x[(size_t)i] = y[(size_t)i] * invY;
        }
        if (!(lambda > 0.0) || !std::isfinite(lambda)) return 0.0;
        return std::sqrt(lambda) * 1.02;
    }

    double estimateCriticalTimeStep(int iterations = 400, bool warmStart = false) const
    {
        const double omega = estimateMaxAngularFrequency(iterations, warmStart);
        return (omega > 0.0 && std::isfinite(omega)) ? 2.0 / omega : std::numeric_limits<double>::infinity();
    }

    double invariantEnergy() const
    {
        std::vector<Vec3> F;
        m_model.computeForces(m_Pprev, F);
        double ke = 0.0, pe = 0.0;
        for (size_t i = 0; i < m_P.size(); ++i)
        {
            const Vec3 v = (m_P[i] - m_Pprev[i]) * (1.0 / m_dt);
            ke += 0.5 * m_mass[i] * dot(v, v);
            pe -= 0.5 * dot(m_P[i] - m_model.restPositions()[i], F[i]);
        }
        return ke + pe;
    }
    double totalEnergy() const
    {
        double ke = 0.0;
        for (size_t i = 0; i < m_P.size(); ++i)
        {
            const Vec3 v = (m_P[i] - m_Pprev[i]) * (1.0 / m_dt);
            ke += 0.5 * m_mass[i] * dot(v, v);
        }
        return ke + m_model.elasticEnergy(m_P);
    }
    double maxDisplacement() const
    {
        double d = 0.0;
        for (size_t i = 0; i < m_P.size(); ++i)
        {
            const double len = (m_P[i] - m_model.restPositions()[i]).length();
            if (!std::isfinite(len)) return std::numeric_limits<double>::infinity();
            d = (std::max)(d, len);
        }
        return d;
    }
    double totalMass() const { double s = 0.0; for (double m : m_mass) s += m; return s; }

    const std::vector<Vec3>& positions() const { return m_P; }
    const std::vector<double>& masses() const { return m_mass; }
    const TensorMassModel& model() const { return m_model; }
    double timeStep() const { return m_dt; }
    bool isFixed(int i) const { checkVertex(i); return m_fixed[(size_t)i] != 0; }

private:
    static void validateTimeStep(double dt)
    {
        if (!(dt > 0.0) || !std::isfinite(dt)) throw std::invalid_argument("time step must be finite and positive");
    }
    void checkVertex(int i) const
    {
        if (i < 0 || (size_t)i >= m_P.size()) throw std::out_of_range("vertex index out of range");
    }
    bool isFree(int i) const { return !m_fixed[(size_t)i] && m_mass[(size_t)i] > 0.0; }
    void recomputeMasses()
    {
        m_mass.assign(m_model.numVertices(), 0.0);
        for (size_t t = 0; t < m_model.numTets(); ++t)
        {
            if (!m_model.isActive((int)t)) continue;
            const double mv = m_rho * m_model.tetVolume((int)t) * 0.25;
            if (!(mv > 0.0) || !std::isfinite(mv)) throw std::runtime_error("invalid tetrahedron volume while assembling lumped mass");
            for (int k = 0; k < 4; ++k) m_mass[(size_t)m_model.tets()[t][k]] += mv;
        }
    }

    TensorMassModel m_model;
    double m_rho = 1000.0, m_damp = 0.0, m_dt = 1e-3;
    // run body(begin,end) over [0,count), on several threads when it is worth it
    template <class Body>
    void runParallel(size_t count, Body& body) const
    {
        if (m_pool && m_threads > 1 && count >= 3000) m_pool->parallelFor(count, std::function<void(size_t, size_t)>(std::ref(body)));
        else body(0, count);
    }

    std::vector<Vec3> m_P, m_Pprev, m_Pnext, m_F, m_ext;
    mutable std::vector<Vec3> m_powerVec;          // last power-iteration vector (warm start)
    int m_threads = 1;
    std::shared_ptr<WorkerPool> m_pool;
    std::vector<double> m_mass;
    std::vector<char> m_fixed;
};

} // namespace tmass
