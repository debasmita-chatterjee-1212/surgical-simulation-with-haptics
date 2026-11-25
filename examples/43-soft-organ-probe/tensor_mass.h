//==============================================================================
// tensor_mass.h
//
// Step 1 of the surgical simulator: tetrahedral mesh + precomputed stiffness
// tensors for the linear-elastic "tensor-mass" model (thesis section 4.4.3,
// after Cotin, Delingette & Ayache 2000).
//
// Header-only, no dependencies (no CHAI3D, no Eigen), so it can be tested on
// its own and later dropped into the CHAI3D application unchanged.
//
// Conventions
//   - P0[i]  : rest position of vertex i
//   - u[i]   : displacement P[i] - P0[i]
//   - For tetrahedron T with vertices (a,b,c,d) and constant shape-function
//     gradients b_k, the 3x3 stiffness block between local vertices k and l is
//         K_kl = V * ( lambda * b_k b_l^T + mu * b_l b_k^T + mu * (b_k.b_l) I )
//     which is the same K^T_jk of the thesis, written out explicitly.
//   - Blocks are stored per VERTEX ([Kii]) and per EDGE ([Kij], one tensor per
//     edge, since [Kji] = [Kij]^T), exactly as in thesis Fig. 4.7.
//   - Restoring force on vertex i:  F_i = -( [Kii] u_i + sum_j [Kij] u_j )
//==============================================================================
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <thread>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace tmass {

//------------------------------------------------------------------------------
// minimal linear algebra
//------------------------------------------------------------------------------
struct Vec3
{
    double x = 0.0, y = 0.0, z = 0.0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    double length() const { return std::sqrt(x * x + y * y + z * z); }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

struct Mat3
{
    double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    static Mat3 identity()
    {
        Mat3 r;
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1.0;
        return r;
    }
    // a b^T
    static Mat3 outer(const Vec3& a, const Vec3& b)
    {
        Mat3 r;
        const double av[3] = {a.x, a.y, a.z}, bv[3] = {b.x, b.y, b.z};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) r.m[i][j] = av[i] * bv[j];
        return r;
    }
    Mat3 operator+(const Mat3& o) const { Mat3 r; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) r.m[i][j] = m[i][j] + o.m[i][j]; return r; }
    Mat3 operator*(double s) const { Mat3 r; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) r.m[i][j] = m[i][j] * s; return r; }
    Mat3& operator+=(const Mat3& o) { for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i][j] += o.m[i][j]; return *this; }
    Mat3& operator-=(const Mat3& o) { for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) m[i][j] -= o.m[i][j]; return *this; }
    Vec3 operator*(const Vec3& v) const
    {
        return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
    }
    // transpose(M) * v without building the transpose
    Vec3 tmul(const Vec3& v) const
    {
        return {m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z,
                m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z,
                m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z};
    }
    Mat3 transposed() const { Mat3 r; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) r.m[i][j] = m[j][i]; return r; }
    double det() const
    {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }
    Mat3 inverse() const
    {
        const double d = det();
        if (std::fabs(d) < 1e-300) throw std::runtime_error("degenerate (zero-volume) tetrahedron");
        Mat3 r;
        r.m[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d;
        r.m[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / d;
        r.m[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d;
        r.m[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / d;
        r.m[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d;
        r.m[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / d;
        r.m[2][0] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d;
        r.m[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / d;
        r.m[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;
        return r;
    }
    double maxAbs() const
    {
        double a = 0.0;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) a = (std::max)(a, std::fabs(m[i][j]));
        return a;
    }
};

//------------------------------------------------------------------------------
// material helpers: Young's modulus E [Pa] and Poisson ratio nu -> Lame (lambda, mu)
//------------------------------------------------------------------------------
// small persistent worker pool: runs fn(begin, end) on contiguous chunks of
// [0, count). The calling thread does chunk 0 itself. Workers spin briefly and
// then sleep when idle, so the pool costs almost nothing between frames.
//------------------------------------------------------------------------------
class WorkerPool
{
public:
    explicit WorkerPool(int threads = 1) { setThreads(threads); }
    ~WorkerPool() { stop(); }
    WorkerPool(const WorkerPool&) = delete;
    WorkerPool& operator=(const WorkerPool&) = delete;

    void setThreads(int n)
    {
        stop();
        m_n = (std::max)(1, n);
        m_quit.store(false);
        m_gen.store(0);
        m_done.store(0);
        for (int k = 1; k < m_n; ++k) m_workers.emplace_back([this, k] { workerLoop(k); });
    }
    int threads() const { return m_n; }

    void parallelFor(size_t count, const std::function<void(size_t, size_t)>& fn)
    {
        if (m_n <= 1 || count < (size_t)(2 * m_n)) { fn(0, count); return; }
        m_fn = &fn;
        m_count = count;
        m_done.store(0, std::memory_order_relaxed);
        m_gen.fetch_add(1, std::memory_order_release);        // wake the workers
        runChunk(0, fn);
        while (m_done.load(std::memory_order_acquire) < m_n - 1) std::this_thread::yield();
    }

private:
    void runChunk(int k, const std::function<void(size_t, size_t)>& fn)
    {
        const size_t b = m_count * (size_t)k / (size_t)m_n;
        const size_t e = m_count * (size_t)(k + 1) / (size_t)m_n;
        if (e > b) fn(b, e);
    }
    void workerLoop(int k)
    {
        unsigned seen = 0;
        auto idleSince = std::chrono::steady_clock::now();
        unsigned spins = 0;
        while (!m_quit.load(std::memory_order_acquire))
        {
            const unsigned g = m_gen.load(std::memory_order_acquire);
            if (g != seen)
            {
                seen = g;
                runChunk(k, *m_fn);
                m_done.fetch_add(1, std::memory_order_release);
                idleSince = std::chrono::steady_clock::now();
                spins = 0;
                continue;
            }
            if ((++spins & 255u) == 0u &&
                std::chrono::steady_clock::now() - idleSince > std::chrono::milliseconds(2))
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            else
                std::this_thread::yield();
        }
    }
    void stop()
    {
        m_quit.store(true, std::memory_order_release);
        for (std::thread& t : m_workers) if (t.joinable()) t.join();
        m_workers.clear();
    }

    int m_n = 1;
    std::vector<std::thread> m_workers;
    const std::function<void(size_t, size_t)>* m_fn = nullptr;
    size_t m_count = 0;
    std::atomic<unsigned> m_gen{0};
    std::atomic<int> m_done{0};
    std::atomic<bool> m_quit{false};
};

//------------------------------------------------------------------------------
inline double lameLambda(double E, double nu) { return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); }
inline double lameMu(double E, double nu) { return E / (2.0 * (1.0 + nu)); }

//------------------------------------------------------------------------------
// tensor-mass model
//------------------------------------------------------------------------------
class TensorMassModel
{
public:
    using Tet = std::array<int, 4>;

    TensorMassModel() = default;

    // Build the model: validates the mesh, precomputes every [Kii] and [Kij].
    void init(const std::vector<Vec3>& restPositions, const std::vector<Tet>& tetrahedra,
              double lambda, double mu)
    {
        m_rest = restPositions;
        m_tets = tetrahedra;
        m_lambda = lambda;
        m_mu = mu;
        m_active.assign(m_tets.size(), 1);
        m_Kvert.assign(m_rest.size(), Mat3());
        m_Kedge.clear();
        m_edgeIndex.clear();
        m_edges.clear();

        for (size_t t = 0; t < m_tets.size(); t++)
        {
            for (int k = 0; k < 4; k++)
                if (m_tets[t][k] < 0 || m_tets[t][k] >= (int)m_rest.size())
                    throw std::runtime_error("tetrahedron references a vertex that does not exist");
            // register the six edges of this tetrahedron
            for (int k = 0; k < 4; k++)
                for (int l = k + 1; l < 4; l++) edgeSlot(m_tets[t][k], m_tets[t][l]);
            addTetContribution((int)t, +1.0);
        }
        buildAdjacency();
    }

    // Restoring elastic force on every vertex for current positions P (size = #vertices).
    // Cost is linear in the number of edges (thesis chapter 5).
    void computeForces(const std::vector<Vec3>& P, std::vector<Vec3>& F) const
    {
        const size_t n = m_rest.size();
        F.assign(n, Vec3());
        std::vector<Vec3> u(n);
        for (size_t i = 0; i < n; i++) u[i] = P[i] - m_rest[i];

        for (size_t i = 0; i < n; i++) F[i] -= m_Kvert[i] * u[i];
        for (size_t e = 0; e < m_edges.size(); e++)
        {
            const int i = m_edges[e].first, j = m_edges[e].second;
            const Mat3& Kij = m_Kedge[e];
            F[i] -= Kij * u[j];
            F[j] -= Kij.transposed() * u[i];   // [Kji] = [Kij]^T
        }
    }

    // Restoring force on ONE vertex, gathered from its own tensors only (no writes to
    // other vertices), so every vertex can be computed on a different thread.
    Vec3 forceAt(int i, const std::vector<Vec3>& P) const
    {
        return -rowKu(i, [&](int j) { return P[(size_t)j] - m_rest[(size_t)j]; });
    }
    // (K u)_i for a displacement field u (used by the time-step estimate)
    Vec3 stiffnessRow(int i, const std::vector<Vec3>& u) const
    {
        return rowKu(i, [&](int j) { return u[(size_t)j]; });
    }

    // Independent assembly tetrahedron-by-tetrahedron, used only to validate the
    // edge/vertex tensors above.
    void computeForcesByElement(const std::vector<Vec3>& P, std::vector<Vec3>& F) const
    {
        const size_t n = m_rest.size();
        F.assign(n, Vec3());
        for (size_t t = 0; t < m_tets.size(); t++)
        {
            if (!m_active[t]) continue;
            Mat3 K[4][4];
            elementBlocks((int)t, K);
            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++)
                {
                    const int a = m_tets[t][k], b = m_tets[t][l];
                    F[a] -= K[k][l] * (P[b] - m_rest[b]);
                }
        }
    }

    // Elastic energy W = 1/2 u^T K u = -1/2 sum_i u_i . F_i
    double elasticEnergy(const std::vector<Vec3>& P) const
    {
        std::vector<Vec3> F;
        computeForces(P, F);
        double w = 0.0;
        for (size_t i = 0; i < F.size(); i++) w -= 0.5 * dot(P[i] - m_rest[i], F[i]);
        return w;
    }

    // Cutting: subtract this tetrahedron's contribution from the vertex and edge
    // tensors. Local and cheap (4 vertex tensors + 6 edge tensors), and the
    // material then behaves exactly as if the tetrahedron had been absent from
    // the start (thesis section 4.4.3). Returns false if it was already removed.
    bool removeTetrahedron(int t)
    {
        if (t < 0 || t >= (int)m_tets.size() || !m_active[t]) return false;
        addTetContribution(t, -1.0);
        m_active[t] = 0;
        return true;
    }

    double tetVolume(int t) const
    {
        return std::fabs(restMatrix(t).det()) / 6.0;
    }

    size_t numVertices() const { return m_rest.size(); }
    size_t numTets() const { return m_tets.size(); }
    size_t numActiveTets() const { return (size_t)std::count(m_active.begin(), m_active.end(), (char)1); }
    size_t numEdges() const { return m_edges.size(); }
    const std::vector<Vec3>& restPositions() const { return m_rest; }
    const std::vector<Tet>& tets() const { return m_tets; }
    bool isActive(int t) const { return m_active[t] != 0; }

private:
    // matrix whose columns are the three edge vectors leaving vertex 0 of the tet
    Mat3 restMatrix(int t) const
    {
        const Vec3 e1 = m_rest[m_tets[t][1]] - m_rest[m_tets[t][0]];
        const Vec3 e2 = m_rest[m_tets[t][2]] - m_rest[m_tets[t][0]];
        const Vec3 e3 = m_rest[m_tets[t][3]] - m_rest[m_tets[t][0]];
        Mat3 D;
        D.m[0][0] = e1.x; D.m[0][1] = e2.x; D.m[0][2] = e3.x;
        D.m[1][0] = e1.y; D.m[1][1] = e2.y; D.m[1][2] = e3.y;
        D.m[2][0] = e1.z; D.m[2][1] = e2.z; D.m[2][2] = e3.z;
        return D;
    }

    // the 4x4 grid of 3x3 stiffness blocks of one tetrahedron
    void elementBlocks(int t, Mat3 K[4][4]) const
    {
        const Mat3 D = restMatrix(t);
        const double V = std::fabs(D.det()) / 6.0;
        const Mat3 Dinv = D.inverse();

        // gradients of the barycentric coordinates: rows of D^-1 for vertices 1..3
        Vec3 g[4];
        for (int k = 1; k < 4; k++) g[k] = Vec3(Dinv.m[k - 1][0], Dinv.m[k - 1][1], Dinv.m[k - 1][2]);
        g[0] = -(g[1] + g[2] + g[3]);

        for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
                K[k][l] = (Mat3::outer(g[k], g[l]) * m_lambda
                         + Mat3::outer(g[l], g[k]) * m_mu
                         + Mat3::identity() * (m_mu * dot(g[k], g[l]))) * V;
    }

    // add (sign=+1) or subtract (sign=-1) one tetrahedron from the stored tensors
    void addTetContribution(int t, double sign)
    {
        Mat3 K[4][4];
        elementBlocks(t, K);
        for (int k = 0; k < 4; k++)
        {
            const int a = m_tets[t][k];
            m_Kvert[a] += K[k][k] * sign;
            for (int l = k + 1; l < 4; l++)
            {
                const int b = m_tets[t][l];
                const size_t e = edgeSlot(a, b);
                // the stored tensor is [K_ij] with i = m_edges[e].first (the smaller index)
                if (a == m_edges[e].first) m_Kedge[e] += K[k][l] * sign;
                else                       m_Kedge[e] += K[l][k] * sign;   // = K[k][l]^T
            }
        }
    }

    // index of the edge (a,b); created on first use
    size_t edgeSlot(int a, int b)
    {
        const int i = (std::min)(a, b), j = (std::max)(a, b);
        const uint64_t key = ((uint64_t)(uint32_t)i << 32) | (uint32_t)j;
        auto it = m_edgeIndex.find(key);
        if (it != m_edgeIndex.end()) return it->second;
        const size_t idx = m_edges.size();
        m_edgeIndex[key] = idx;
        m_edges.push_back({i, j});
        m_Kedge.push_back(Mat3());
        return idx;
    }

    // CSR adjacency: for each vertex, the edges that touch it (the edge list never
    // changes after init; cutting only changes the tensor values)
    void buildAdjacency()
    {
        const size_t n = m_rest.size();
        m_adjStart.assign(n + 1, 0);
        for (const auto& e : m_edges) { m_adjStart[(size_t)e.first + 1]++; m_adjStart[(size_t)e.second + 1]++; }
        for (size_t i = 0; i < n; i++) m_adjStart[i + 1] += m_adjStart[i];
        m_adjEdge.assign(m_adjStart[n], 0);
        m_adjOther.assign(m_adjStart[n], 0);
        std::vector<int> fill(m_adjStart.begin(), m_adjStart.end() - 1);
        for (size_t e = 0; e < m_edges.size(); e++)
        {
            const int a = m_edges[e].first, b = m_edges[e].second;
            m_adjEdge[(size_t)fill[(size_t)a]] = (int)e; m_adjOther[(size_t)fill[(size_t)a]++] = b;
            m_adjEdge[(size_t)fill[(size_t)b]] = (int)e; m_adjOther[(size_t)fill[(size_t)b]++] = a;
        }
    }
    template <class U>
    Vec3 rowKu(int i, const U& u) const
    {
        Vec3 f = m_Kvert[(size_t)i] * u(i);
        for (int k = m_adjStart[(size_t)i]; k < m_adjStart[(size_t)i + 1]; k++)
        {
            const int e = m_adjEdge[(size_t)k], j = m_adjOther[(size_t)k];
            if (m_edges[(size_t)e].first == i) f += m_Kedge[(size_t)e] * u(j);   // [Kij]
            else                               f += m_Kedge[(size_t)e].tmul(u(j)); // [Kji] = [Kij]^T
        }
        return f;
    }

    std::vector<Vec3> m_rest;
    std::vector<Tet> m_tets;
    std::vector<char> m_active;
    double m_lambda = 0.0, m_mu = 0.0;

    std::vector<Mat3> m_Kvert;                      // [Kii], one per vertex
    std::vector<Mat3> m_Kedge;                      // [Kij], one per edge (i < j)
    std::vector<std::pair<int, int>> m_edges;       // edge -> (i, j), i < j
    std::unordered_map<uint64_t, size_t> m_edgeIndex;
    std::vector<int> m_adjStart, m_adjEdge, m_adjOther;   // CSR vertex -> incident edges
};

//------------------------------------------------------------------------------
// test mesh: nx*ny*nz cubes of side h, each cut into 6 tetrahedra (Kuhn split,
// which is conforming across neighbouring cubes)
//------------------------------------------------------------------------------
inline void makeBlockMesh(int nx, int ny, int nz, double h,
                          std::vector<Vec3>& verts, std::vector<TensorMassModel::Tet>& tets)
{
    verts.clear();
    tets.clear();
    auto id = [&](int i, int j, int k) { return (k * (ny + 1) + j) * (nx + 1) + i; };
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++) verts.push_back(Vec3(i * h, j * h, k * h));

    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                for (int p = 0; p < 6; p++)
                {
                    int c[3] = {i, j, k};
                    TensorMassModel::Tet T;
                    T[0] = id(c[0], c[1], c[2]);
                    for (int s = 0; s < 3; s++)
                    {
                        c[perms[p][s]] += 1;
                        T[s + 1] = id(c[0], c[1], c[2]);
                    }
                    tets.push_back(T);
                }
}

} // namespace tmass
