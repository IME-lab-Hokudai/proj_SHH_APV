#pragma once
// SH deringing (P.-P. Sloan, "Deringing Spherical Harmonics", SIGGRAPH Asia 2017 Technical Briefs).
//
// With nine coefficients, the irradiance of a probe close to a small, bright emitter rings: it goes below zero on
// normals facing away from the emitter (about -4% of the peak for a point source). Sloan windows the stored radiance
// with the least amount of windowing that makes the irradiance non-negative:
//   - one-parameter family of per-band scales (Table 1 of the paper); band 0 is never scaled, so the mean irradiance
//     is preserved; the window is applied to the radiance coefficients, before the cosine convolution (A_l);
//   - the least window is found per colour channel by binary search with linear interpolation between table columns,
//     and the strongest of the three is applied to all three channels (no colour shift);
//   - the strongest table window makes the irradiance of a delta (point) light non-negative, so it does so for any
//     non-negative radiance (a sum of deltas), i.e. for every probe.
// Table 1 lists rows l1, l2, l3 for w = inf, 16.7, 11.3, 10, 9, 7, 5.6. Rows l2 and l3 are exactly sinc(pi l / w) at
// l = 2, 3; used as the scales of bands 1 and 2 they satisfy the paper's guarantee (delta irradiance at w = 5.6:
// min / peak = +1.45%), while rows l1, l2 do not (-0.79%), so bands 1 and 2 take rows l2 and l3 here.
// Non-negativity is tested with the exact minimum of the irradiance over the sphere (a quadratic in the normal:
// 3x3 eigenvalues and a bisection) instead of the paper's conservative bound, so the window is never stronger than
// needed. Beyond the table (t in (6, 7]) bands 1 and 2 fade linearly to zero, which only the shader uses for
// interpolated SH that is not a projection of non-negative radiance (band 0 alone is non-negative when it is).
// Same basis and A_l as evaluateIrradiance in XAtlasSceneShader.slang (Condon-Shortley phase, z-up basis on Falcor's
// (z, x, y)); the shader side is SHNonNegative.slang.
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace SHNonNegative
{
using Vec9 = std::array<double, 9>;

inline Vec9 basis(double nx, double ny, double nz)
{
    const double x = nz, y = nx, z = ny;
    return { 0.2820947917738781,
             -0.4886025119029200 * y,
             0.4886025119029200 * z,
             -0.4886025119029200 * x,
             1.0925484305920792 * x * y,
             -1.0925484305920792 * y * z,
             0.9461746957575601 * z * z - 0.3153915652525200,
             -1.0925484305920792 * x * z,
             0.546274 * (x * x - y * y) };
}

inline double dot9(const Vec9& a, const Vec9& b)
{
    double s = 0.0;
    for (int j = 0; j < 9; ++j) s += a[j] * b[j];
    return s;
}

// E(n) = k0 + b.n + n^T Q n in Falcor world coordinates, for u = A_l c.
inline void quadratic(const Vec9& u, double& k0, double b[3], double Q[3][3])
{
    k0 = 0.2820947917738781 * u[0] - 0.3153915652525200 * u[6];
    b[0] = -0.4886025119029200 * u[1];
    b[1] = 0.4886025119029200 * u[2];
    b[2] = -0.4886025119029200 * u[3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Q[i][j] = 0.0;
    auto sym = [&](int i, int j, double v) { Q[i][j] += 0.5 * v; Q[j][i] += 0.5 * v; };
    sym(2, 0, 1.0925484305920792 * u[4]);  // x y   = n.z n.x
    sym(0, 1, -1.0925484305920792 * u[5]); // y z   = n.x n.y
    Q[1][1] += 0.9461746957575601 * u[6];  // z^2   = n.y^2
    sym(2, 1, -1.0925484305920792 * u[7]); // x z   = n.z n.y
    Q[2][2] += 0.546274 * u[8];            // x^2   = n.z^2
    Q[0][0] -= 0.546274 * u[8];            // -y^2  = -n.x^2
}

// Symmetric 3x3 eigen-decomposition (cyclic Jacobi): Q = V diag(w) V^T, columns of V are the eigenvectors.
inline void eigenSym3(const double Qin[3][3], double w[3], double V[3][3])
{
    double a[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) { a[i][j] = Qin[i][j]; V[i][j] = i == j ? 1.0 : 0.0; }
    for (int sweep = 0; sweep < 50; ++sweep)
    {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off < 1e-30 * (a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2]) + 1e-300) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q)
            {
                if (a[p][q] == 0.0) continue;
                const double th = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
                const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) // a = J^T a J
                {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
                }
            }
    }
    for (int i = 0; i < 3; ++i) w[i] = a[i][i];
}

// Global minimum of E(n) = k0 + b.n + n^T Q n on the unit sphere (trust-region subproblem): n = -(Q + mu I)^-1 b / 2
// with mu >= -lambda_min(Q) chosen so that |n| = 1; in the hard case (b has no component along the lowest eigenvector)
// the remaining length goes along that eigenvector. Writes the minimiser to n, returns E there.
inline double globalMinimum(const Vec9& u, double n[3])
{
    double k0, b[3], Q[3][3], w[3], V[3][3];
    quadratic(u, k0, b, Q);
    eigenSym3(Q, w, V);
    int lo = 0;
    for (int i = 1; i < 3; ++i) if (w[i] < w[lo]) lo = i;
    double beta[3]; // b in the eigenbasis
    for (int i = 0; i < 3; ++i) beta[i] = b[0] * V[0][i] + b[1] * V[1][i] + b[2] * V[2][i];
    const double scale = std::fabs(w[0]) + std::fabs(w[1]) + std::fabs(w[2]) + std::fabs(b[0]) + std::fabs(b[1]) + std::fabs(b[2]) + 1e-300;
    auto len2 = [&](double mu)
    {
        double s = 0.0;
        for (int i = 0; i < 3; ++i) { const double d = 2.0 * (w[i] + mu); const double x = beta[i] / d; s += x * x; }
        return s;
    };
    double y[3];
    const double muLo = -w[lo];
    const double tiny = 1e-12 * scale;
    // |n(mu)| decreases in mu on (-lambda_min, inf) and |n| <= |b| / (2 (lambda_min + mu)): bisect for |n| = 1.
    double lower = muLo + tiny, upper = muLo + 0.5 * std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) + tiny;
    for (int it = 0; it < 200; ++it)
    {
        const double m = 0.5 * (lower + upper);
        (len2(m) > 1.0 ? lower : upper) = m;
        if (upper - lower < 1e-15 * (1.0 + std::fabs(upper))) break;
    }
    for (int i = 0; i < 3; ++i) y[i] = -beta[i] / (2.0 * (w[i] + upper));
    // (Nearly) hard case: b has (almost) no component along the lowest eigenvector, so |y| < 1 even at mu -> -lambda_min.
    // Keep the other components and give the lowest one the rest of the unit length, with the sign that lowers E.
    double r2 = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
    if (r2 < 1.0)
    {
        const double rest = r2 - y[lo] * y[lo];
        y[lo] = (beta[lo] > 0.0 ? -1.0 : 1.0) * std::sqrt(std::max(0.0, 1.0 - rest));
    }
    double len = 0.0;
    for (int i = 0; i < 3; ++i) { n[i] = V[i][0] * y[0] + V[i][1] * y[1] + V[i][2] * y[2]; len += n[i] * n[i]; }
    len = std::sqrt(len);
    for (int i = 0; i < 3; ++i) n[i] /= len;
    double e = k0 + b[0] * n[0] + b[1] * n[1] + b[2] * n[2];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) e += n[i] * Q[i][j] * n[j];
    return e;
}

// Per-band window scales (band 1, band 2) at parameter t: 0 = no window, 6 = the strongest column of Sloan's Table 1
// (w = 5.6), linear interpolation in between; (6, 7] fades both to zero.
inline void windowScales(double t, double& s1, double& s2)
{
    static const double kBand1[7] = { 1.0, 0.9766, 0.9493, 0.9355, 0.9207, 0.8710, 0.8030 }; // Table 1 row l2
    static const double kBand2[7] = { 1.0, 0.9478, 0.8880, 0.8584, 0.8270, 0.7241, 0.5904 }; // Table 1 row l3
    t = std::clamp(t, 0.0, 7.0);
    if (t >= 6.0)
    {
        const double f = 7.0 - t;
        s1 = kBand1[6] * f;
        s2 = kBand2[6] * f;
        return;
    }
    const int i = std::min(int(t), 5);
    const double f = t - i;
    s1 = kBand1[i] + f * (kBand1[i + 1] - kBand1[i]);
    s2 = kBand2[i] + f * (kBand2[i + 1] - kBand2[i]);
}

// Minimum over the sphere of one channel's irradiance with the window (s1, s2).
inline double windowedIrradianceMin(const Vec9& c, double s1, double s2)
{
    static const double A[9] = { 3.141593, 2.094395, 2.094395, 2.094395, 0.785398, 0.785398, 0.785398, 0.785398, 0.785398 };
    Vec9 u;
    for (int j = 0; j < 9; ++j) u[j] = A[j] * c[j] * (j == 0 ? 1.0 : j < 4 ? s1 : s2);
    double n[3];
    return globalMinimum(u, n);
}

// Least window parameter t in [0, tMax] for which this channel's irradiance is >= 0 (0 if it already is).
inline double leastWindow(const Vec9& c, double tMax = 6.0)
{
    const double tol = 1e-6 * std::fabs(0.2820947917738781 * 3.141593 * c[0]);
    double s1, s2;
    if (windowedIrradianceMin(c, 1.0, 1.0) >= -tol) return 0.0;
    windowScales(tMax, s1, s2);
    if (windowedIrradianceMin(c, s1, s2) < -tol) return tMax; // cannot be reached within tMax: strongest allowed
    double lo = 0.0, hi = tMax;
    for (int it = 0; it < 30; ++it)
    {
        const double m = 0.5 * (lo + hi);
        windowScales(m, s1, s2);
        (windowedIrradianceMin(c, s1, s2) >= -tol ? hi : lo) = m;
    }
    return hi;
}

// Sloan deringing of one probe: the least window per channel, the strongest of the three for all channels.
// coeffs: nine RGB coefficients (any type with .x .y .z). Returns false if no window was needed; otherwise sets the
// band-1 and band-2 scales and applies them to coeffs. The caller scales the probe's gradients by the same factors
// (the windowed field near the probe).
template<typename V3>
inline bool deringProbe(V3* coeffs, double& s1, double& s2)
{
    double t = 0.0;
    for (int ch = 0; ch < 3; ++ch)
    {
        Vec9 c;
        for (int j = 0; j < 9; ++j) c[j] = ch == 0 ? coeffs[j].x : ch == 1 ? coeffs[j].y : coeffs[j].z;
        t = std::max(t, leastWindow(c));
    }
    s1 = s2 = 1.0;
    if (t <= 0.0) return false;
    windowScales(t, s1, s2);
    for (int j = 1; j < 9; ++j)
    {
        const float s = float(j < 4 ? s1 : s2);
        coeffs[j].x *= s;
        coeffs[j].y *= s;
        coeffs[j].z *= s;
    }
    return true;
}
} // namespace SHNonNegative
