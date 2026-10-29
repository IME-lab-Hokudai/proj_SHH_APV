// FDVerificationCore.h
// Finite-difference verification of the analytic SH gradient / Hessian (E1a, E1c).
// Header-only, no Falcor dependency, so the maths can be unit-tested on its own.
//
// All positions and directions are in the pipeline's polar frame (Z-up):
//   polar = (world.z, world.x, world.y)
//
// Two ways to evaluate the SH coefficients away from the traced position x0:
//   * frozen  (E1a): keep every hit point s_i, its normal n_i and its radiance L_i,
//                    and move only the probe (Krivanek et al. 2005, Sec. 4.3.1).
//                    Its derivatives are exactly what the analytic formulas compute.
//   * retraced (E1c): trace the probe again at the shifted position and use the
//                    pipeline's own estimator. Its derivatives include visibility changes.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace fdv
{

// ---------------------------------------------------------------------------
// Minimal double-precision 3-vector
// ---------------------------------------------------------------------------
template<typename T>
struct V3
{
    T x{}, y{}, z{};
    V3() = default;
    V3(T a, T b, T c) : x(a), y(b), z(c) {}
    T& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    const T& operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};
template<typename T> V3<T> operator+(const V3<T>& a, const V3<T>& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
template<typename T> V3<T> operator-(const V3<T>& a, const V3<T>& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
template<typename T, typename S> V3<T> operator*(const V3<T>& a, const S& s) { return { a.x * s, a.y * s, a.z * s }; }
template<typename T> T dot(const V3<T>& a, const V3<T>& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

using Vec3d = V3<double>;

// One record of the probe trace, converted to double.
// Direction samples: L = radiance leaving the first hit toward the probe, s/n = hit point/normal,
//                    hit = false for a ray miss (no environment light: L = 0).
// Light samples:     L = emitted radiance Le, s/n = point/normal on the emitter,
//                    dA = area weight 1/(M pdfArea), hit = false if invalid or occluded.
struct Sample
{
    Vec3d L, s, n;
    double dA = 0.0;
    bool hit = false;
};

// 9 SH coefficients, each RGB.
using CoeffsRGB = std::array<Vec3d, 9>;

constexpr double kPi = 3.14159265358979323846;
const Vec3d kLuma{ 0.2126, 0.7152, 0.0722 };

// ---------------------------------------------------------------------------
// Real SH up to l = 2, same convention and order (index l(l+1)+m) as
// SHGradientAndHessianL2 / SHGradientL2 in envMap_SH.cpp. Templated so the
// unit test can push automatic-differentiation types through it.
// ---------------------------------------------------------------------------
template<typename T>
void shBasisL2(const V3<T>& w, std::array<T, 9>& y)
{
    const T x = w.x, yy = w.y, z = w.z;
    y[0] = T(0.2820947917738781);
    y[1] = T(-0.4886025119029199) * yy;
    y[2] = T(0.4886025119029199) * z;
    y[3] = T(-0.4886025119029199) * x;
    y[4] = T(1.0925484305920792) * x * yy;
    y[5] = T(-1.0925484305920792) * yy * z;
    y[6] = T(0.9461746957575601) * z * z - T(0.3153915652525200);
    y[7] = T(-1.0925484305920792) * x * z;
    y[8] = T(0.5462742152960396) * (x * x - yy * yy);
}

inline double value(double v) { return v; }
using std::sqrt;

struct FrozenStats
{
    uint32_t grazing = 0;     // direction hits with |cos xi| < 1e-4 at x0 (ill-conditioned patches)
    uint32_t backFacing = 0;  // patches that face away from the evaluation point (cos xi' <= 0)
};

// ---------------------------------------------------------------------------
// Frozen-patch coefficients at xe from samples traced at x0.
//   direction hit i:  Omega_i(xe) = (4 pi / N) * (r0^2 cos xi'(xe)) / (r'^2 cos xi0)
//                     i.e. the patch area dA_i = (4 pi / N) r0^2 / cos xi0 is held fixed
//   direction miss:   point at infinity, same direction and weight 4 pi / N everywhere
//   light sample:     Omega(xe) = dA cos xi'(xe) / r'^2, area dA held fixed
// At xe == x0 this reproduces calculateSHCoeffs + accumulateDirectLightSamples.
// Templated on the scalar of xe so automatic differentiation can be used in tests.
// ---------------------------------------------------------------------------
template<typename T>
void frozenCoeffsT(const Sample* dirSamples, uint32_t dirCount, const Vec3d* dirs,
                   const Sample* lightSamples, uint32_t lightCount,
                   const Vec3d& x0, const V3<T>& xe,
                   std::array<V3<T>, 9>& out, FrozenStats* stats = nullptr)
{
    for (auto& c : out) c = V3<T>(T(0), T(0), T(0));
    std::array<T, 9> y;
    const double omegaUniform = 4.0 * kPi / double(dirCount);

    for (uint32_t i = 0; i < dirCount; ++i)
    {
        const Sample& sd = dirSamples[i];
        if (!sd.hit)
        {
            // Without an environment light L is 0; kept general for safety.
            const V3<T> w(T(dirs[i].x), T(dirs[i].y), T(dirs[i].z));
            shBasisL2(w, y);
            for (int b = 0; b < 9; ++b)
                out[b] = out[b] + V3<T>(T(sd.L.x), T(sd.L.y), T(sd.L.z)) * (y[b] * T(omegaUniform));
            continue;
        }
        const Vec3d q0 = sd.s - x0;
        const double r0sq = dot(q0, q0);
        const double r0 = std::sqrt(r0sq);
        const double cos0 = -dot(sd.n, q0) / r0;
        if (stats && std::fabs(cos0) < 1e-4) ++stats->grazing;
        if (cos0 == 0.0) continue;
        const double patchArea = omegaUniform * r0sq / cos0; // dA_i, held fixed

        const V3<T> q(T(sd.s.x) - xe.x, T(sd.s.y) - xe.y, T(sd.s.z) - xe.z);
        const T rsq = q.x * q.x + q.y * q.y + q.z * q.z;
        const T r = sqrt(rsq);
        const T cosE = -(T(sd.n.x) * q.x + T(sd.n.y) * q.y + T(sd.n.z) * q.z) / r;
        if (stats && value(cosE) <= 0.0) ++stats->backFacing;
        const T omega = T(patchArea) * cosE / rsq;
        const V3<T> w(q.x / r, q.y / r, q.z / r);
        shBasisL2(w, y);
        for (int b = 0; b < 9; ++b)
            out[b] = out[b] + V3<T>(T(sd.L.x), T(sd.L.y), T(sd.L.z)) * (y[b] * omega);
    } 

    for (uint32_t i = 0; i < lightCount; ++i)
    {
        const Sample& sd = lightSamples[i];
        if (!sd.hit || sd.dA <= 0.0) continue;
        // The pipeline drops light points behind the emitter at x0; keep the same set.
        const Vec3d q0 = sd.s - x0;
        const double r0 = std::sqrt(dot(q0, q0));
        if (r0 < 1e-6 || -dot(sd.n, q0) / r0 <= 0.0) continue;

        const V3<T> q(T(sd.s.x) - xe.x, T(sd.s.y) - xe.y, T(sd.s.z) - xe.z);
        const T rsq = q.x * q.x + q.y * q.y + q.z * q.z;
        const T r = sqrt(rsq);
        const T cosE = -(T(sd.n.x) * q.x + T(sd.n.y) * q.y + T(sd.n.z) * q.z) / r;
        if (value(cosE) <= 0.0) { if (stats) ++stats->backFacing; continue; }
        const T omega = T(sd.dA) * cosE / rsq;
        const V3<T> w(q.x / r, q.y / r, q.z / r);
        shBasisL2(w, y);
        for (int b = 0; b < 9; ++b)
            out[b] = out[b] + V3<T>(T(sd.L.x), T(sd.L.y), T(sd.L.z)) * (y[b] * omega);
    }
}

inline void frozenCoeffs(const Sample* dirSamples, uint32_t dirCount, const Vec3d* dirs,
                         const Sample* lightSamples, uint32_t lightCount,
                         const Vec3d& x0, const Vec3d& xe, CoeffsRGB& out, FrozenStats* stats = nullptr)
{
    frozenCoeffsT<double>(dirSamples, dirCount, dirs, lightSamples, lightCount, x0, xe, out, stats);
}

// Coefficients exactly as the pipeline estimates them at the traced position.
inline void tracedCoeffs(const Sample* dirSamples, uint32_t dirCount, const Vec3d* dirs,
                         const Sample* lightSamples, uint32_t lightCount,
                         const Vec3d& x, CoeffsRGB& out)
{
    frozenCoeffs(dirSamples, dirCount, dirs, lightSamples, lightCount, x, x, out, nullptr);
}

// ---------------------------------------------------------------------------
// 19-point central stencil.
//   0        centre
//   1..6     +h e_j, -h e_j          (j = 0,1,2 -> index 1+2j, 2+2j)
//   7..18    (+-h e_j, +-h e_k), j<k, order ++, +-, -+, --
// ---------------------------------------------------------------------------
constexpr int kStencilSize = 19;

inline std::array<Vec3d, kStencilSize> stencilOffsets(double h)
{
    std::array<Vec3d, kStencilSize> o{};
    int idx = 1;
    for (int j = 0; j < 3; ++j)
    {
        Vec3d p, m;
        p[j] = h; m[j] = -h;
        o[idx++] = p;
        o[idx++] = m;
    }
    const int pairs[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
    const double sg[4][2] = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
    for (auto& pr : pairs)
        for (auto& s : sg)
        {
            Vec3d d;
            d[pr[0]] = s[0] * h;
            d[pr[1]] = s[1] * h;
            o[idx++] = d;
        }
    return o;
}

// Central-difference gradient and Hessian of a scalar sampled on the stencil.
// Truncation error O(h^2) for both.
inline void fdDerivatives(const std::array<double, kStencilSize>& f, double h, double g[3], double H[3][3])
{
    const double h2 = h * h;
    for (int j = 0; j < 3; ++j)
    {
        const double fp = f[1 + 2 * j], fm = f[2 + 2 * j];
        g[j] = (fp - fm) / (2.0 * h);
        H[j][j] = (fp - 2.0 * f[0] + fm) / h2;
    }
    const int pairs[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
    for (int p = 0; p < 3; ++p)
    {
        const int base = 7 + 4 * p;
        const double v = (f[base] - f[base + 1] - f[base + 2] + f[base + 3]) / (4.0 * h2);
        H[pairs[p][0]][pairs[p][1]] = v;
        H[pairs[p][1]][pairs[p][0]] = v;
    }
}

// Derivatives of all coefficients from the 19 stencil evaluations.
//   gradRGB[b][j] = d f_b / d x_j  (RGB)
//   hessLum[b][j][k] = d^2 (luma . f_b) / d x_j d x_k   (what the pipeline stores)
struct FDDerivs
{
    std::array<std::array<Vec3d, 3>, 9> gradRGB{};
    std::array<std::array<std::array<double, 3>, 3>, 9> hessLum{};
};

inline FDDerivs fdFromStencil(const std::array<CoeffsRGB, kStencilSize>& c, double h)
{
    FDDerivs d;
    std::array<double, kStencilSize> f;
    double g[3], H[3][3];
    for (int b = 0; b < 9; ++b)
    {
        for (int ch = 0; ch < 3; ++ch)
        {
            for (int s = 0; s < kStencilSize; ++s) f[s] = c[s][b][ch];
            fdDerivatives(f, h, g, H);
            for (int j = 0; j < 3; ++j) d.gradRGB[b][j][ch] = g[j];
        }
        for (int s = 0; s < kStencilSize; ++s) f[s] = dot(c[s][b], kLuma);
        fdDerivatives(f, h, g, H);
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) d.hessLum[b][j][k] = H[j][k];
    }
    return d;
}

// Frozen-patch FD (E1a): one trace at x0, 19 reprojections.
inline FDDerivs frozenFD(const Sample* dirSamples, uint32_t dirCount, const Vec3d* dirs,
                         const Sample* lightSamples, uint32_t lightCount,
                         const Vec3d& x0, double h, FrozenStats* stats = nullptr)
{
    const auto offs = stencilOffsets(h);
    std::array<CoeffsRGB, kStencilSize> c;
    for (int s = 0; s < kStencilSize; ++s)
        frozenCoeffs(dirSamples, dirCount, dirs, lightSamples, lightCount, x0, x0 + offs[s], c[s], s == 0 ? stats : nullptr);
    return fdFromStencil(c, h);
}

// ---------------------------------------------------------------------------
// Ray flips between the centre trace (probe at x0) and a shifted trace (x0 + d) (E1c).
// A direction sample "flips" if its hit/miss status changes, or if its hit point
// jumps further than a continuous surface allows. On one continuous surface the
// hit point of a fixed direction moves by about |d| / cos(xi); a jump to another
// depth (occluder <-> background) moves it by a finite distance. Creases between
// walls (box corners) are continuous and are not counted.
// A light sample flips if its visibility changes.
// ---------------------------------------------------------------------------
struct FlipCount
{
    uint32_t dirFlips = 0, dirHits = 0;
    uint32_t lightFlips = 0, lightValid = 0;
};

inline FlipCount countFlips(const Sample* dirA, const Sample* dirB, uint32_t dirCount,
                            const Sample* lightA, const Sample* lightB, uint32_t lightCount,
                            const Vec3d& x0, double shiftLength, double jumpFactor = 10.0)
{
    FlipCount fc;
    for (uint32_t i = 0; i < dirCount; ++i)
    {
        const Sample& a = dirA[i];
        const Sample& b = dirB[i];
        if (a.hit || b.hit) ++fc.dirHits;
        if (a.hit != b.hit) { ++fc.dirFlips; continue; }
        if (!a.hit) continue;
        const Vec3d q0 = a.s - x0;
        const double r0 = std::sqrt(dot(q0, q0));
        const double cos0 = std::max(std::fabs(dot(a.n, q0)) / r0, 0.05);
        const Vec3d ds = b.s - a.s;
        if (std::sqrt(dot(ds, ds)) > jumpFactor * shiftLength / cos0) ++fc.dirFlips;
    }
    for (uint32_t i = 0; i < lightCount; ++i)
    {
        if (lightA[i].dA <= 0.0 && lightB[i].dA <= 0.0) continue; // never a valid light point
        ++fc.lightValid;
        if (lightA[i].hit != lightB[i].hit) ++fc.lightFlips;
    }
    return fc;
}

// ---------------------------------------------------------------------------
// Error between analytic and FD derivatives.
// Relative errors use Frobenius norms over all entries, per quantity and per band.
// ---------------------------------------------------------------------------
struct AnalyticDerivs
{
    std::array<std::array<Vec3d, 3>, 9> gradRGB{};                // [b][j] rgb
    std::array<std::array<std::array<double, 3>, 3>, 9> hessLum{}; // [b][j][k]
};

struct ErrorStats
{
    double gradNormAn = 0, gradErr = 0;   // Frobenius norms
    double hessNormAn = 0, hessErr = 0;
    double hessNormBand[3] = { 0, 0, 0 }, hessErrBand[3] = { 0, 0, 0 };
    double gradRel() const { return gradNormAn > 0 ? gradErr / gradNormAn : 0; }
    double hessRel() const { return hessNormAn > 0 ? hessErr / hessNormAn : 0; }
    double hessRelBand(int l) const { return hessNormBand[l] > 0 ? hessErrBand[l] / hessNormBand[l] : 0; }
};

inline int bandOf(int b) { return b == 0 ? 0 : (b < 4 ? 1 : 2); }

inline ErrorStats compare(const AnalyticDerivs& an, const FDDerivs& fd)
{
    ErrorStats e;
    for (int b = 0; b < 9; ++b)
    {
        for (int j = 0; j < 3; ++j)
            for (int ch = 0; ch < 3; ++ch)
            {
                const double a = an.gradRGB[b][j][ch], d = fd.gradRGB[b][j][ch] - a;
                e.gradNormAn += a * a;
                e.gradErr += d * d;
            }
        const int l = bandOf(b);
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
            {
                const double a = an.hessLum[b][j][k], d = fd.hessLum[b][j][k] - a;
                e.hessNormAn += a * a;
                e.hessErr += d * d;
                e.hessNormBand[l] += a * a;
                e.hessErrBand[l] += d * d;
            }
    }
    e.gradNormAn = std::sqrt(e.gradNormAn);
    e.gradErr = std::sqrt(e.gradErr);
    e.hessNormAn = std::sqrt(e.hessNormAn);
    e.hessErr = std::sqrt(e.hessErr);
    for (int l = 0; l < 3; ++l)
    {
        e.hessNormBand[l] = std::sqrt(e.hessNormBand[l]);
        e.hessErrBand[l] = std::sqrt(e.hessErrBand[l]);
    }
    return e;
}

} // namespace fdv
