// SH deringing at shading points (SH_NONNEG), the runtime half of SHNonNegative.h.
//
// The probes are windowed with Sloan's method when the grid is loaded ("Deringing Spherical Harmonics", SIGGRAPH Asia
// 2017; SHNonNegative.h), so every probe's irradiance is non-negative. Hermite interpolation between them can still
// overshoot below zero (next to a small, bright emitter the probe gradients are huge). At a shading point we therefore
// scale bands 1 and 2 of the interpolated SH by the largest factor a in [0, 1] that keeps its irradiance non-negative,
// the same factor for R, G and B. With E(n) = E0 + f(n) (E0: the band-0 part, f: bands 1-2), the scaled irradiance is
// E0 + a f(n), so a = E0 / (E0 - min E) exactly, from one exact minimum per channel: no search. A cheap lower bound
// (Gershgorin) skips the minimum where the irradiance clearly cannot go negative.

// E(n) = k0 + dot(b, n) + dot(n, mul(Q, n)) in Falcor world coordinates for channel ch, with band weights A (as used
// by evaluateIrradiance, including any SH window). Basis as evaluateIrradiance: (x, y, z) = (n.z, n.x, n.y).
void shIrradianceQuadratic(float3 sh[9], int ch, float3 A, out float k0, out float3 b, out float3x3 Q)
{
    float u[9];
    [unroll]
    for (int j = 0; j < 9; ++j) u[j] = sh[j][ch] * (j == 0 ? A.x : j < 4 ? A.y : A.z);
    k0 = 0.2820947917738781f * u[0] - 0.3153915652525200f * u[6];
    b = float3(-0.4886025119029200f * u[1], 0.4886025119029200f * u[2], -0.4886025119029200f * u[3]);
    const float xy = 0.5f * 1.0925484305920792f;
    Q = float3x3(0.f);
    Q[2][0] += xy * u[4]; Q[0][2] += xy * u[4];   // x y  = n.z n.x
    Q[0][1] -= xy * u[5]; Q[1][0] -= xy * u[5];   // y z  = n.x n.y (negative sign in the basis)
    Q[1][1] += 0.9461746957575601f * u[6];        // z^2  = n.y^2
    Q[2][1] -= xy * u[7]; Q[1][2] -= xy * u[7];   // x z  = n.z n.y (negative sign in the basis)
    Q[2][2] += 0.546274f * u[8];                  // x^2  = n.z^2
    Q[0][0] -= 0.546274f * u[8];                  // -y^2 = -n.x^2
}

// One Jacobi rotation in the (p, q) plane of the symmetric a; accumulates the eigenvectors in V.
void shJacobiRotate(inout float3x3 a, inout float3x3 V, int p, int q)
{
    float apq = a[p][q];
    if (abs(apq) < 1e-30f) return;
    float th = 0.5f * (a[q][q] - a[p][p]) / apq;
    float t = (th >= 0.f ? 1.f : -1.f) / (abs(th) + sqrt(th * th + 1.f));
    float c = 1.f / sqrt(t * t + 1.f), s = t * c;
    [unroll]
    for (int k = 0; k < 3; ++k)
    {
        float akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
    }
    [unroll]
    for (int k = 0; k < 3; ++k)
    {
        float apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
    }
    [unroll]
    for (int k = 0; k < 3; ++k)
    {
        float vkp = V[k][p], vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
    }
}

// Minimum over the unit sphere of E(n) = k0 + b.n + n^T Q n (trust-region subproblem).
float shQuadraticMinOnSphere(float k0, float3 b, float3x3 Q)
{
    float3x3 a = Q;
    float3x3 V = float3x3(1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f);
    for (int sweep = 0; sweep < 4; ++sweep)
    {
        shJacobiRotate(a, V, 0, 1);
        shJacobiRotate(a, V, 0, 2);
        shJacobiRotate(a, V, 1, 2);
    }
    float3 w = float3(a[0][0], a[1][1], a[2][2]);
    float3 beta = mul(b, V); // b in the eigenbasis (columns of V)
    float wLo = min(w.x, min(w.y, w.z));
    float scale = abs(w.x) + abs(w.y) + abs(w.z) + abs(b.x) + abs(b.y) + abs(b.z) + 1e-30f;
    // |y(mu)|^2 with y_i = -beta_i / (2 (w_i + mu)) decreases on (-wLo, inf); |y| <= |b| / (2 (wLo + mu)).
    float lo = -wLo + 1e-6f * scale, hi = -wLo + 0.5f * length(b) + 1e-6f * scale;
    float3 y;
    for (int it = 0; it < 28; ++it)
    {
        float mu = 0.5f * (lo + hi);
        y = -beta / (2.f * (w + mu));
        if (dot(y, y) > 1.f) lo = mu; else hi = mu;
    }
    y = -beta / (2.f * (w + hi));
    // (Nearly) hard case: b has (almost) no component along the lowest eigenvector, so |y| < 1 even at mu -> -wLo. Keep
    // the other components and give the lowest one the rest of the unit length, with the sign that lowers E.
    float r2 = dot(y, y);
    if (r2 < 1.f)
    {
        int iLo = w.x == wLo ? 0 : (w.y == wLo ? 1 : 2);
        float rest = r2 - y[iLo] * y[iLo];
        y[iLo] = (beta[iLo] > 0.f ? -1.f : 1.f) * sqrt(max(0.f, 1.f - rest));
    }
    float3 n = normalize(mul(V, y));
    return k0 + dot(b, n) + dot(n, mul(Q, n));
}

// Lower bound of E(n) = k0 + b.n + n^T Q n on the unit sphere: k0 - |b| + (Gershgorin bound of Q's lowest eigenvalue).
float shQuadraticLowerBound(float k0, float3 b, float3x3 Q)
{
    float g0 = Q[0][0] - abs(Q[0][1]) - abs(Q[0][2]);
    float g1 = Q[1][1] - abs(Q[1][0]) - abs(Q[1][2]);
    float g2 = Q[2][2] - abs(Q[2][0]) - abs(Q[2][1]);
    return k0 - length(b) + min(g0, min(g1, g2));
}

// Deringing of interpolated SH at a shading point, in place (A: band weights of evaluateIrradiance): scale bands 1-2 by
// the largest a in [0, 1] that keeps the irradiance >= 0 in every channel. Returns a (1 = unchanged).
float shDeringShadingPoint(inout float3 sh[9], float3 A)
{
    float a = 1.f;
    [unroll]
    for (int ch = 0; ch < 3; ++ch)
    {
        float k0; float3 b; float3x3 Q;
        shIrradianceQuadratic(sh, ch, A, k0, b, Q);
        if (shQuadraticLowerBound(k0, b, Q) >= 0.f) continue; // cannot go negative
        float m = shQuadraticMinOnSphere(k0, b, Q);
        if (m >= 0.f) continue;
        float e0 = 0.2820947917738781f * A.x * sh[0][ch]; // band-0 (mean) irradiance
        a = min(a, e0 > 0.f ? e0 / (e0 - m) : 0.f);
    }
    if (a < 1.f)
    {
        [unroll]
        for (int j = 1; j < 9; ++j) sh[j] *= a;
    }
    return a;
}
