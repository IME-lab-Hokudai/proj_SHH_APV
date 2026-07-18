#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
// ------------------------------------------------------------------
// CPU PORT OF SHADER LOGIC
// ------------------------------------------------------------------
float3 evaluateIrradiance(float3 n, const std::vector<float3>& shCoeffs)
{
    // Real Spherical Harmonics up to band 2
    // Coordinate convention: x = forward, y = right, z = up
    // so we need to remap falcor coords (Y-up) to polar coords (Z-up)
    float x = n.z;
    float y = n.x;
    float z = n.y;

    // ============================================================
    // Real Spherical Harmonics up to l = 2
    // Convention: Y_l^m, index = l(l + 1) + m
    // Coordinate system: (x, y, z)
    // ============================================================
    //
    //  Index layout table:
    //
    //   l |   m   | index |       Expression
    //  ----|-------|--------|------------------------------------
    //   0 |   0   |   0    | Y_0^0  = 0.282095
    //   1 |  -1   |   1    | Y_1^-1 = 0.488603 * y
    //   1 |   0   |   2    | Y_1^0  = 0.488603 * z
    //   1 |   1   |   3    | Y_1^1  = 0.488603 * x
    //   2 |  -2   |   4    | Y_2^-2 = 1.092548 * xy
    //   2 |  -1   |   5    | Y_2^-1 = 1.092548 * yz
    //   2 |   0   |   6    | Y_2^0  = 0.946175 * z^2 - 0.315392
    //   2 |   1   |   7    | Y_2^1  = 1.092548 * xz
    //   2 |   2   |   8    | Y_2^2  = 0.546274 * (x^2 - y^2)
    // ============================================================

    // with Condon–Shortley phase(-1) ^ m to match Iwasaki sensei code convention
    //  Band 0
    float Y0 = 0.2820947917738781f; // l=0, m=0

    // Band 1
    float Y1 = -0.4886025119029200f * y; // l=1, m=-1
    float Y2 = 0.4886025119029200f * z;  // l=1, m=0
    float Y3 = -0.4886025119029200f * x; // l=1, m=1

    // Band 2
    float Y4 = 1.0925484305920792f * x * y;                       // l=2, m=-2
    float Y5 = -1.0925484305920792f * y * z;                      // l=2, m=-1
    float Y6 = 0.9461746957575601f * z * z - 0.3153915652525200f; // l=2, m=0
    float Y7 = -1.0925484305920792f * x * z;                      // l=2, m=1
    float Y8 = 0.546274f * (x * x - y * y);                       // l=2, m=2

    // Cosine lobe coefficients for irradiance (Peter-Pike Sloan 2002)
    float A[3] = { 3.141593f, 2.094395f, 0.785398f };

    // Dot product: coeffs · basis
    float3 result =
        shCoeffs[0] * Y0 * A[0] +
        shCoeffs[1] * Y1 * A[1] +
        shCoeffs[2] * Y2 * A[1] +
        shCoeffs[3] * Y3 * A[1] +
        shCoeffs[4] * Y4 * A[2] +
        shCoeffs[5] * Y5 * A[2] +
        shCoeffs[6] * Y6 * A[2] +
        shCoeffs[7] * Y7 * A[2] +
        shCoeffs[8] * Y8 * A[2];

    // Clamp negative lighting
    // return max(result, float3(0.f));
    return result;
}
// 1. Standard Cubic Hermite Basis
float hermite1D(float t, float v0, float s0, float v1, float s1)
{
    float t2 = t * t;
    float t3 = t2 * t;
    float h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
    float h10 = t3 - 2.0f * t2 + t;
    float h01 = -2.0f * t3 + 3.0f * t2;
    float h11 = t3 - t2;
    return h00 * v0 + h10 * s0 + h01 * v1 + h11 * s1;
}
//

// Derivative of the Cubic Hermite Basis (d/dt)
float hermite1D_deriv(float t, float v0, float s0, float v1, float s1)
{
    float t2 = t * t;

    // Derivatives of basis functions
    // h00 = 2t^3 - 3t^2 + 1  ->  6t^2 - 6t
    float dh00 = 6.0f * t2 - 6.0f * t;

    // h10 = t^3 - 2t^2 + t   ->  3t^2 - 4t + 1
    float dh10 = 3.0f * t2 - 4.0f * t + 1.0f;

    // h01 = -2t^3 + 3t^2     ->  -6t^2 + 6t
    float dh01 = -6.0f * t2 + 6.0f * t;

    // h11 = t^3 - t^2        ->  3t^2 - 2t
    float dh11 = 3.0f * t2 - 2.0f * t;

    return dh00 * v0 + dh10 * s0 + dh01 * v1 + dh11 * s1;
}

// 2. Full Tricubic Hermite Interpolation on CPU
// This mirrors 'fetchHermiteInterpolatedSH' from your shader exactly.
void AdaptiveProbeVolume::interpolateHermite_CPU(
    int coarseProbeIdx,
    float3 pos,
    std::vector<float3>& outCoeffs,
    std::vector<GradSHCoeff>& outGrads) // We must also interpolate gradients!
{
    const auto& p = mProbes[coarseProbeIdx];
    float3 size = p.maxPoint - p.minPoint;
    float3 t = (pos - p.minPoint) / size;

    // Clamp t
    t.x = std::max(0.0f, std::min(1.0f, t.x));
    t.y = std::max(0.0f, std::min(1.0f, t.y));
    t.z = std::max(0.0f, std::min(1.0f, t.z));

    int cIdx[8];
    for (int k = 0; k < 8; ++k) cIdx[k] = p.corners[k];

    size_t numBands = mCorners[cIdx[0]].shCoeffs.size();
    outCoeffs.resize(numBands);
    outGrads.resize(numBands);

    for (size_t band = 0; band < numBands; ++band)
    {
        // 1. Fetch Inputs (Same as before)
        float3 v[8];
        float3 g_FalcorX[8], g_FalcorY[8], g_FalcorZ[8];

        for (int i = 0; i < 8; ++i) {
            const auto& c = mCorners[cIdx[i]];
            v[i] = c.shCoeffs[band];
            if (band < c.shGradients.size()) {
                // Swizzle: FalcorX=Y(.y), FalcorY=Z(.z), FalcorZ=X(.x)
                g_FalcorX[i] = float3(c.shGradients[band].r.y, c.shGradients[band].g.y, c.shGradients[band].b.y);
                g_FalcorY[i] = float3(c.shGradients[band].r.z, c.shGradients[band].g.z, c.shGradients[band].b.z);
                g_FalcorZ[i] = float3(c.shGradients[band].r.x, c.shGradients[band].g.x, c.shGradients[band].b.x);
            }
            else {
                g_FalcorX[i] = g_FalcorY[i] = g_FalcorZ[i] = float3(0.0f);
            }
        }

        // --- INTERPOLATION PASSES ---

        // PASS 1: Z-Axis (Interpolate Value AND Partial Derivatives)
        // We need q (value) and dq_dz (derivative)
        float3 q[4], dq_dz[4];
        float3 q_dX[4], q_dY[4]; // Linear interp of transverse gradients

        for (int i = 0; i < 4; ++i) {
            int i0 = 2 * i; int i1 = 2 * i + 1;

            // Value
            q[i] = float3(
                hermite1D(t.z, v[i0].x, g_FalcorZ[i0].x * size.z, v[i1].x, g_FalcorZ[i1].x * size.z),
                hermite1D(t.z, v[i0].y, g_FalcorZ[i0].y * size.z, v[i1].y, g_FalcorZ[i1].y * size.z),
                hermite1D(t.z, v[i0].z, g_FalcorZ[i0].z * size.z, v[i1].z, g_FalcorZ[i1].z * size.z)
            );

            // Derivative d/dz (Note: divide by size.z to get back to World Space!)
            dq_dz[i] = float3(
                hermite1D_deriv(t.z, v[i0].x, g_FalcorZ[i0].x * size.z, v[i1].x, g_FalcorZ[i1].x * size.z),
                hermite1D_deriv(t.z, v[i0].y, g_FalcorZ[i0].y * size.z, v[i1].y, g_FalcorZ[i1].y * size.z),
                hermite1D_deriv(t.z, v[i0].z, g_FalcorZ[i0].z * size.z, v[i1].z, g_FalcorZ[i1].z * size.z)
            ) / size.z;

            // Transverse Gradients (Linear approximation is enough for cross-terms)
            q_dX[i] = lerp(g_FalcorX[i0], g_FalcorX[i1], t.z);
            q_dY[i] = lerp(g_FalcorY[i0], g_FalcorY[i1], t.z);
        }

        // PASS 2: Y-Axis
        float3 r[2], dr_dy[2], dr_dz[2];
        float3 r_dX[2];

        for (int i = 0; i < 2; ++i) {
            int i0 = 2 * i; int i1 = 2 * i + 1;

            // Value
            r[i] = float3(
                hermite1D(t.y, q[i0].x, q_dY[i0].x * size.y, q[i1].x, q_dY[i1].x * size.y),
                hermite1D(t.y, q[i0].y, q_dY[i0].y * size.y, q[i1].y, q_dY[i1].y * size.y),
                hermite1D(t.y, q[i0].z, q_dY[i0].z * size.y, q[i1].z, q_dY[i1].z * size.y)
            );

            // Derivative d/dy
            dr_dy[i] = float3(
                hermite1D_deriv(t.y, q[i0].x, q_dY[i0].x * size.y, q[i1].x, q_dY[i1].x * size.y),
                hermite1D_deriv(t.y, q[i0].y, q_dY[i0].y * size.y, q[i1].y, q_dY[i1].y * size.y),
                hermite1D_deriv(t.y, q[i0].z, q_dY[i0].z * size.y, q[i1].z, q_dY[i1].z * size.y)
            ) / size.y;

            // Propagate d/dz (Linear Interp of previous stage)
            dr_dz[i] = lerp(dq_dz[i0], dq_dz[i1], t.y);

            // Transverse X
            r_dX[i] = lerp(q_dX[i0], q_dX[i1], t.y);
        }

        // PASS 3: X-Axis (Final Gather)

        // Final Value
        outCoeffs[band] = float3(
            hermite1D(t.x, r[0].x, r_dX[0].x * size.x, r[1].x, r_dX[1].x * size.x),
            hermite1D(t.x, r[0].y, r_dX[0].y * size.y, r[1].y, r_dX[1].y * size.x),
            hermite1D(t.x, r[0].z, r_dX[0].z * size.z, r[1].z, r_dX[1].z * size.x)
        );

        // Final Gradients
        // d/dx comes from Hermite deriv of this stage
        float3 gradX = float3(
            hermite1D_deriv(t.x, r[0].x, r_dX[0].x * size.x, r[1].x, r_dX[1].x * size.x),
            hermite1D_deriv(t.x, r[0].y, r_dX[0].y * size.y, r[1].y, r_dX[1].y * size.x),
            hermite1D_deriv(t.x, r[0].z, r_dX[0].z * size.z, r[1].z, r_dX[1].z * size.x)
        ) / size.x;

        // d/dy and d/dz come from linear interpolation of previous stages
        float3 gradY = lerp(dr_dy[0], dr_dy[1], t.x);
        float3 gradZ = lerp(dr_dz[0], dr_dz[1], t.x);

        // --- MAP BACK TO YOUR SWIZZLE ---
        // Your shader expects: 
        // .x = Falcor Z (gradZ)
        // .y = Falcor X (gradX)
        // .z = Falcor Y (gradY)
        outGrads[band].r = float3(gradZ.x, gradX.x, gradY.x);
        outGrads[band].g = float3(gradZ.y, gradX.y, gradY.y);
        outGrads[band].b = float3(gradZ.z, gradX.z, gradY.z);
    }
}

void AdaptiveProbeVolume::interpolateHermite_CPU(
    int coarseProbeIdx,
    float3 pos,
    std::vector<float3>& outCoeffs) const
{
    const auto& p = mProbes[coarseProbeIdx];

    float3 size = p.maxPoint - p.minPoint;
    float3 t = (pos - p.minPoint) / size;

    t.x = std::max(0.0f, std::min(1.0f, t.x));
    t.y = std::max(0.0f, std::min(1.0f, t.y));
    t.z = std::max(0.0f, std::min(1.0f, t.z));

    int cIdx[8];
    for (int k = 0; k < 8; ++k)
    {
        cIdx[k] = p.corners[k];
    }

    if (cIdx[0] < 0 || cIdx[0] >= int(mCorners.size()))
    {
        outCoeffs.clear();
        return;
    }

    const size_t numBands = mCorners[cIdx[0]].shCoeffs.size();

    outCoeffs.clear();
    outCoeffs.resize(numBands, float3(0.0f));

    for (size_t band = 0; band < numBands; ++band)
    {
        float3 v[8];

        // Gradients in Falcor/world axes.
        float3 gX[8];
        float3 gY[8];
        float3 gZ[8];

        for (int i = 0; i < 8; ++i)
        {
            if (cIdx[i] < 0 || cIdx[i] >= int(mCorners.size()))
            {
                outCoeffs.clear();
                return;
            }

            const Corner& c = mCorners[cIdx[i]];

            if (band >= c.shCoeffs.size())
            {
                outCoeffs.clear();
                return;
            }

            v[i] = c.shCoeffs[band];

            if (band < c.shGradients.size())
            {
                // Stored gradient convention:
                // .x = world Z
                // .y = world X
                // .z = world Y
                //
                // Therefore:
                // world X = .y
                // world Y = .z
                // world Z = .x
                gX[i] = float3(
                    c.shGradients[band].r.y,
                    c.shGradients[band].g.y,
                    c.shGradients[band].b.y
                );

                gY[i] = float3(
                    c.shGradients[band].r.z,
                    c.shGradients[band].g.z,
                    c.shGradients[band].b.z
                );

                gZ[i] = float3(
                    c.shGradients[band].r.x,
                    c.shGradients[band].g.x,
                    c.shGradients[band].b.x
                );
            }
            else
            {
                gX[i] = float3(0.0f);
                gY[i] = float3(0.0f);
                gZ[i] = float3(0.0f);
            }
        }

        // ------------------------------------------------------------
        // Corner order:
        // 0 = (0,0,0)
        // 1 = (0,0,1)
        // 2 = (0,1,0)
        // 3 = (0,1,1)
        // 4 = (1,0,0)
        // 5 = (1,0,1)
        // 6 = (1,1,0)
        // 7 = (1,1,1)
        //
        // Same separable value-only Hermite logic as uniform version.
        // ------------------------------------------------------------

        // Pass 1: Hermite along Z.
        float3 q[4];
        float3 q_dX[4];
        float3 q_dY[4];

        for (int i = 0; i < 4; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            q[i] = float3(
                hermite1D(t.z, v[i0].x, gZ[i0].x * size.z, v[i1].x, gZ[i1].x * size.z),
                hermite1D(t.z, v[i0].y, gZ[i0].y * size.z, v[i1].y, gZ[i1].y * size.z),
                hermite1D(t.z, v[i0].z, gZ[i0].z * size.z, v[i1].z, gZ[i1].z * size.z)
            );

            // Propagate transverse gradients by linear interpolation,
            // same as the uniform CPU version.
            q_dX[i] = lerp(gX[i0], gX[i1], t.z);
            q_dY[i] = lerp(gY[i0], gY[i1], t.z);
        }

        // Pass 2: Hermite along Y.
        float3 r[2];
        float3 r_dX[2];

        for (int i = 0; i < 2; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            r[i] = float3(
                hermite1D(t.y, q[i0].x, q_dY[i0].x * size.y, q[i1].x, q_dY[i1].x * size.y),
                hermite1D(t.y, q[i0].y, q_dY[i0].y * size.y, q[i1].y, q_dY[i1].y * size.y),
                hermite1D(t.y, q[i0].z, q_dY[i0].z * size.y, q[i1].z, q_dY[i1].z * size.y)
            );

            // Propagate remaining transverse X gradient.
            r_dX[i] = lerp(q_dX[i0], q_dX[i1], t.y);
        }

        // Pass 3: Hermite along X.
        outCoeffs[band] = float3(
            hermite1D(t.x, r[0].x, r_dX[0].x * size.x, r[1].x, r_dX[1].x * size.x),
            hermite1D(t.x, r[0].y, r_dX[0].y * size.x, r[1].y, r_dX[1].y * size.x),
            hermite1D(t.x, r[0].z, r_dX[0].z * size.x, r[1].z, r_dX[1].z * size.x)
        );
    }
}

void AdaptiveProbeVolume::constrainHangingNodesHermite()
{
    const int FACE_CORNERS[6][4] = {
            {0, 1, 2, 3}, // Face 0: -X (Left)  <-- WAS {0, 2, 4, 6} (WRONG: This was Back)
            {4, 5, 6, 7}, // Face 1: +X (Right) <-- WAS {1, 3, 5, 7} (WRONG: This was Front)
            {0, 1, 4, 5}, // Face 2: -Y (Down)  (Correct)
            {2, 3, 6, 7}, // Face 3: +Y (Up)    (Correct)
            {0, 2, 4, 6}, // Face 4: -Z (Back)  <-- WAS {0, 1, 2, 3} (WRONG: This was Left)
            {1, 3, 5, 7}  // Face 5: +Z (Front) <-- WAS {4, 5, 6, 7} (WRONG: This was Right)
    };

    for (size_t i = 0; i < mProbes.size(); ++i)
    {
        Probe& p = mProbes[i];
        if (!p.isLeaf) continue;

        for (int face = 0; face < 6; ++face)
        {
            int neighborIdx = p.coarseNeighbors[face];
            if (neighborIdx != -1)
            {
                Probe& coarseProbe = mProbes[neighborIdx];

                for (int k = 0; k < 4; ++k)
                {
                    int localCornerIdx = FACE_CORNERS[face][k];
                    int globalCornerIdx = p.corners[localCornerIdx];
                    Corner& c = mCorners[globalCornerIdx];

                    // 1. Calculate Hermite-Interpolated Value from Coarse Neighbor
                    std::vector<float3> constrainedCoeffs;
                    std::vector<GradSHCoeff> constrainedGrads; // We only overwrite coeffs for now

                    interpolateHermite_CPU(neighborIdx, c.position, constrainedCoeffs, constrainedGrads);

                    // 2. Overwrite ONLY the Coefficients (Values)
                    // We keep the Ray-Traced Gradients because they contain the "High Freq" detail
                    // of the fine voxel, but we force the "DC offset" to match the neighbor.
                    c.shCoeffs = constrainedCoeffs;
                    c.shGradients = constrainedGrads;
                }
            }
        }
    }
}
static float smoothstep01(float edge0, float edge1, float x)
{
    float denom = std::max(edge1 - edge0, 1e-8f);
    float t = (x - edge0) / denom;

    t = std::max(0.0f, std::min(1.0f, t));

    return t * t * (3.0f - 2.0f * t);
}

static std::string replaceExtensionOrAppend(
    const std::string& baseName,
    const std::string& suffixWithExtension
)
{
    std::string result = baseName;

    size_t dotPos = result.find_last_of('.');

    if (dotPos != std::string::npos)
    {
        result.replace(
            dotPos,
            std::string::npos,
            suffixWithExtension
        );
    }
    else
    {
        result += suffixWithExtension;
    }

    return result;
}

static double percentileForPaper(std::vector<float> values, double percentile)
{
    if (values.empty())
        return 0.0;

    std::sort(values.begin(), values.end());

    double pos = percentile * 0.01 * double(values.size() - 1);
    size_t i0 = size_t(std::floor(pos));
    size_t i1 = std::min(i0 + 1, values.size() - 1);

    double t = pos - double(i0);

    return double(values[i0]) * (1.0 - t) + double(values[i1]) * t;
}

static double percentForPaper(uint64_t n, uint64_t d)
{
    return d > 0 ? 100.0 * double(n) / double(d) : 0.0;
}

// ------------------------------------------------------------------
// Helper: Analytic Eigenvalues for 3x3 Symmetric Matrix
// ------------------------------------------------------------------
static void computeEigenvalues3x3(const float3x3& A, float& e1, float& e2, float& e3)
{
    const float ONE_THIRD = 1.0f / 3.0f;
    float m = (A[0][0] + A[1][1] + A[2][2]) * ONE_THIRD;
    float k00 = A[0][0] - m, k11 = A[1][1] - m, k22 = A[2][2] - m;
    float k01 = A[0][1], k02 = A[0][2], k12 = A[1][2];
    float q = 0.5f * (k00 * (k11 * k22 - k12 * k12) - k01 * (k01 * k22 - k12 * k02) + k02 * (k01 * k12 - k11 * k02));
    float p = (k00 * k00 + k11 * k11 + k22 * k22 + 2.0f * (k01 * k01 + k02 * k02 + k12 * k12)) / 6.0f;

    if (p < 1e-20f) { e1 = e2 = e3 = m; return; }

    float p_sqrt = std::sqrt(p);
    float det_val = q / (p * p_sqrt);
    det_val = std::max(-1.0f, std::min(1.0f, det_val));
    float phi = ONE_THIRD * std::acos(det_val);
    float two_sqrt_p = 2.0f * p_sqrt;
    float s = std::sin(phi), c = std::cos(phi);

    e1 = m + two_sqrt_p * c;
    e2 = m - two_sqrt_p * (c * 0.5f + s * 0.8660254f);
    e3 = m - two_sqrt_p * (c * 0.5f - s * 0.8660254f);
}

static int basisIndexToBand(int basisIdx)
{
    if (basisIdx == 0) return 0;
    if (basisIdx <= 3) return 1;
    return 2;
}

static float diffuseCosineA(int l)
{
    switch (l)
    {
    case 0: return 3.141593f;
    case 1: return 2.094395f; // 2pi / 3
    case 2: return 0.785398f; // pi / 4
    default: return 0.0f;
    }
}

static float maxAbsEigenvalue3x3_Local(const float3x3& h)
{
    float e1, e2, e3;
    computeEigenvalues3x3(h, e1, e2, e3);
    return std::max({ std::abs(e1), std::abs(e2), std::abs(e3) });
}

static float computeCoeffResidualL2_SHSpace(
    const std::vector<float3>& actualCoeffs,
    const std::vector<float3>& predictedCoeffs
)
{
    const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

    const size_t basisCount = 9;

    float sumSq = 0.0f;

    for (size_t basisIdx = 0; basisIdx < basisCount; ++basisIdx)
    {
        float actualLum = dot(actualCoeffs[basisIdx], kLuma);
        float predLum = dot(predictedCoeffs[basisIdx], kLuma);

        float diff = actualLum - predLum;
        sumSq += diff * diff;
    }

    return std::sqrt(sumSq);
}

static float computeCoeffVecL2_SHSpace(const std::vector<float3>& coeffs)
{
    const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

    float sumSq = 0.0f;
    for (const auto& coeff : coeffs)
    {
        float lum = dot(coeff, kLuma);
        sumSq += lum * lum;
    }

    return std::sqrt(sumSq);
}

static constexpr float kDiffuseA[9] =
{
    3.14159265358979324f, // l = 0

    2.09439510239319549f, // l = 1
    2.09439510239319549f,
    2.09439510239319549f,

    0.78539816339744831f, // l = 2
    0.78539816339744831f,
    0.78539816339744831f,
    0.78539816339744831f,
    0.78539816339744831f
};

static void evalSH9ForNormal(float3 n, float Y[9])
{
    // Same coordinate convention as evaluateIrradiance().
    //float x = n.z;
    //float y = n.x;
    //float z = n.y;
    float x = n.x;
    float y = n.y;
    float z = n.z;

    Y[0] = 0.2820947917738781f;

    Y[1] = -0.4886025119029200f * y;
    Y[2] = 0.4886025119029200f * z;
    Y[3] = -0.4886025119029200f * x;

    Y[4] = 1.0925484305920792f * x * y;
    Y[5] = -1.0925484305920792f * y * z;
    Y[6] = 0.9461746957575601f * z * z - 0.3153915652525200f;
    Y[7] = -1.0925484305920792f * x * z;
    Y[8] = 0.5462742152960396f * (x * x - y * y);
}

static const std::vector<float3>& sampledNormalsForIrrMetric()
{
    static const std::vector<float3> normals =
    {
        float3(1,  0,  0),
        float3(-1,  0,  0),
        float3(0,  1,  0),
        float3(0, -1,  0),
        float3(0,  0,  1),
        float3(0,  0, -1),

        math::normalize(float3(1,  1,  1)),
        math::normalize(float3(1,  1, -1)),
        math::normalize(float3(1, -1,  1)),
        math::normalize(float3(1, -1, -1)),
        math::normalize(float3(-1,  1,  1)),
        math::normalize(float3(-1,  1, -1)),
        math::normalize(float3(-1, -1,  1)),
        math::normalize(float3(-1, -1, -1))
    };

    return normals;
}

static float computeLambdaVecL2_SHSpace(const std::vector<float3x3>& hessians)
{
    float sumSq = 0.0f;

    for (const auto& h : hessians)
    {
        float lambda = maxAbsEigenvalue3x3_Local(h);
        sumSq += lambda * lambda;
    }

    return std::sqrt(sumSq);
}


static float computeLambdaVecL2_IrradianceSpace(const std::vector<float3x3>& hessians)
{
    if (hessians.empty())
        return 0.0f;

    const auto& normals = sampledNormalsForIrrMetric();

    float sumLambdaSq = 0.0f;

    for (float3 n : normals)
    {
        float Y[9];
        evalSH9ForNormal(n, Y);

        float3x3 HI = float3x3::zeros();

        const size_t basisCount = std::min<size_t>(hessians.size(), 9);

        for (size_t basisIdx = 0; basisIdx < basisCount; ++basisIdx)
        {
            float w = kDiffuseA[basisIdx] * Y[basisIdx];
            HI = HI + hessians[basisIdx] * w;
        }

        float lambda = maxAbsEigenvalue3x3_Local(HI);
        sumLambdaSq += lambda * lambda;
    }

    return std::sqrt(sumLambdaSq / float(normals.size()));
}

// ------------------------------------------------------------------
// Lifecycle
// ------------------------------------------------------------------

ref<AdaptiveProbeVolume> AdaptiveProbeVolume::create(ref<Device> pDevice)
{
    return ref<AdaptiveProbeVolume>(new AdaptiveProbeVolume(pDevice));
}


AdaptiveProbeVolume::AdaptiveProbeVolume(ref<Device> pDevice) : mpDevice(pDevice) {}

// ------------------------------------------------------------------
// Build Logic
// ------------------------------------------------------------------

void AdaptiveProbeVolume::startBuild(const ref<Scene>& pScene, float errorThreshold, bool useRelativeError)
{
    mProbes.clear();
    mCorners.clear();
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();
    mCornerLookup.clear();
    mCurrentThreshold = errorThreshold;
    mUseRelativeError = useRelativeError;
    resetResidualPaperStats();
    auto bounds = pScene->getSceneBounds();

    float boundsScale = 0.96f;

    // Scale the scene bounds about its center.
    float3 center = 0.5f * (bounds.minPoint + bounds.maxPoint);
    float3 halfExtent = 0.5f * (bounds.maxPoint - bounds.minPoint);
    halfExtent *= boundsScale;

    float3 scaledMin = center - halfExtent;
    float3 scaledMax = center + halfExtent;

    // 1. Create Root Probe
    //Probe root;
    //root.level = 0;
    //root.minPoint = bounds.minPoint;
    //root.maxPoint = bounds.maxPoint;

    // 1. Create Root Probe
    Probe root;
    root.level = 0;
    root.minPoint = scaledMin;
    root.maxPoint = scaledMax;

    float3 size = root.maxPoint - root.minPoint;

    // 2. Create the initial 8 Corners
    for (int i = 0; i < 8; ++i)
    {
        float3 offset = float3((i & 4) ? size.x : 0, (i & 2) ? size.y : 0, (i & 1) ? size.z : 0);

        Corner c;
        c.position = root.minPoint + offset;

        mCorners.push_back(c);
        int idx = (int)mCorners.size() - 1;

        root.corners[i] = idx;
        mPendingNewCorners.push_back(idx);
    }

    mProbes.push_back(root);
    mProbesPendingCheck.push_back(0);
}

void AdaptiveProbeVolume::getPendingPositions(std::vector<float3>& positions) const
{
    positions.clear();
    positions.reserve(mPendingNewCorners.size());
    for (int idx : mPendingNewCorners)
    {
        positions.push_back(mCorners[idx].position);
    }
}

void AdaptiveProbeVolume::setCornerData(
    uint32_t batchIndex,
    const std::vector<float3>& coeffs,
    const std::vector<GradSHCoeff>& grads,
    const std::vector<float3x3>& hessians
    // const std::vector<float>& distMeans,
    // const std::vector<float>& distMeanSqs
)
{
    // 1. Locate corner.
    if (batchIndex >= mPendingNewCorners.size()) return;

    int cornerIdx = mPendingNewCorners[batchIndex];
    Corner& c = mCorners[cornerIdx];

    // -----------------------------------------------------------
    // PART A: Runtime data
    // -----------------------------------------------------------
    c.shCoeffs = coeffs;
    c.shGradients = grads;
    c.shHessians = hessians;
    c.maxLambdaVecL2 = computeLambdaVecL2_SHSpace(hessians);

    // -----------------------------------------------------------
    // Residual scale correction + paper statistics
    // -----------------------------------------------------------
    c.residualObservedError = 0.0f;
    c.residualRatio = 1.0f;
    c.residualCorrectionScale = 1.0f;
    c.residualCorrectionValid = false;
    bool hasResidualParent =
        c.residualParentProbeIdx >= 0 &&
        c.residualParentProbeIdx < int(mProbes.size());

    bool predictedIsValid =
        c.residualPredictedError > mResidualCorrectionEps;

    bool validResidual = false;

    float observed = 0.0f;
    float predicted = 0.0f;
    float rho = 1.0f;
    float scale = 1.0f;

    int parentLevel = -1;

    if (hasResidualParent)
    {
        parentLevel = mProbes[c.residualParentProbeIdx].level;
    }

    float logRatioAbs = 0.0f;

    // For statistics:
    // disagreementConfidence = raw ratio gain before signal damping.
    // signalConfidence       = confidence that observed/predicted signal is not tiny.
    // residualConfidence     = final gain used in scale.
    float disagreementConfidence = 0.0f;
    float signalConfidence = 0.0f;
    float residualConfidence = 0.0f;

    if (mUseResidualCorrection && hasResidualParent && predictedIsValid)
    {
        std::vector<float3> predictedCoeffs;

        // Predict this new point from the parent cell.
        interpolateHermite_CPU(
            c.residualParentProbeIdx,
            c.position,
            predictedCoeffs
        );

        if (!predictedCoeffs.empty())
        {
            observed = computeCoeffResidualL2_SHSpace(
                coeffs,
                predictedCoeffs
            );

            predicted = std::max(
                c.residualPredictedError,
                mResidualCorrectionEps
            );

            rho = observed / predicted;

            float safeRho = std::max(rho, 1e-6f);

            // Symmetric disagreement measure.
            // rho = 0.5 and rho = 2.0 both mean factor-2 disagreement.
            logRatioAbs = std::abs(std::log(safeRho));

            // ------------------------------------------------------------
            // Dead-zone bounds.
            //
            // rho inside [rhoIgnoreMin, rhoIgnoreMax]:
            //     scale = 1
            //
            // rho < rhoIgnoreMin:
            //     Hessian prediction is larger than measured residual.
            //     Prune/reduce Hessian score.
            //
            // rho > rhoIgnoreMax:
            //     Hessian prediction is smaller than measured residual.
            //     Refine/increase Hessian score.
            // ------------------------------------------------------------
            float rhoIgnoreMin = std::min(
                mResidualConfidenceTau0,
                mResidualConfidenceTau1
            );

            float rhoIgnoreMax = std::max(
                mResidualConfidenceTau0,
                mResidualConfidenceTau1
            );

            rhoIgnoreMin = std::max(rhoIgnoreMin, 1e-6f);
            rhoIgnoreMax = std::max(rhoIgnoreMax, rhoIgnoreMin + 1e-6f);

            // Suppress only the extra gain above 1 when signal is tiny.
            // Base correction strength is still applied outside the dead zone.
            float signal = std::max(observed, predicted);

            signalConfidence =
                signal / (signal + mResidualConfidenceEps);

            // ------------------------------------------------------------
            // Gain limits derived from scale limits.
            //
            // scale = 1 - pruneStrength  * gain
            // scale = 1 + refineStrength * gain
            //
            // Therefore:
            // maxPruneGain  = (1 - minScale) / pruneStrength
            // maxRefineGain = (maxScale - 1) / refineStrength
            //
            // This removes the need for a separate hard-coded gainMax.
            // ------------------------------------------------------------
            const float eps = 1e-6f;

            float maxPruneGain =
                (1.0f - mResidualCorrectionMinScale) /
                std::max(mResidualPruneStrength, eps);

            float maxRefineGain =
                (mResidualCorrectionMaxScale - 1.0f) /
                std::max(mResidualRefineStrength, eps);

            // Ensure the gain cap does not break the rule:
            // outside dead zone, gain should be at least 1.
            maxPruneGain = std::max(1.0f, maxPruneGain);
            maxRefineGain = std::max(1.0f, maxRefineGain);

            float rawGain = 0.0f;
            float correctionGain = 0.0f;

            scale = 1.0f;

            if (rho < rhoIgnoreMin)
            {
                // Hessian overestimates the measured residual.
                //
                // rawGain >= 1 outside the dead zone.
                rawGain =
                    rhoIgnoreMin / std::max(rho, 1e-6f);

                rawGain =
                    std::min(rawGain, maxPruneGain);

                // Important:
                // Do not multiply the whole gain by signalConfidence.
                // Only damp the extra part above 1.
                //
                // Therefore correctionGain >= 1 outside the dead zone.
                correctionGain =
                    1.0f + signalConfidence * (rawGain - 1.0f);

                correctionGain =
                    std::max(1.0f, correctionGain);

                scale =
                    1.0f -
                    mResidualPruneStrength * correctionGain;
            }
            else if (rho > rhoIgnoreMax)
            {
                // Hessian underestimates the measured residual.
                //
                // rawGain >= 1 outside the dead zone.
                rawGain =
                    rho / std::max(rhoIgnoreMax, 1e-6f);

                rawGain =
                    std::min(rawGain, maxRefineGain);

                // Same rule:
                // base refine strength is always applied outside dead zone.
                correctionGain =
                    1.0f + signalConfidence * (rawGain - 1.0f);

                correctionGain =
                    std::max(1.0f, correctionGain);

                scale =
                    1.0f +
                    mResidualRefineStrength * correctionGain;
            }
            else
            {
                // Inside dead zone.
                rawGain = 0.0f;
                correctionGain = 0.0f;
                scale = 1.0f;
            }

            disagreementConfidence = rawGain;
            residualConfidence = correctionGain;

            // Final safety clamp.
            // With derived gain limits this should normally not be active,
            // except for numerical safety or inconsistent parameter settings.
            scale = std::max(
                mResidualCorrectionMinScale,
                std::min(mResidualCorrectionMaxScale, scale)
            );

            c.residualObservedError = observed;
            c.residualRatio = rho;
            c.residualCorrectionScale = scale;

            validResidual = true;
            c.residualCorrectionValid = true;
        }
    }

    // Record directly during construction.
    recordResidualSampleForPaper(
        parentLevel,
        hasResidualParent,
        validResidual,
        observed,
        predicted,
        rho,
        scale,
        logRatioAbs,
        disagreementConfidence,
        signalConfidence,
        residualConfidence
    );
}

void AdaptiveProbeVolume::setCornerMetricData(uint32_t batchIndex, float coeffVecL2, float maxLambdaVecL2)
{
    if (batchIndex >= mPendingNewCorners.size()) return;

    int cornerIdx = mPendingNewCorners[batchIndex];
    Corner& c = mCorners[cornerIdx];

    c.coeffVecL2 = coeffVecL2;
    c.maxLambdaVecL2 = maxLambdaVecL2;

    // Do not fill runtime data here.
    c.shCoeffs.clear();
    c.shGradients.clear();
}

void AdaptiveProbeVolume::scheduleLeafCornersForRefinementRecheck()
{
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();

    std::vector<uint8_t> used(mCorners.size(), 0);

    for (int pIdx = 0; pIdx < (int)mProbes.size(); ++pIdx)
    {
        const Probe& p = mProbes[pIdx];
        if (!p.isLeaf) continue;
        if (p.level >= mMaxLevel) continue;

        mProbesPendingCheck.push_back(pIdx);

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = p.corners[k];
            if (cIdx < 0) continue;

            if (!used[cIdx])
            {
                used[cIdx] = 1;
                mPendingNewCorners.push_back(cIdx);
            }
        }
    }
}

void AdaptiveProbeVolume::clearPendingBatch()
{
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();
}

void AdaptiveProbeVolume::scheduleAllLeafCornersForRuntimeBake()
{
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();

    std::vector<uint8_t> used(mCorners.size(), 0);

    for (int pIdx = 0; pIdx < (int)mProbes.size(); ++pIdx)
    {
        const Probe& p = mProbes[pIdx];
        if (!p.isLeaf) continue;

        // Include all leaves so runtime data is recomputed,
        // including max-level leaves.
        mProbesPendingCheck.push_back(pIdx);

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = p.corners[k];
            if (cIdx < 0) continue;

            if (!used[cIdx])
            {
                used[cIdx] = 1;
                mPendingNewCorners.push_back(cIdx);
            }
        }
    }
}

void AdaptiveProbeVolume::setCornerRuntimeData(uint32_t batchIndex, const std::vector<float3>& coeffs, const std::vector<GradSHCoeff>& grads)
{
    if (batchIndex >= mPendingNewCorners.size()) return;

    int cornerIdx = mPendingNewCorners[batchIndex];
    Corner& c = mCorners[cornerIdx];

    // Runtime data for Hermite interpolation.
    c.shCoeffs = coeffs;
    c.shGradients = grads;

    // Validity / brightness only.
    // Do NOT recompute Hessian here. Final bake does not refine topology.
    //const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

    //float sumSqL = 0.0f;
    //for (const auto& coeff : coeffs)
    //{
    //    float lum = dot(coeff, kLuma);
    //    sumSqL += lum * lum;
    //}

    //c.coeffVecL2 = std::sqrt(sumSqL);
    //c.isValid = (c.coeffVecL2 >= 0.0001f);

    // Preserve c.maxLambdaVecL2 from previous metric pass.
}



void AdaptiveProbeVolume::finishBatch()
{
    // Clear pending list; we will fill it with the NEXT generation's corners
    mPendingNewCorners.clear();

    std::vector<int> nextProbesToCheck;

    for (int probeIdx : mProbesPendingCheck)
    {
        if (mProbes[probeIdx].level >= mMaxLevel) continue;

        // --------------------------------------------------
        // Step 4: Calculate Error for this Probe
        // --------------------------------------------------
        float3 diag = mProbes[probeIdx].maxPoint - mProbes[probeIdx].minPoint;
        float distSq = dot(diag, diag);

        //float totalError = 0.0f;
        //for (int k = 0; k < 8; ++k)
        //{
        //    int cIdx = mProbes[probeIdx].corners[k];
        //    const Corner& c = mCorners[cIdx];
        //    // E_abs = 0.5 * Lambda * dx^2
        //    float e = (c.maxLambdaVecL2 * distSq) * 0.5f;

        //    totalError += e;
        //}
        //float avgProbeError = 0.0f;
        //avgProbeError = totalError / 8.0f;
        //float totalCorrectedError = 0.0f;
        //float totalHessianPredictedError = 0.0f;

        //for (int k = 0; k < 8; ++k)
        //{
        //    int cIdx = mProbes[probeIdx].corners[k];
        //    const Corner& c = mCorners[cIdx];

        //    // Original Hessian prediction:
        //    // E_abs = 0.5 * ||lambda|| * ||Delta x||^2
        //    float hessianError = 0.5f * c.maxLambdaVecL2 * distSq;

        //    // Store the uncorrected prediction separately.
        //    // This is used as parent predicted error for newly created child points.
        //    totalHessianPredictedError += hessianError;

        //    float correctedError = hessianError;

        //    if (mUseResidualCorrection)
        //    {
        //        correctedError *= c.residualCorrectionScale;
        //    }

        //    //if (mUseRelativeError)
        //    //{
        //    //    float norm = std::max(c.coeffVecL2, 1e-5f);
        //    //    correctedError /= norm;
        //    //}

        //    totalCorrectedError += correctedError;
        //}

        //float avgProbeError = totalCorrectedError / 8.0f;

        //// This is the parent Hessian prediction before residual scaling.
        //// New child points will use this as residualPredictedError.
        //float avgHessianPredictedError = totalHessianPredictedError / 8.0f;
        float totalHessianPredictedError = 0.0f;

        float maxResidualScale = 0.0f;
        bool hasValidResidualScale = false;

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = mProbes[probeIdx].corners[k];
            const Corner& c = mCorners[cIdx];

            // Original Hessian prediction:
            // E_abs = 0.5 * ||lambda|| * ||Delta x||^2
            float hessianError =
                0.5f *
                c.maxLambdaVecL2 *
                distSq;

            totalHessianPredictedError +=
                hessianError;

            // Use only spawned/evaluated residual corners.
            // Inherited/root corners may have scale = 1 only because they have
            // no residual observation, not because they are actually neutral.
            if (mUseResidualCorrection &&
                c.residualCorrectionValid)
            {
                maxResidualScale =
                    std::max(
                        maxResidualScale,
                        c.residualCorrectionScale
                    );

                hasValidResidualScale = true;
            }
        }

        float avgHessianPredictedError =
            totalHessianPredictedError / 8.0f;

        float cellResidualScale =
            hasValidResidualScale
            ? maxResidualScale
            : 1.0f;

        float avgProbeError =
            avgHessianPredictedError *
            cellResidualScale;
        recordResidualDecisionForPaper(
            mProbes[probeIdx].level,
            avgHessianPredictedError,
            avgProbeError
        );
        //float maxProbeError = 0.0f;
        //for (int k = 0; k < 8; ++k)
        //{
        //    int cIdx = mProbes[probeIdx].corners[k];
        //    const Corner& c = mCorners[cIdx];

        //    // E_abs = 0.5 * Lambda * dx^2
        //    float e = 0.5f * c.maxLambdaVecL2 * distSq;

        //    //if (mUseRelativeError)
        //    //{
        //    //    float norm = std::max(c.coeffVecL2, 1e-5f);
        //    //    e /= norm;
        //    //}

        //    maxProbeError = std::max(maxProbeError, e);
        //}

        // --------------------------------------------------
        // Step 5: Subdivide
        // --------------------------------------------------
        if (avgProbeError > mCurrentThreshold)
        //if (maxProbeError > mCurrentThreshold)
        {
            mProbes[probeIdx].isLeaf = false;

            float3 minP = mProbes[probeIdx].minPoint;
            float3 maxP = mProbes[probeIdx].maxPoint;
            float3 center = (minP + maxP) * 0.5f;
            // 1. Setup the Lookup Logic ---------------------------------------------
            float quantizationScale = 10000.0f; // Precision ~0.1mm

            auto getCornerKey = [&](float3 p) -> CornerKey {
                return {
                    (int)(floor(p.x * quantizationScale + 0.5f)),
                    (int)(floor(p.y * quantizationScale + 0.5f)),
                    (int)(floor(p.z * quantizationScale + 0.5f))
                };
                };
            // Helper: Create new corner and add to pending list
            /*auto addNewCorner = [&](float3 pos) -> int {
                Corner c; c.position = pos;
                mCorners.push_back(c);
                int idx = (int)mCorners.size() - 1;
                mPendingNewCorners.push_back(idx);
                return idx;
                };*/
            auto addNewCorner = [&](float3 pos) -> int {
                CornerKey key = getCornerKey(pos);

                // CHECK: Does this corner already exist?
                auto it = mCornerLookup.find(key);
                if (it != mCornerLookup.end()) {
                    // YES: Return existing index.
                    // DO NOT add to mPendingNewCorners (It's already been processed!)
                    return it->second;
                }

                // NO: It's new. Create it.
                Corner c;
                c.position = pos;
                // ------------------------------------------------------------
                // Residual correction metadata
                // ------------------------------------------------------------
                // This point was created because parent probeIdx was subdivided.
                // After tracing the point, setCornerData() will compare:
                //
                // actual SH at this point
                // vs.
                // parent Hermite interpolation at this point.
                c.residualParentProbeIdx = probeIdx;

                // Parent Hessian-predicted error before residual scaling.
                //c.residualPredictedError = avgHessianPredictedError;
                c.residualPredictedError = computeParentPointHessianPrediction(
                    probeIdx,
                    pos
                );

                c.residualObservedError = 0.0f;
                c.residualRatio = 1.0f;
                c.residualCorrectionScale = 1.0f;

                mCorners.push_back(c);
                int idx = (int)mCorners.size() - 1;

                // Register in map
                mCornerLookup[key] = idx;

                // Schedule for Ray Tracing
                mPendingNewCorners.push_back(idx);

                return idx;
                };

            // A. Generate 19 NEW Corners
            int c_center = addNewCorner(center);

            int c_fX0 = addNewCorner({ minP.x, center.y, center.z }); int c_fX1 = addNewCorner({ maxP.x, center.y, center.z });
            int c_fY0 = addNewCorner({ center.x, minP.y, center.z }); int c_fY1 = addNewCorner({ center.x, maxP.y, center.z });
            int c_fZ0 = addNewCorner({ center.x, center.y, minP.z }); int c_fZ1 = addNewCorner({ center.x, center.y, maxP.z });

            int c_eX_Y0Z0 = addNewCorner({ center.x, minP.y, minP.z }); int c_eX_Y1Z0 = addNewCorner({ center.x, maxP.y, minP.z });
            int c_eX_Y0Z1 = addNewCorner({ center.x, minP.y, maxP.z }); int c_eX_Y1Z1 = addNewCorner({ center.x, maxP.y, maxP.z });

            int c_eY_X0Z0 = addNewCorner({ minP.x, center.y, minP.z }); int c_eY_X1Z0 = addNewCorner({ maxP.x, center.y, minP.z });
            int c_eY_X0Z1 = addNewCorner({ minP.x, center.y, maxP.z }); int c_eY_X1Z1 = addNewCorner({ maxP.x, center.y, maxP.z });

            int c_eZ_X0Y0 = addNewCorner({ minP.x, minP.y, center.z }); int c_eZ_X1Y0 = addNewCorner({ maxP.x, minP.y, center.z });
            int c_eZ_X0Y1 = addNewCorner({ minP.x, maxP.y, center.z }); int c_eZ_X1Y1 = addNewCorner({ maxP.x, maxP.y, center.z });

            // B. Retrieve 8 OLD Corners
            const int* P = mProbes[probeIdx].corners;

            // C. Map to Grid
            int grid[3][3][3];
            grid[0][0][0] = P[0]; grid[0][0][2] = P[1]; grid[0][2][0] = P[2]; grid[0][2][2] = P[3];
            grid[2][0][0] = P[4]; grid[2][0][2] = P[5]; grid[2][2][0] = P[6]; grid[2][2][2] = P[7];
            grid[1][1][1] = c_center;
            grid[0][1][1] = c_fX0; grid[2][1][1] = c_fX1;
            grid[1][0][1] = c_fY0; grid[1][2][1] = c_fY1;
            grid[1][1][0] = c_fZ0; grid[1][1][2] = c_fZ1;
            grid[1][0][0] = c_eX_Y0Z0; grid[1][2][0] = c_eX_Y1Z0; grid[1][0][2] = c_eX_Y0Z1; grid[1][2][2] = c_eX_Y1Z1;
            grid[0][1][0] = c_eY_X0Z0; grid[2][1][0] = c_eY_X1Z0; grid[0][1][2] = c_eY_X0Z1; grid[2][1][2] = c_eY_X1Z1;
            grid[0][0][1] = c_eZ_X0Y0; grid[2][0][1] = c_eZ_X1Y0; grid[0][2][1] = c_eZ_X0Y1; grid[2][2][1] = c_eZ_X1Y1;

            // D. Create 8 Child Probes
            for (int k = 0; k < 8; ++k)
            {
                Probe child;
                child.level = mProbes[probeIdx].level + 1;

                float3 childSize = (maxP - minP) * 0.5f;
                float3 offset = float3((k & 4) ? childSize.x : 0, (k & 2) ? childSize.y : 0, (k & 1) ? childSize.z : 0);
                child.minPoint = minP + offset;
                child.maxPoint = child.minPoint + childSize;

                int startX = (k & 4) ? 1 : 0;
                int startY = (k & 2) ? 1 : 0;
                int startZ = (k & 1) ? 1 : 0;

                for (int c = 0; c < 8; ++c)
                {
                    int dx = (c & 4) ? 1 : 0;
                    int dy = (c & 2) ? 1 : 0;
                    int dz = (c & 1) ? 1 : 0;
                    child.corners[c] = grid[startX + dx][startY + dy][startZ + dz];
                }

                mProbes.push_back(child);
                mProbes[probeIdx].children[k] = (int)mProbes.size() - 1;
                nextProbesToCheck.push_back((int)mProbes.size() - 1);
            }
        }
    }
    mProbesPendingCheck = nextProbesToCheck;
}
// ------------------------------------------------------------------
// Neighbor Computation (Fixes T-Junctions)
// ------------------------------------------------------------------
void AdaptiveProbeVolume::computeNeighbors()
{
    // Directions corresponding to the 6 faces: -X, +X, -Y, +Y, -Z, +Z
    const float3 DIRECTIONS[6] = {
        float3(-1, 0, 0), float3(1, 0, 0),
        float3(0, -1, 0), float3(0, 1, 0),
        float3(0, 0, -1), float3(0, 0, 1)
    };

    for (size_t i = 0; i < mProbes.size(); ++i)
    {
        Probe& p = mProbes[i];

        // Initialize neighbors to -1 (No neighbor / No blend needed)
        for (int k = 0; k < 6; ++k) p.coarseNeighbors[k] = -1;

        // We only care about LEAF nodes
        if (p.children[0] != -1) continue;

        float3 size = p.maxPoint - p.minPoint;
        float3 center = (p.minPoint + p.maxPoint) * 0.5f;
        int totalCoarseNeighbors = 0;
        // Check all 6 faces
        for (int face = 0; face < 6; ++face)
        {
            // Step 5% of the size past the wall to ensure we enter the neighbor's voxel
            float3 dir = DIRECTIONS[face];
            float3 wallPos = center + dir * (size * 0.5f);
            float3 nudge = dir * (size * 0.1f); // Push 10% outwards
            float3 checkPos = wallPos + nudge;

            int neighborIdx = traverseOctreeCPU(checkPos);

            // If we found a valid neighbor (that is not ourselves)
            if (neighborIdx != -1 && neighborIdx != (int)i)
            {
                const Probe& neighbor = mProbes[neighborIdx];
                float3 neighborSize = neighbor.maxPoint - neighbor.minPoint;

                // THE RULE: Only store if neighbor is LARGER (Coarser)
                // Use 1.1x tolerance to avoid floating point equality issues
                if (neighborSize.x > size.x * 1.1f)
                {
                    p.coarseNeighbors[face] = neighborIdx;
                    totalCoarseNeighbors++;
                }
            }
        }
        std::cout << "Found " << totalCoarseNeighbors << " coarse neighbor links." << std::endl;
    }
}
//

void AdaptiveProbeVolume::uploadToGPU()
{
    // Ensure neighbors are computed before uploading!
    //computeNeighbors();
    //constrainHangingNodes();
    //constrainHangingNodesHermite();
    // 1. Pack Probes (Tree Topology)
    std::vector<GPUProbe> gpuProbes;
    gpuProbes.reserve(mProbes.size());

    for (const auto& p : mProbes)
    {
        GPUProbe gp;
        gp.minPoint = p.minPoint;
        gp.maxPoint = p.maxPoint;
        gp.pad1 = 0; gp.pad2 = 0;

        for (int i = 0; i < 8; ++i) gp.children[i] = p.children[i];
        for (int i = 0; i < 8; ++i) gp.corners[i] = p.corners[i];

        // --------------------------------------------------------
        // NEW: Copy Coarse Neighbor Indices
        // --------------------------------------------------------
        for (int k = 0; k < 6; ++k) gp.coarseNeighbors[k] = p.coarseNeighbors[k];
        gp.pad3[0] = 0; // Padding maintenance
        gp.pad3[1] = 0;

        gpuProbes.push_back(gp);
    }

    mpProbeBuffer = mpDevice->createStructuredBuffer(
        sizeof(GPUProbe),
        (uint32_t)gpuProbes.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        gpuProbes.data()
    );

    // 2. Pack Corners (Physics Data)
    std::vector<GPUCorner> gpuCorners;
    gpuCorners.reserve(mCorners.size());

    for (const auto& c : mCorners)
    {
        GPUCorner gc;
        // Safety: ensure we don't read out of bounds if corner has fewer bands
        int numBands = std::min((int)c.shCoeffs.size(), 9);

        for (int i = 0; i < 9; ++i)
        {
            if (i < numBands)
            {
                // PACKING:
                //gc.coeffs[i] = float4(c.shCoeffs[i], c.distMean[i]);
                gc.coeffs[i] = float4(c.shCoeffs[i], 0.0f);

                // GradR RGB = Gradient X, A = Mean Squared Distance
                //gc.gradR[i] = float4(c.shGradients[i].r, c.distMeanSq[i]);
                gc.gradR[i] = float4(c.shGradients[i].r, 0.0f);
                gc.gradG[i] = float4(c.shGradients[i].g, 0.0f);
                gc.gradB[i] = float4(c.shGradients[i].b, 0.0f);
            }
            else
            {
                gc.coeffs[i] = float4(0);
                gc.gradR[i] = float4(0);
                gc.gradG[i] = float4(0);
                gc.gradB[i] = float4(0);
            }
        }
        gpuCorners.push_back(gc);
    }

    mpCornerBuffer = mpDevice->createStructuredBuffer(
        sizeof(GPUCorner),
        (uint32_t)gpuCorners.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        gpuCorners.data()
    );
}

void AdaptiveProbeVolume::printDebugInfo(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) return;

    const auto mem = calculateMemoryFootprint();

    uint leafCounter = 0;
    for (const auto& p : mProbes)
    {
        if (p.isLeaf) leafCounter++;
    }

    out << "================================================================================\n";
    out << " ADAPTIVE PROBE VOLUME HIERARCHY (SEEDED BRANCH)\n";
    out << "================================================================================\n";
    out << " Probes (total):        " << mProbes.size() << "\n";
    out << " Probes (leaf):         " << leafCounter << "\n";
    out << " Corners:               " << mCorners.size() << "\n";
    out << " Threshold:             " << mCurrentThreshold << "\n";
    out << " Max Level:             " << mMaxLevel << "\n";
    out << " Metric:                " << (mUseRelativeError ? "Relative (E_rel)" : "Absolute (E_abs)") << "\n";
    out << " Build time (ms):       " << mBuildTimeMs << "\n";
    out << "\n";

    out << " Seeded Grid Enabled:   " << (mUseSeedGrid ? "Yes" : "No") << "\n";
    if (mUseSeedGrid)
    {
        out << " Seed Min Point:        "
            << mSeedMinPoint.x << ", "
            << mSeedMinPoint.y << ", "
            << mSeedMinPoint.z << "\n";

        out << " Seed Cell Size:        "
            << mSeedCellSize.x << ", "
            << mSeedCellSize.y << ", "
            << mSeedCellSize.z << "\n";

        out << " Seed Resolution:       "
            << mSeedResolution.x << " x "
            << mSeedResolution.y << " x "
            << mSeedResolution.z << "\n";

        out << " Seed Probe Indices:    " << mSeedProbeIndices.size() << "\n";
    }

    out << "\n";
    out << " GPU Corner Bytes:      " << mem.gpuCornersBytes << "\n";
    out << " GPU Probe Bytes:       " << mem.gpuProbesBytes << "\n";
    out << " Total Bytes:           " << mem.totalBytes << "\n";
    out << "================================================================================\n\n";

    auto computeProbeError = [&](const Probe& probe, float& avgError, float& maxLambdaInProbe, float& maxCoeffVecL2)
        {
            float3 diag = probe.maxPoint - probe.minPoint;
            float distSq = dot(diag, diag);

            //float totalError = 0.0f;
            //uint validCount = 0;
            //maxLambdaInProbe = 0.0f;

            //for (int k = 0; k < 8; ++k)
            //{
            //    int cIdx = probe.corners[k];
            //    if (cIdx < 0) continue;

            //    const Corner& c = mCorners[cIdx];
            //    maxLambdaInProbe = std::max(maxLambdaInProbe, c.maxLambdaVecL2);

            //    if (!c.isValid) continue;

            //    validCount++;
            //    float e = 0.5f * c.maxLambdaVecL2 * distSq;

            //    if (mUseRelativeError)
            //    {
            //        float norm = std::max(c.coeffVecL2, 1e-5f);
            //        totalError += e / norm;
            //    }
            //    else
            //    {
            //        totalError += e;
            //    }
            //}

            //avgError = (validCount > 0) ? (totalError / float(validCount)) : 0.0f;

            for (int k = 0; k < 8; ++k)
            {
                int cIdx = probe.corners[k];
                const Corner& c = mCorners[cIdx];

                // E_abs = 0.5 * Lambda * dx^2
                float e = 0.5f * c.maxLambdaVecL2 * distSq;

                if (mUseRelativeError)
                {
                    float norm = std::max(c.coeffVecL2, 1e-5f);
                    e /= norm;
                }

                if (e > avgError) {
                    avgError = e;
                    maxLambdaInProbe =c.maxLambdaVecL2;
                    maxCoeffVecL2 = c.coeffVecL2;
                }
            }
        };

    std::function<void(int, std::string, bool)> printProbeTree =
        [&](int probeIdx, std::string prefix, bool isLast)
        {
            const Probe& probe = mProbes[probeIdx];

            float avgError = 0.0f;
            float maxLambdaInProbe = 0.0f;
            float maxCoeffVecL2 = 0.0f;
            computeProbeError(probe, avgError, maxLambdaInProbe, maxCoeffVecL2);

            float3 diag = probe.maxPoint - probe.minPoint;
            float size = std::sqrt(dot(diag, diag));

            out << prefix << (isLast ? "L-- " : "|-- ");
            out << "[Idx " << probeIdx << "] ";
            out << "[Lvl " << probe.level << "] ";

            if (!probe.isLeaf)
            {
                out << "SUBDIVIDED ";
                out << "(Err: " << std::fixed << std::setprecision(4) << avgError
                    << " > " << mCurrentThreshold << ") ";
                out << "Size: " << std::setprecision(3) << size << " ";
                out << "MaxLambdaVecL2: " << std::setprecision(4) << maxLambdaInProbe << " ";
                out << "CoeffVecL2: " << std::setprecision(4) << maxCoeffVecL2 << "\n";

                std::vector<int> childrenIdx;
                for (int i = 0; i < 8; ++i)
                {
                    if (probe.children[i] != -1) childrenIdx.push_back(probe.children[i]);
                }

                for (size_t i = 0; i < childrenIdx.size(); ++i)
                {
                    std::string newPrefix = prefix + (isLast ? "    " : "|   ");
                    printProbeTree(childrenIdx[i], newPrefix, i == childrenIdx.size() - 1);
                }
            }
            else
            {
                out << "LEAF ";
                if (probe.level >= mMaxLevel)
                    out << "(Max Depth) ";
                else
                    out << "(Err: " << std::fixed << std::setprecision(4) << avgError
                    << " <= " << mCurrentThreshold << ") ";

                out << "Size: " << std::setprecision(3) << size << " ";
                out << "MaxLambdaVecL2:  " << std::setprecision(4) << maxLambdaInProbe << " ";
                out << "CoeffVecL2: " << std::setprecision(4) << maxCoeffVecL2 << "\n";
            }
        };

    if (mUseSeedGrid && !mSeedProbeIndices.empty())
    {
        out << "TOP-LEVEL SEEDED PROBES\n";
        out << "--------------------------------------------------------------------------------\n";

        for (uint z = 0; z < mSeedResolution.z; ++z)
        {
            for (uint y = 0; y < mSeedResolution.y; ++y)
            {
                for (uint x = 0; x < mSeedResolution.x; ++x)
                {
                    uint flat = (z * mSeedResolution.y + y) * mSeedResolution.x + x;
                    if (flat >= mSeedProbeIndices.size()) continue;

                    int rootProbeIdx = mSeedProbeIndices[flat];
                    if (rootProbeIdx < 0 || rootProbeIdx >= (int)mProbes.size()) continue;

                    out << "\n";
                    out << "Seed Cell [" << x << ", " << y << ", " << z << "] -> Probe " << rootProbeIdx << "\n";
                    printProbeTree(rootProbeIdx, "", true);
                }
            }
        }
    }
    else
    {
        out << "TREE\n";
        out << "--------------------------------------------------------------------------------\n";
        if (!mProbes.empty())
            printProbeTree(0, "", true);
        else
            out << "Tree is empty.\n";
    }

    out.close();
}

//

void AdaptiveProbeVolume::saveToFile(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out)
    {
        logError("Failed to open file for saving: " + filename);
        return;
    }

    out << std::fixed << std::setprecision(8);

    const auto mem = calculateMemoryFootprint();

    // ------------------------------------------------------------------
   // Leaf-level statistics
   // ------------------------------------------------------------------
    std::vector<uint32_t> leafCountByLevel(std::max(0, mMaxLevel) + 1, 0);

    uint32_t totalLeafCount = 0;
    uint32_t maxLevelLeafCount = 0;
    uint64_t weightedLeafLevelSum = 0;

    for (const auto& p : mProbes)
    {
        if (!p.isLeaf) continue;

        totalLeafCount++;

        if (p.level >= 0 && p.level <= mMaxLevel)
        {
            leafCountByLevel[p.level]++;
            weightedLeafLevelSum += uint64_t(p.level);
        }

        if (p.level >= mMaxLevel)
        {
            maxLevelLeafCount++;
        }
    }

    double avgLeafLevel = 0.0;
    if (totalLeafCount > 0)
    {
        avgLeafLevel = double(weightedLeafLevelSum) / double(totalLeafCount);
    }

    // ------------------------------------------------------------------
    // Header
    // ------------------------------------------------------------------
    out << "ADAPTIVE_GRID_V5_SEEDED\n";

    // ------------------------------------------------------------------
    // Build settings
    // ------------------------------------------------------------------
    out << "SETTINGS\n";
    out << mCurrentThreshold << " "
        << mMaxLevel << " "
        << (mUseRelativeError ? 1 : 0) << "\n";

    // ------------------------------------------------------------------
    // Build stats
    // ------------------------------------------------------------------
    out << "STATS\n";
    out << mBuildTimeMs << "\n";

    // ------------------------------------------------------------------
    // Leaf-level statistics
    // ------------------------------------------------------------------
    out << "LEAF_LEVEL_STATS\n";
    out << totalLeafCount << " "
        << maxLevelLeafCount << " "
        << avgLeafLevel << "\n";

    out << leafCountByLevel.size() << "\n";

    for (size_t level = 0; level < leafCountByLevel.size(); ++level)
    {
        out << level << " " << leafCountByLevel[level] << "\n";
    }

    // ------------------------------------------------------------------
    // Seed grid metadata
    // ------------------------------------------------------------------
    out << "SEED_GRID\n";
    out << (mUseSeedGrid ? 1 : 0) << "\n";

    out << mSeedMinPoint.x << " "
        << mSeedMinPoint.y << " "
        << mSeedMinPoint.z << "\n";

    out << mSeedCellSize.x << " "
        << mSeedCellSize.y << " "
        << mSeedCellSize.z << "\n";

    out << mSeedResolution.x << " "
        << mSeedResolution.y << " "
        << mSeedResolution.z << "\n";

    out << mSeedProbeIndices.size() << "\n";
    for (size_t i = 0; i < mSeedProbeIndices.size(); ++i)
    {
        out << mSeedProbeIndices[i] << " ";
    }
    out << "\n";

    // ------------------------------------------------------------------
    // Memory footprint
    // ------------------------------------------------------------------
    out << "MEMORY\n";
    out << mem.totalMB << "\n";

    // ------------------------------------------------------------------
    // Corners
    // ------------------------------------------------------------------
    out << "NUM_CORNERS " << mCorners.size() << "\n";

    for (const auto& c : mCorners)
    {
        out << c.position.x << " "
            << c.position.y << " "
            << c.position.z << " ";

        out << c.maxLambdaVecL2 << " "
            << c.coeffVecL2 << " "
            << 1 << "\n";

        out << c.shCoeffs.size() << "\n";
        for (const auto& val : c.shCoeffs)
        {
            out << val.x << " "
                << val.y << " "
                << val.z << " ";
        }
        out << "\n";

        out << c.shGradients.size() << "\n";
        for (const auto& g : c.shGradients)
        {
            out << g.r.x << " " << g.r.y << " " << g.r.z << " ";
            out << g.g.x << " " << g.g.y << " " << g.g.z << " ";
            out << g.b.x << " " << g.b.y << " " << g.b.z << " ";
        }
        out << "\n";
    }

    // ------------------------------------------------------------------
    // Probes
    // ------------------------------------------------------------------
    out << "NUM_PROBES " << mProbes.size() << "\n";

    for (const auto& p : mProbes)
    {
        out << (p.isLeaf ? 1 : 0) << " "
            << p.level << "\n";

        out << p.minPoint.x << " "
            << p.minPoint.y << " "
            << p.minPoint.z << " ";

        out << p.maxPoint.x << " "
            << p.maxPoint.y << " "
            << p.maxPoint.z << "\n";

        for (int i = 0; i < 8; ++i)
        {
            out << p.corners[i] << " ";
        }
        out << "\n";

        for (int i = 0; i < 8; ++i)
        {
            out << p.children[i] << " ";
        }
        out << "\n";

        for (int i = 0; i < 6; ++i)
        {
            out << p.coarseNeighbors[i] << " ";
        }
        out << "\n";
    }

    out.close();

    logInfo("Successfully saved AdaptiveProbeVolume (V5 seeded) to " + filename);
}

void AdaptiveProbeVolume::loadFromFile(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
    {
        logError("Failed to open file for loading: " + filename);
        return;
    }

    mProbes.clear();
    mCorners.clear();
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();
    mCornerLookup.clear();

    mUseSeedGrid = false;
    mSeedMinPoint = float3(0.f);
    mSeedCellSize = float3(0.f);
    mSeedResolution = uint3(0);
    mSeedProbeIndices.clear();

    std::string header;
    in >> header;

    if (header != "ADAPTIVE_GRID_V5_SEEDED")
    {
        logError("Invalid file format or version mismatch (Expected V5 seeded): " + filename);
        return;
    }

    std::string tag;

    // ------------------------------------------------------------------
    // Settings
    // ------------------------------------------------------------------
    in >> tag;
    if (tag != "SETTINGS")
    {
        logError("Missing SETTINGS block in: " + filename);
        return;
    }

    int useRelErrInt = 0;
    in >> mCurrentThreshold >> mMaxLevel >> useRelErrInt;
    mUseRelativeError = (useRelErrInt != 0);

    // ------------------------------------------------------------------
    // Stats
    // ------------------------------------------------------------------
    in >> tag;
    if (tag != "STATS")
    {
        logError("Missing STATS block in: " + filename);
        return;
    }

    in >> mBuildTimeMs;

    // ------------------------------------------------------------------
    // Next block can be either:
    //   LEAF_LEVEL_STATS  -> new files
    //   SEED_GRID         -> old files
    // ------------------------------------------------------------------
    in >> tag;

    if (tag == "LEAF_LEVEL_STATS")
    {
        uint32_t totalLeafCountFromFile = 0;
        uint32_t maxLevelLeafCountFromFile = 0;
        double avgLeafLevelFromFile = 0.0;

        in >> totalLeafCountFromFile
            >> maxLevelLeafCountFromFile
            >> avgLeafLevelFromFile;

        size_t leafLevelStatCount = 0;
        in >> leafLevelStatCount;

        for (size_t i = 0; i < leafLevelStatCount; ++i)
        {
            uint32_t level = 0;
            uint32_t count = 0;
            in >> level >> count;
        }

        // Read next expected block.
        in >> tag;
    }

    // ------------------------------------------------------------------
    // Seed grid metadata
    // ------------------------------------------------------------------
    if (tag != "SEED_GRID")
    {
        logError("Missing SEED_GRID block in: " + filename);
        return;
    }

    int useSeedInt = 0;
    in >> useSeedInt;
    mUseSeedGrid = (useSeedInt != 0);

    in >> mSeedMinPoint.x
        >> mSeedMinPoint.y
        >> mSeedMinPoint.z;

    in >> mSeedCellSize.x
        >> mSeedCellSize.y
        >> mSeedCellSize.z;

    in >> mSeedResolution.x
        >> mSeedResolution.y
        >> mSeedResolution.z;

    size_t seedCount = 0;
    in >> seedCount;

    mSeedProbeIndices.resize(seedCount);
    for (size_t i = 0; i < seedCount; ++i)
    {
        in >> mSeedProbeIndices[i];
    }

    // ------------------------------------------------------------------
    // Memory block
    // ------------------------------------------------------------------
    in >> tag;
    if (tag != "MEMORY")
    {
        logError("Missing MEMORY block in: " + filename);
        return;
    }

    // Current save writes totalMB, so read as double.
    double totalMBFromFile = 0.0;
    in >> totalMBFromFile;

    // ------------------------------------------------------------------
    // Corners
    // ------------------------------------------------------------------
    size_t numCorners = 0;
    in >> tag >> numCorners;

    if (tag != "NUM_CORNERS")
    {
        logError("Missing NUM_CORNERS block in: " + filename);
        return;
    }

    mCorners.resize(numCorners);

    for (size_t i = 0; i < numCorners; ++i)
    {
        Corner& c = mCorners[i];

        int isValidInt = 0;

        in >> c.position.x
            >> c.position.y
            >> c.position.z;

        in >> c.maxLambdaVecL2
            >> c.coeffVecL2
            >> isValidInt;

        //c.isValid = (isValidInt != 0);

        size_t numCoeffs = 0;
        in >> numCoeffs;

        c.shCoeffs.resize(numCoeffs);
        for (size_t k = 0; k < numCoeffs; ++k)
        {
            in >> c.shCoeffs[k].x
                >> c.shCoeffs[k].y
                >> c.shCoeffs[k].z;
        }

        size_t numGrads = 0;
        in >> numGrads;

        c.shGradients.resize(numGrads);
        for (size_t k = 0; k < numGrads; ++k)
        {
            in >> c.shGradients[k].r.x
                >> c.shGradients[k].r.y
                >> c.shGradients[k].r.z;

            in >> c.shGradients[k].g.x
                >> c.shGradients[k].g.y
                >> c.shGradients[k].g.z;

            in >> c.shGradients[k].b.x
                >> c.shGradients[k].b.y
                >> c.shGradients[k].b.z;
        }

        // Rebuild lookup table for future corner sharing.
        CornerKey key{
            (int)std::floor(c.position.x * 10000.0f + 0.5f),
            (int)std::floor(c.position.y * 10000.0f + 0.5f),
            (int)std::floor(c.position.z * 10000.0f + 0.5f)
        };

        mCornerLookup[key] = (int)i;
    }

    // ------------------------------------------------------------------
    // Probes
    // ------------------------------------------------------------------
    size_t numProbes = 0;
    in >> tag >> numProbes;

    if (tag != "NUM_PROBES")
    {
        logError("Missing NUM_PROBES block in: " + filename);
        return;
    }

    mProbes.resize(numProbes);

    for (size_t i = 0; i < numProbes; ++i)
    {
        Probe& p = mProbes[i];

        int isLeafInt = 0;
        in >> isLeafInt >> p.level;
        p.isLeaf = (isLeafInt != 0);

        in >> p.minPoint.x
            >> p.minPoint.y
            >> p.minPoint.z;

        in >> p.maxPoint.x
            >> p.maxPoint.y
            >> p.maxPoint.z;

        for (int k = 0; k < 8; ++k)
        {
            in >> p.corners[k];
        }

        for (int k = 0; k < 8; ++k)
        {
            in >> p.children[k];
        }

        for (int k = 0; k < 6; ++k)
        {
            in >> p.coarseNeighbors[k];
        }
    }

    in.close();

    logInfo("Successfully loaded AdaptiveProbeVolume (V5 seeded) from " + filename);
}


void AdaptiveProbeVolume::startBuildSeeded(
    const ref<Scene>& pScene,
    uint3 seedResolution,
    float errorThreshold,
    bool useRelativeError
)
{
    resetResidualPaperStats();
    uint32_t cellsPerAxis = 1u << mMaxLevel;
    uint64_t maxNodes = 0;
    uint64_t levelCount = 1;

    for (int  l = 0; l <= mMaxLevel; ++l)
    {
        maxNodes += levelCount;
        levelCount *= 8;
    }

    uint64_t maxCorners =
        uint64_t(cellsPerAxis + 1) *
        uint64_t(cellsPerAxis + 1) *
        uint64_t(cellsPerAxis + 1);

    mProbes.reserve(size_t(maxNodes));
    mCorners.reserve(size_t(maxCorners));
    mPendingNewCorners.reserve(size_t(maxCorners));
    mProbesPendingCheck.reserve(size_t(maxNodes));
    mCornerLookup.reserve(size_t(maxCorners * 2));

    //mProbes.clear();
    //mCorners.clear();
    //mPendingNewCorners.clear();
    //mProbesPendingCheck.clear();
    //mCornerLookup.clear();

    mCurrentThreshold = errorThreshold;
    mUseRelativeError = useRelativeError;

    auto bounds = pScene->getSceneBounds();

    float boundsScale = 0.98f;
    float3 center = 0.5f * (bounds.minPoint + bounds.maxPoint);
    float3 halfExtent = 0.5f * (bounds.maxPoint - bounds.minPoint);
    halfExtent *= boundsScale;

    float3 scaledMin = center - halfExtent;
    float3 scaledMax = center + halfExtent;

    float3 totalSize = scaledMax - scaledMin;
    float3 cellSize = totalSize / float3(seedResolution);

    mUseSeedGrid = true;
    mSeedMinPoint = scaledMin;
    mSeedCellSize = cellSize;
    mSeedResolution = seedResolution;
    mSeedProbeIndices.clear();
    mSeedProbeIndices.resize(seedResolution.x * seedResolution.y * seedResolution.z);

    auto getCornerIndex = [&](uint32_t ix, uint32_t iy, uint32_t iz) -> int
        {
            float3 p = scaledMin + float3(ix, iy, iz) * cellSize;

            CornerKey key{
                (int)std::floor(p.x * 10000.0f + 0.5f),
                (int)std::floor(p.y * 10000.0f + 0.5f),
                (int)std::floor(p.z * 10000.0f + 0.5f)
            };

            auto it = mCornerLookup.find(key);
            if (it != mCornerLookup.end()) return it->second;

            Corner c;
            c.position = p;
            mCorners.push_back(c);
            int idx = (int)mCorners.size() - 1;
            mCornerLookup[key] = idx;
            mPendingNewCorners.push_back(idx);
            return idx;
        };

    for (uint32_t z = 0; z < seedResolution.z; ++z)
    {
        for (uint32_t y = 0; y < seedResolution.y; ++y)
        {
            for (uint32_t x = 0; x < seedResolution.x; ++x)
            {
                Probe p;
                p.isLeaf = true;
                p.level = 0;
                p.minPoint = scaledMin + float3(x, y, z) * cellSize;
                p.maxPoint = p.minPoint + cellSize;

                p.corners[0] = getCornerIndex(x, y, z);
                p.corners[1] = getCornerIndex(x, y, z + 1);
                p.corners[2] = getCornerIndex(x, y + 1, z);
                p.corners[3] = getCornerIndex(x, y + 1, z + 1);
                p.corners[4] = getCornerIndex(x + 1, y, z);
                p.corners[5] = getCornerIndex(x + 1, y, z + 1);
                p.corners[6] = getCornerIndex(x + 1, y + 1, z);
                p.corners[7] = getCornerIndex(x + 1, y + 1, z + 1);

                for (int i = 0; i < 6; ++i) p.coarseNeighbors[i] = -1;

                uint flat = (z * seedResolution.y + y) * seedResolution.x + x;
                mSeedProbeIndices[flat] = (int)mProbes.size(); // before push_back or after adjust accordingly

                mProbes.push_back(p);
                mProbesPendingCheck.push_back((int)mProbes.size() - 1);
            }
        }
    }
}

void AdaptiveProbeVolume::getPendingPositionsRange(
    uint32_t start,
    uint32_t count,
    std::vector<float3>& positions
) const
{
    positions.clear();
    uint32_t end = std::min(start + count, (uint32_t)mPendingNewCorners.size());
    positions.reserve(end - start);

    for (uint32_t i = start; i < end; ++i)
    {
        positions.push_back(mCorners[mPendingNewCorners[i]].position);
    }
}

void AdaptiveProbeVolume::setCornerDataRange(
    uint32_t start,
    const std::vector<std::vector<float3>>& coeffsBatch,
    const std::vector<std::vector<GradSHCoeff>>& gradsBatch,
    const std::vector<std::vector<float3x3>>& hessiansBatch
)
{
    uint32_t count = (uint32_t)coeffsBatch.size();
    for (uint32_t i = 0; i < count; ++i)
    {
        setCornerData(start + i, coeffsBatch[i], gradsBatch[i], hessiansBatch[i]);
    }
}

float3 AdaptiveProbeVolume::evaluateIrradianceHermiteCPU(float3 posW, float3 normalW) const
{
    std::vector<float3> coeffs;

    int probeIdx = traverseOctreeCPU(posW);
    if (probeIdx < 0)
        return float3(0.0f);

    interpolateHermite_CPU(probeIdx, posW, coeffs);

    if (coeffs.size() < 9)
        return float3(0.0f);

    return evaluateIrradiance(normalW, coeffs);
}

int AdaptiveProbeVolume::findSeedProbeCPU(float3 pos) const
{
    if (!mUseSeedGrid) return 0;
    if (mSeedProbeIndices.empty()) return -1;

    // Convert world position to seed-grid coordinates
    float3 rel = (pos - mSeedMinPoint) / mSeedCellSize;

    int ix = (int)std::floor(rel.x);
    int iy = (int)std::floor(rel.y);
    int iz = (int)std::floor(rel.z);

    // Bounds check in seed grid
    if (ix < 0 || iy < 0 || iz < 0 ||
        ix >= (int)mSeedResolution.x ||
        iy >= (int)mSeedResolution.y ||
        iz >= (int)mSeedResolution.z)
    {
        return -1;
    }

    uint idx = (uint(iz) * mSeedResolution.y + uint(iy)) * mSeedResolution.x + uint(ix);
    if (idx >= mSeedProbeIndices.size()) return -1;

    return mSeedProbeIndices[idx];
}

int AdaptiveProbeVolume::traverseOctreeCPU(float3 pos) const
{
    if (mProbes.empty()) return -1;

    int currentIdx = 0;

    if (mUseSeedGrid)
    {
        currentIdx = findSeedProbeCPU(pos);
        if (currentIdx < 0) return -1;
    }
    else
    {
        // Old single-root behavior
        const Probe& root = mProbes[0];
        if (pos.x < root.minPoint.x || pos.y < root.minPoint.y || pos.z < root.minPoint.z ||
            pos.x > root.maxPoint.x || pos.y > root.maxPoint.y || pos.z > root.maxPoint.z)
        {
            return -1;
        }
        currentIdx = 0;
    }

    // Local octree descent from the seed probe (or root)
    for (int i = 0; i < mMaxLevel + 8; ++i)
    {
        const Probe& p = mProbes[currentIdx];

        // Extra safety: if somehow pos is outside this probe, fail gracefully
        if (pos.x < p.minPoint.x || pos.y < p.minPoint.y || pos.z < p.minPoint.z ||
            pos.x > p.maxPoint.x || pos.y > p.maxPoint.y || pos.z > p.maxPoint.z)
        {
            return -1;
        }

        if (p.children[0] == -1)
        {
            return currentIdx;
        }

        float3 center = 0.5f * (p.minPoint + p.maxPoint);
        int octant = 0;
        if (pos.x >= center.x) octant |= 4;
        if (pos.y >= center.y) octant |= 2;
        if (pos.z >= center.z) octant |= 1;

        int childIdx = p.children[octant];
        if (childIdx == -1) return currentIdx;

        currentIdx = childIdx;
    }

    return currentIdx;
}

AdaptiveProbeVolume::MemoryFootprintInfo AdaptiveProbeVolume::calculateMemoryFootprint() const
{
    MemoryFootprintInfo info{};

    // GPU packed layout
    info.gpuProbesBytes = uint64_t(mProbes.size()) * sizeof(GPUProbe);
    info.gpuCornersBytes = uint64_t(mCorners.size()) * sizeof(GPUCorner);
    const double invMB = 1.0 / 1e6;               // decimal
    info.totalBytes =
    info.gpuCornersBytes +
    info.gpuProbesBytes;
    info.totalMB = info.totalBytes * invMB;
    return info;
}

void AdaptiveProbeVolume::printCoarseStageDebugInfo(const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out.is_open()) return;

    uint32_t leafCount = 0;
    uint32_t maxLevelLeafCount = 0;

    // --------------------------------------------------
    // Header summary
    // --------------------------------------------------
    out << "=== Coarse Stage Debug Info ===\n";
    out << "Total Probes: " << mProbes.size() << "\n";
    out << "Total Corners: " << mCorners.size() << "\n";

    for (const auto& p : mProbes)
    {
        if (!p.isLeaf) continue;
        leafCount++;
        if (p.level >= mMaxLevel)
            maxLevelLeafCount++;
    }

    out << "Leaf Probes: " << leafCount << "\n";
    out << "Max-Level Leaves: " << maxLevelLeafCount << "\n";
    out << "Pending Probes: " << mProbesPendingCheck.size() << "\n";
    out << "Pending Corners: " << mPendingNewCorners.size() << "\n";
    out << "Threshold: " << mCurrentThreshold << "\n";
    out << "====================================\n\n";

    // --------------------------------------------------
    // Per-probe dump (THIS is what you need)
    // --------------------------------------------------
    for (size_t i = 0; i < mProbes.size(); ++i)
    {
        const Probe& p = mProbes[i];

        out << "Probe " << i << "\n";
        out << "  Level: " << p.level << "\n";
        out << "  IsLeaf: " << p.isLeaf << "\n";

        float3 size = p.maxPoint - p.minPoint;
        out << "  Size: (" << size.x << ", " << size.y << ", " << size.z << ")\n";

        // ---- Critical: error metric ----
        float maxError = 0.0f;

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = p.corners[k];
            if (cIdx < 0) continue;

            const Corner& c = mCorners[cIdx];

            float e = (c.maxLambdaVecL2 * dot(size, size)) * 0.5f;

            if (mUseRelativeError)
            {
                float norm = std::max(c.coeffVecL2, 1e-5f);
                e /= norm;
            }

            maxError = std::max(maxError, e);
        }

        out << "  MaxError: " << maxError << "\n";

        // ---- Corners ----
        out << "  Corners:\n";
        for (int k = 0; k < 8; ++k)
        {
            int cIdx = p.corners[k];
            if (cIdx < 0) continue;

            const Corner& c = mCorners[cIdx];

            out << "    Corner " << k << " idx " << cIdx
                << " LambdaL2: " << c.maxLambdaVecL2
                << " CoeffL2: " << c.coeffVecL2 << "\n";
        }

        out << "\n";
    }

    out.close();
}

void AdaptiveProbeVolume::finishBatchCoarseLimited(int maxEvalLevel)
{
    mPendingNewCorners.clear();

    std::vector<int> nextProbesToCheck;

    for (int probeIdx : mProbesPendingCheck)
    {
        if (mProbes[probeIdx].level >= mMaxLevel) continue;

        float3 diag = mProbes[probeIdx].maxPoint - mProbes[probeIdx].minPoint;
        float distSq = dot(diag, diag);

        float maxProbeError = 0.0f;

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = mProbes[probeIdx].corners[k];
            const Corner& c = mCorners[cIdx];

            float e = 0.5f * c.maxLambdaVecL2 * distSq;

            if (mUseRelativeError)
            {
                float norm = std::max(c.coeffVecL2, 1e-5f);
                e /= norm;
            }

            maxProbeError = std::max(maxProbeError, e);
        }

        if (maxProbeError <= mCurrentThreshold)
            continue;

        bool enqueueChildrenForCoarse = (mProbes[probeIdx].level < maxEvalLevel);

        float3 minP = mProbes[probeIdx].minPoint;
        float3 maxP = mProbes[probeIdx].maxPoint;
        float3 center = (minP + maxP) * 0.5f;

        mProbes[probeIdx].isLeaf = false;

        float quantizationScale = 10000.0f;

        auto getCornerKey = [&](float3 p) -> CornerKey {
            return {
                (int)(floor(p.x * quantizationScale + 0.5f)),
                (int)(floor(p.y * quantizationScale + 0.5f)),
                (int)(floor(p.z * quantizationScale + 0.5f))
            };
            };

        auto addNewCorner = [&](float3 pos) -> int {
            CornerKey key = getCornerKey(pos);

            auto it = mCornerLookup.find(key);
            if (it != mCornerLookup.end())
                return it->second;

            Corner c;
            c.position = pos;
            mCorners.push_back(c);

            int idx = (int)mCorners.size() - 1;
            mCornerLookup[key] = idx;

            if (enqueueChildrenForCoarse)
                mPendingNewCorners.push_back(idx);

            return idx;
            };

        int c_center = addNewCorner(center);

        int c_fX0 = addNewCorner({ minP.x, center.y, center.z });
        int c_fX1 = addNewCorner({ maxP.x, center.y, center.z });
        int c_fY0 = addNewCorner({ center.x, minP.y, center.z });
        int c_fY1 = addNewCorner({ center.x, maxP.y, center.z });
        int c_fZ0 = addNewCorner({ center.x, center.y, minP.z });
        int c_fZ1 = addNewCorner({ center.x, center.y, maxP.z });

        int c_eX_Y0Z0 = addNewCorner({ center.x, minP.y, minP.z });
        int c_eX_Y1Z0 = addNewCorner({ center.x, maxP.y, minP.z });
        int c_eX_Y0Z1 = addNewCorner({ center.x, minP.y, maxP.z });
        int c_eX_Y1Z1 = addNewCorner({ center.x, maxP.y, maxP.z });

        int c_eY_X0Z0 = addNewCorner({ minP.x, center.y, minP.z });
        int c_eY_X1Z0 = addNewCorner({ maxP.x, center.y, minP.z });
        int c_eY_X0Z1 = addNewCorner({ minP.x, center.y, maxP.z });
        int c_eY_X1Z1 = addNewCorner({ maxP.x, center.y, maxP.z });

        int c_eZ_X0Y0 = addNewCorner({ minP.x, minP.y, center.z });
        int c_eZ_X1Y0 = addNewCorner({ maxP.x, minP.y, center.z });
        int c_eZ_X0Y1 = addNewCorner({ minP.x, maxP.y, center.z });
        int c_eZ_X1Y1 = addNewCorner({ maxP.x, maxP.y, center.z });

        const int* P = mProbes[probeIdx].corners;

        int grid[3][3][3];
        grid[0][0][0] = P[0]; grid[0][0][2] = P[1]; grid[0][2][0] = P[2]; grid[0][2][2] = P[3];
        grid[2][0][0] = P[4]; grid[2][0][2] = P[5]; grid[2][2][0] = P[6]; grid[2][2][2] = P[7];

        grid[1][1][1] = c_center;

        grid[0][1][1] = c_fX0; grid[2][1][1] = c_fX1;
        grid[1][0][1] = c_fY0; grid[1][2][1] = c_fY1;
        grid[1][1][0] = c_fZ0; grid[1][1][2] = c_fZ1;

        grid[1][0][0] = c_eX_Y0Z0; grid[1][2][0] = c_eX_Y1Z0;
        grid[1][0][2] = c_eX_Y0Z1; grid[1][2][2] = c_eX_Y1Z1;

        grid[0][1][0] = c_eY_X0Z0; grid[2][1][0] = c_eY_X1Z0;
        grid[0][1][2] = c_eY_X0Z1; grid[2][1][2] = c_eY_X1Z1;

        grid[0][0][1] = c_eZ_X0Y0; grid[2][0][1] = c_eZ_X1Y0;
        grid[0][2][1] = c_eZ_X0Y1; grid[2][2][1] = c_eZ_X1Y1;

        const int parentLevel = mProbes[probeIdx].level;

        for (int k = 0; k < 8; ++k)
        {
            Probe child;
            child.level = parentLevel + 1;

            float3 childSize = (maxP - minP) * 0.5f;
            float3 offset = float3(
                (k & 4) ? childSize.x : 0,
                (k & 2) ? childSize.y : 0,
                (k & 1) ? childSize.z : 0
            );

            child.minPoint = minP + offset;
            child.maxPoint = child.minPoint + childSize;

            int startX = (k & 4) ? 1 : 0;
            int startY = (k & 2) ? 1 : 0;
            int startZ = (k & 1) ? 1 : 0;

            for (int c = 0; c < 8; ++c)
            {
                int dx = (c & 4) ? 1 : 0;
                int dy = (c & 2) ? 1 : 0;
                int dz = (c & 1) ? 1 : 0;

                child.corners[c] = grid[startX + dx][startY + dy][startZ + dz];
            }

            mProbes.push_back(child);

            int childIdx = (int)mProbes.size() - 1;
            mProbes[probeIdx].children[k] = childIdx;

            if (enqueueChildrenForCoarse)
                nextProbesToCheck.push_back(childIdx);
        }
    }

    mProbesPendingCheck = nextProbesToCheck;
}
void AdaptiveProbeVolume::resetResidualPaperStats()
{
    mResidualPaperStats.reset();
}

void AdaptiveProbeVolume::recordResidualSampleForPaper(
    int parentLevel,
    bool hasParent,
    bool validResidual,
    float observedResidual,
    float predictedParentError,
    float rho,
    float scale,
    float logRatioAbs,
    float disagreementConfidence,
    float signalConfidence,
    float residualConfidence
)
{
    auto recordOne = [&](ResidualSamplePaperStats& s)
        {
            s.tracedCorners++;

            if (!hasParent)
            {
                s.noParentCorners++;
                return;
            }

            s.parentCorners++;

            if (!validResidual)
            {
                s.invalidResidualCorners++;
                return;
            }

            s.validResidualCorners++;

            s.observedResidual.add(observedResidual);
            s.predictedParentError.add(predictedParentError);
            s.residualRatio.add(rho);
            s.correctionScale.add(scale);

            const float scaleOneEps = 1e-5f;

            // Tolerance for "near min/max scale" statistics.
            // 1e-3 means 0.500228 counts as near 0.5.
            const float scaleLimitEps = 1e-3f;

            if (scale < 1.0f - scaleOneEps)
                s.scaleBelowOne++;
            else if (scale > 1.0f + scaleOneEps)
                s.scaleAboveOne++;
            else
                s.scaleEqualOne++;

            // Count values that are effectively at the min/max scale limit.
            // This is better interpreted as "near min/max limit", not necessarily
            // "the final clamp operation was exactly triggered".
            if (scale <= mResidualCorrectionMinScale + scaleLimitEps)
                s.scaleAtMinClamp++;

            if (scale >= mResidualCorrectionMaxScale - scaleLimitEps)
                s.scaleAtMaxClamp++;
            s.logRatioAbs.add(logRatioAbs);
            s.disagreementConfidence.add(disagreementConfidence);
            s.signalConfidence.add(signalConfidence);
            s.residualConfidence.add(residualConfidence);

            // ------------------------------------------------------------
 // Dead-zone residual statistics.
 //
 // In the new formulation:
 //     mResidualConfidenceTau0 and mResidualConfidenceTau1
 //     are rho ignore-range bounds, e.g. 0.9 and 1.1.
 //
 // rho < rhoIgnoreMin:
 //     Hessian overestimates residual -> pruning side
 //
 // rho inside [rhoIgnoreMin, rhoIgnoreMax]:
 //     residual correction ignored
 //
 // rho > rhoIgnoreMax:
 //     Hessian underestimates residual -> refinement side
 // ------------------------------------------------------------

            float rhoIgnoreMin = std::min(
                mResidualConfidenceTau0,
                mResidualConfidenceTau1
            );

            float rhoIgnoreMax = std::max(
                mResidualConfidenceTau0,
                mResidualConfidenceTau1
            );

            rhoIgnoreMin = std::max(rhoIgnoreMin, 1e-6f);
            rhoIgnoreMax = std::max(rhoIgnoreMax, rhoIgnoreMin + 1e-6f);

            if (rho < rhoIgnoreMin)
            {
                // Prune-side correction.
                s.overEstimateSignal++;
            }
            else if (rho > rhoIgnoreMax)
            {
                // Refine-side correction.
                s.underEstimateSignal++;
            }
            else
            {
                // Inside dead zone.
                // No residual correction should be applied.
                s.nearConsistentSignal++;
            }

            // Reinterpret these existing counters:
            //
            // lowConfidenceSamples:
            //     correction gain is effectively zero.
            //     Usually inside the ignore range, or suppressed by weak signal.
            //
            // highConfidenceSamples:
            //     correction gain reaches at least the base strength.
            //     That means residualConfidence >= 1.
            const float gainEps = 1e-5f;

            if (residualConfidence <= gainEps)
                s.lowConfidenceSamples++;

            if (residualConfidence >= 1.0f - gainEps)
                s.highConfidenceSamples++;
        };

    recordOne(mResidualPaperStats.allSamples);

    if (parentLevel >= 0 && parentLevel < kPaperStatsMaxLevel)
        recordOne(mResidualPaperStats.samplesByParentLevel[parentLevel]);
}

void AdaptiveProbeVolume::recordResidualDecisionForPaper(
    int cellLevel,
    float originalHessianError,
    float correctedResidualError
)
{
    auto recordOne = [&](ResidualDecisionPaperStats& d)
        {
            d.testedCells++;

            bool originalSubdivide =
                originalHessianError > mCurrentThreshold;

            bool correctedSubdivide =
                correctedResidualError > mCurrentThreshold;

            if (originalSubdivide && correctedSubdivide)
                d.bothSubdivide++;
            else if (!originalSubdivide && !correctedSubdivide)
                d.bothLeaf++;
            else if (originalSubdivide && !correctedSubdivide)
                d.prunedByResidual++;
            else if (!originalSubdivide && correctedSubdivide)
                d.addedByResidual++;

            d.originalHessianError.add(originalHessianError);
            d.correctedResidualError.add(correctedResidualError);

            float ratio = originalHessianError > 1e-8f
                ? correctedResidualError / originalHessianError
                : 1.0f;

            d.correctedOverOriginalRatio.add(ratio);
        };

    recordOne(mResidualPaperStats.allDecisions);

    if (cellLevel >= 0 && cellLevel < kPaperStatsMaxLevel)
        recordOne(mResidualPaperStats.decisionsByCellLevel[cellLevel]);
}

void AdaptiveProbeVolume::printResidualPaperStats(std::ostream& out) const
{
    auto printScalar = [&](const char* name, const PaperScalarStats& s)
        {
            out << "   " << name << "\n";
            out << "      count:       " << s.count << "\n";
            out << "      mean:        " << s.mean() << "\n";
            out << "      stddev:      " << s.stddev() << "\n";
            out << "      median:      " << percentileForPaper(s.values, 50.0) << "\n";
            out << "      p05:         " << percentileForPaper(s.values, 5.0) << "\n";
            out << "      p95:         " << percentileForPaper(s.values, 95.0) << "\n";
            out << "      min:         " << s.minOrZero() << "\n";
            out << "      max:         " << s.maxOrZero() << "\n";
        };

    auto printSample = [&](const char* title, const ResidualSamplePaperStats& s)
        {
            out << "\n";
            out << title << "\n";
            out << "--------------------------------------------------------------------------------\n";

            out << " Traced corners:              " << s.tracedCorners << "\n";
            out << " No-parent corners:           " << s.noParentCorners
                << " (" << percentForPaper(s.noParentCorners, s.tracedCorners) << "%)\n";
            out << " Parent corners:              " << s.parentCorners
                << " (" << percentForPaper(s.parentCorners, s.tracedCorners) << "%)\n";
            out << " Valid residual corners:      " << s.validResidualCorners
                << " (" << percentForPaper(s.validResidualCorners, s.tracedCorners) << "%)\n";
            out << " Invalid residual corners:    " << s.invalidResidualCorners
                << " (" << percentForPaper(s.invalidResidualCorners, s.tracedCorners) << "%)\n";

            out << " Scale < 1:                   " << s.scaleBelowOne
                << " (" << percentForPaper(s.scaleBelowOne, s.validResidualCorners) << "%)\n";
            out << " Scale = 1:                   " << s.scaleEqualOne
                << " (" << percentForPaper(s.scaleEqualOne, s.validResidualCorners) << "%)\n";
            out << " Scale > 1:                   " << s.scaleAboveOne
                << " (" << percentForPaper(s.scaleAboveOne, s.validResidualCorners) << "%)\n";

            out << " Scale near min limit:          " << s.scaleAtMinClamp
                << " (" << percentForPaper(s.scaleAtMinClamp, s.validResidualCorners) << "%)\n";
            out << " Scale at max clamp:          " << s.scaleAtMaxClamp
                << " (" << percentForPaper(s.scaleAtMaxClamp, s.validResidualCorners) << "%)\n";

            printScalar("Observed residual R_abs", s.observedResidual);
            printScalar("Predicted parent error E_pred", s.predictedParentError);
            printScalar("Residual ratio rho", s.residualRatio);
            printScalar("Correction scale s", s.correctionScale);

            out << " Prune-side signal rho < rhoIgnoreMin: "
                << s.overEstimateSignal
                << " (" << percentForPaper(s.overEstimateSignal, s.validResidualCorners) << "%)\n";

            out << " Refine-side signal rho > rhoIgnoreMax: "
                << s.underEstimateSignal
                << " (" << percentForPaper(s.underEstimateSignal, s.validResidualCorners) << "%)\n";

            out << " Ignored dead-zone samples: "
                << s.nearConsistentSignal
                << " (" << percentForPaper(s.nearConsistentSignal, s.validResidualCorners) << "%)\n";

            out << " Zero/ignored correction samples: "
                << s.lowConfidenceSamples
                << " (" << percentForPaper(s.lowConfidenceSamples, s.validResidualCorners) << "%)\n";

            out << " Base-or-stronger correction samples: "
                << s.highConfidenceSamples
                << " (" << percentForPaper(s.highConfidenceSamples, s.validResidualCorners) << "%)\n";

            printScalar("Abs log ratio", s.logRatioAbs);
            printScalar("Correction confidence/gain", s.disagreementConfidence);
            printScalar("Signal confidence", s.signalConfidence);
            printScalar("Final correction gain", s.residualConfidence);
        };

    auto printDecision = [&](const char* title, const ResidualDecisionPaperStats& d)
        {
            out << "\n";
            out << title << "\n";
            out << "--------------------------------------------------------------------------------\n";

            uint64_t originalSubdivide =
                d.bothSubdivide + d.prunedByResidual;

            uint64_t correctedSubdivide =
                d.bothSubdivide + d.addedByResidual;

            out << " Tested cells:                " << d.testedCells << "\n";

            out << " Original-Hessian subdivide:  " << originalSubdivide
                << " (" << percentForPaper(originalSubdivide, d.testedCells) << "%)\n";

            out << " Corrected subdivide:         " << correctedSubdivide
                << " (" << percentForPaper(correctedSubdivide, d.testedCells) << "%)\n";

            out << " Both subdivide:              " << d.bothSubdivide
                << " (" << percentForPaper(d.bothSubdivide, d.testedCells) << "%)\n";

            out << " Both leaf:                   " << d.bothLeaf
                << " (" << percentForPaper(d.bothLeaf, d.testedCells) << "%)\n";

            out << " Pruned by residual:          " << d.prunedByResidual
                << " (" << percentForPaper(d.prunedByResidual, d.testedCells) << "%)\n";

            out << " Added by residual:           " << d.addedByResidual
                << " (" << percentForPaper(d.addedByResidual, d.testedCells) << "%)\n";

            out << " Pruning rate among original subdivisions: "
                << percentForPaper(d.prunedByResidual, originalSubdivide) << "%\n";

            out << " Added-refinement rate among original leaves: "
                << percentForPaper(d.addedByResidual, d.bothLeaf + d.addedByResidual) << "%\n";

            printScalar("Original Hessian error", d.originalHessianError);
            printScalar("Residual-corrected error", d.correctedResidualError);
            printScalar("Corrected / original ratio", d.correctedOverOriginalRatio);
        };

    out << "\n";
    out << "================================================================================\n";
    out << " RESIDUAL PAPER STATISTICS\n";
    out << "================================================================================\n";

    out << " Threshold:                   " << mCurrentThreshold << "\n";
    out << " Residual correction enabled: " << (mUseResidualCorrection ? "Yes" : "No") << "\n";
    //out << " Eta:                         " << mResidualCorrectionEta << "\n";
    out << " Scale min / max:             "
        << mResidualCorrectionMinScale << " / "
        << mResidualCorrectionMaxScale << "\n";

    printSample("ALL RESIDUAL SAMPLES", mResidualPaperStats.allSamples);
    printDecision("ALL CELL DECISIONS", mResidualPaperStats.allDecisions);

    out << "\n";
    out << "================================================================================\n";
    out << " RESIDUAL SAMPLES BY PARENT LEVEL\n";
    out << "================================================================================\n";

    for (int level = 0; level < kPaperStatsMaxLevel; ++level)
    {
        const auto& s = mResidualPaperStats.samplesByParentLevel[level];

        if (s.tracedCorners == 0)
            continue;

        std::string title = "PARENT LEVEL " + std::to_string(level);
        printSample(title.c_str(), s);
    }

    out << "\n";
    out << "================================================================================\n";
    out << " CELL DECISIONS BY CELL LEVEL\n";
    out << "================================================================================\n";

    for (int level = 0; level < kPaperStatsMaxLevel; ++level)
    {
        const auto& d = mResidualPaperStats.decisionsByCellLevel[level];

        if (d.testedCells == 0)
            continue;

        std::string title = "CELL LEVEL " + std::to_string(level);
        printDecision(title.c_str(), d);
    }

    out << "================================================================================\n";
}

void AdaptiveProbeVolume::writeResidualPaperStatsLog(
    const std::string& filename
) const
{
    std::ofstream out(filename);
    if (!out)
        return;

    out << std::fixed << std::setprecision(8);

    printResidualPaperStats(out);
}

void AdaptiveProbeVolume::exportResidualPaperStatsCSV(
    const std::string& summaryCsv,
    const std::string& levelCsv
) const
{
    auto writeScalarColumns = [](std::ofstream& csv, const PaperScalarStats& s)
        {
            csv
                << s.count << ","
                << s.mean() << ","
                << s.stddev() << ","
                << percentileForPaper(s.values, 50.0) << ","
                << percentileForPaper(s.values, 5.0) << ","
                << percentileForPaper(s.values, 95.0) << ","
                << s.minOrZero() << ","
                << s.maxOrZero();
        };

    auto writeSampleColumns = [&](std::ofstream& csv, const ResidualSamplePaperStats& s)
        {
            csv
                << s.tracedCorners << ","
                << s.noParentCorners << ","
                << s.parentCorners << ","
                << s.validResidualCorners << ","
                << s.invalidResidualCorners << ","
                << percentForPaper(s.validResidualCorners, s.tracedCorners) << ","
                << s.scaleBelowOne << ","
                << s.scaleEqualOne << ","
                << s.scaleAboveOne << ","
                << s.scaleAtMinClamp << ","
                << s.scaleAtMaxClamp << ","
                << percentForPaper(s.scaleBelowOne, s.validResidualCorners) << ","
                << percentForPaper(s.scaleAboveOne, s.validResidualCorners) << ","
                << percentForPaper(s.scaleAtMinClamp, s.validResidualCorners) << ","
                << percentForPaper(s.scaleAtMaxClamp, s.validResidualCorners) << ","
                << s.overEstimateSignal << ","
                << s.underEstimateSignal << ","
                << s.nearConsistentSignal << ","
                << s.lowConfidenceSamples << ","
                << s.highConfidenceSamples << ","
                << percentForPaper(s.overEstimateSignal, s.validResidualCorners) << ","
                << percentForPaper(s.underEstimateSignal, s.validResidualCorners) << ","
                << percentForPaper(s.nearConsistentSignal, s.validResidualCorners) << ","
                << percentForPaper(s.lowConfidenceSamples, s.validResidualCorners) << ","
                << percentForPaper(s.highConfidenceSamples, s.validResidualCorners) << ",";
            writeScalarColumns(csv, s.observedResidual);
            csv << ",";
            writeScalarColumns(csv, s.predictedParentError);
            csv << ",";
            writeScalarColumns(csv, s.residualRatio);
            csv << ",";
            writeScalarColumns(csv, s.correctionScale);
        };

    auto writeDecisionColumns = [&](std::ofstream& csv, const ResidualDecisionPaperStats& d)
        {
            uint64_t originalSubdivide =
                d.bothSubdivide + d.prunedByResidual;

            uint64_t correctedSubdivide =
                d.bothSubdivide + d.addedByResidual;

            csv
                << d.testedCells << ","
                << originalSubdivide << ","
                << correctedSubdivide << ","
                << d.bothSubdivide << ","
                << d.bothLeaf << ","
                << d.prunedByResidual << ","
                << d.addedByResidual << ","
                << percentForPaper(originalSubdivide, d.testedCells) << ","
                << percentForPaper(correctedSubdivide, d.testedCells) << ","
                << percentForPaper(d.prunedByResidual, originalSubdivide) << ","
                << percentForPaper(d.addedByResidual, d.bothLeaf + d.addedByResidual) << ",";

            writeScalarColumns(csv, d.originalHessianError);
            csv << ",";
            writeScalarColumns(csv, d.correctedResidualError);
            csv << ",";
            writeScalarColumns(csv, d.correctedOverOriginalRatio);
        };

    MemoryFootprintInfo mem = calculateMemoryFootprint();

    uint64_t totalProbes = uint64_t(mProbes.size());
    uint64_t leafProbes = 0;

    for (const Probe& p : mProbes)
    {
        if (p.isLeaf)
            leafProbes++;
    }

    // ------------------------------------------------------------
    // Summary CSV: one row per run.
    // ------------------------------------------------------------
    {
        std::ofstream csv(summaryCsv);
        if (csv)
        {
            csv << std::fixed << std::setprecision(8);

            csv
                << "threshold,"
                << "residualEnabled,"
                << "pruneStrength,"
                << "refineStrength,"
                << "rhoIgnoreMin,"
                << "rhoIgnoreMax,"
                << "signalConfidenceEps,"
                << "scaleMin,"
                << "scaleMax,"
                << "maxLevel,"
                << "totalProbes,"
                << "leafProbes,"
                << "corners,"
                << "gpuCornerBytes,"
                << "gpuProbeBytes,"
                << "totalBytes,"
                << "totalMB,"
                << "sample_tracedCorners,"
                << "sample_noParentCorners,"
                << "sample_parentCorners,"
                << "sample_validResidualCorners,"
                << "sample_invalidResidualCorners,"
                << "sample_validResidualPercent,"
                << "sample_scaleBelowOne,"
                << "sample_scaleEqualOne,"
                << "sample_scaleAboveOne,"
                << "sample_scaleAtMinClamp,"
                << "sample_scaleAtMaxClamp,"
                << "sample_scaleBelowOnePercent,"
                << "sample_scaleAboveOnePercent,"
                << "sample_scaleAtMinClampPercent,"
                << "sample_scaleAtMaxClampPercent,"
                << "sample_pruneSideSignal,"
                << "sample_refineSideSignal,"
                << "sample_ignoredDeadZone,"
                << "sample_zeroIgnoredCorrection,"
                << "sample_baseOrStrongerCorrection,"
                << "sample_pruneSideSignalPercent,"
                << "sample_refineSideSignalPercent,"
                << "sample_ignoredDeadZonePercent,"
                << "sample_zeroIgnoredCorrectionPercent,"
                << "sample_baseOrStrongerCorrectionPercent,"
                << "observed_count,observed_mean,observed_std,observed_median,observed_p05,observed_p95,observed_min,observed_max,"
                << "predicted_count,predicted_mean,predicted_std,predicted_median,predicted_p05,predicted_p95,predicted_min,predicted_max,"
                << "rho_count,rho_mean,rho_std,rho_median,rho_p05,rho_p95,rho_min,rho_max,"
                << "scale_count,scale_mean,scale_std,scale_median,scale_p05,scale_p95,scale_minValue,scale_maxValue,"
                << "decision_testedCells,"
                << "decision_originalSubdivide,"
                << "decision_correctedSubdivide,"
                << "decision_bothSubdivide,"
                << "decision_bothLeaf,"
                << "decision_prunedByResidual,"
                << "decision_addedByResidual,"
                << "decision_originalSubdividePercent,"
                << "decision_correctedSubdividePercent,"
                << "decision_pruningRateAmongOriginalSubdivisions,"
                << "decision_addedRateAmongOriginalLeaves,"
                << "origErr_count,origErr_mean,origErr_std,origErr_median,origErr_p05,origErr_p95,origErr_min,origErr_max,"
                << "corrErr_count,corrErr_mean,corrErr_std,corrErr_median,corrErr_p05,corrErr_p95,corrErr_min,corrErr_max,"
                << "corrRatio_count,corrRatio_mean,corrRatio_std,corrRatio_median,corrRatio_p05,corrRatio_p95,corrRatio_min,corrRatio_max"
                << "\n";

            csv
                << mCurrentThreshold << ","
                << (mUseResidualCorrection ? 1 : 0) << ","
                << mResidualPruneStrength << ","
                << mResidualRefineStrength << ","
                << mResidualConfidenceTau0 << ","
                << mResidualConfidenceTau1 << ","
                << mResidualConfidenceEps << ","
                << mResidualCorrectionMinScale << ","
                << mResidualCorrectionMaxScale << ","
                << mMaxLevel << ","
                << totalProbes << ","
                << leafProbes << ","
                << mCorners.size() << ","
                << mem.gpuCornersBytes << ","
                << mem.gpuProbesBytes << ","
                << mem.totalBytes << ","
                << mem.totalMB << ",";

            writeSampleColumns(csv, mResidualPaperStats.allSamples);
            csv << ",";
            writeDecisionColumns(csv, mResidualPaperStats.allDecisions);
            csv << "\n";
        }
    }

    // ------------------------------------------------------------
    // Level CSV: one row per level.
    // Useful for plots: pruning rate by level, scale distribution by level.
    // ------------------------------------------------------------
    {
        std::ofstream csv(levelCsv);
        if (csv)
        {
            csv << std::fixed << std::setprecision(8);

            csv
                << "kind,"
                << "level,"
                << "tracedCorners,"
                << "noParentCorners,"
                << "parentCorners,"
                << "validResidualCorners,"
                << "invalidResidualCorners,"
                << "validResidualPercent,"
                << "scaleBelowOne,"
                << "scaleEqualOne,"
                << "scaleAboveOne,"
                << "scaleAtMinClamp,"
                << "scaleAtMaxClamp,"
                << "scaleBelowOnePercent,"
                << "scaleAboveOnePercent,"
                << "scaleAtMinClampPercent,"
                << "scaleAtMaxClampPercent,"
                << "sample_pruneSideSignal,"
                << "sample_refineSideSignal,"
                << "sample_ignoredDeadZone,"
                << "sample_zeroIgnoredCorrection,"
                << "sample_baseOrStrongerCorrection,"
                << "sample_pruneSideSignalPercent,"
                << "sample_refineSideSignalPercent,"
                << "sample_ignoredDeadZonePercent,"
                << "sample_zeroIgnoredCorrectionPercent,"
                << "sample_baseOrStrongerCorrectionPercent,"
                << "observed_count,observed_mean,observed_std,observed_median,observed_p05,observed_p95,observed_min,observed_max,"
                << "predicted_count,predicted_mean,predicted_std,predicted_median,predicted_p05,predicted_p95,predicted_min,predicted_max,"
                << "rho_count,rho_mean,rho_std,rho_median,rho_p05,rho_p95,rho_min,rho_max,"
                << "scale_count,scale_mean,scale_std,scale_median,scale_p05,scale_p95,scale_minValue,scale_maxValue,"
                << "testedCells,"
                << "originalSubdivide,"
                << "correctedSubdivide,"
                << "bothSubdivide,"
                << "bothLeaf,"
                << "prunedByResidual,"
                << "addedByResidual,"
                << "originalSubdividePercent,"
                << "correctedSubdividePercent,"
                << "pruningRateAmongOriginalSubdivisions,"
                << "addedRateAmongOriginalLeaves,"
                << "origErr_count,origErr_mean,origErr_std,origErr_median,origErr_p05,origErr_p95,origErr_min,origErr_max,"
                << "corrErr_count,corrErr_mean,corrErr_std,corrErr_median,corrErr_p05,corrErr_p95,corrErr_min,corrErr_max,"
                << "corrRatio_count,corrRatio_mean,corrRatio_std,corrRatio_median,corrRatio_p05,corrRatio_p95,corrRatio_min,corrRatio_max"
                << "\n";

            for (int level = 0; level < kPaperStatsMaxLevel; ++level)
            {
                const auto& s = mResidualPaperStats.samplesByParentLevel[level];
                const auto& d = mResidualPaperStats.decisionsByCellLevel[level];

                if (s.tracedCorners == 0 && d.testedCells == 0)
                    continue;

                csv << "level," << level << ",";

                writeSampleColumns(csv, s);
                csv << ",";
                writeDecisionColumns(csv, d);
                csv << "\n";
            }
        }
    }
}

void AdaptiveProbeVolume::printResidualMainStats(std::ostream& out) const
{
    const auto& s = mResidualPaperStats.allSamples;
    const auto& d = mResidualPaperStats.allDecisions;

    auto percent = [](uint64_t n, uint64_t total) -> double
        {
            return total > 0
                ? 100.0 * double(n) / double(total)
                : 0.0;
        };

    uint64_t leafCount = 0;
    for (const Probe& p : mProbes)
    {
        if (p.isLeaf)
            leafCount++;
    }

    MemoryFootprintInfo mem = calculateMemoryFootprint();

    uint64_t originalSubdivide =
        d.bothSubdivide + d.prunedByResidual;

    uint64_t correctedSubdivide =
        d.bothSubdivide + d.addedByResidual;

    uint64_t originalLeaf =
        d.bothLeaf + d.addedByResidual;

    double validResidualPercent =
        percent(s.validResidualCorners, s.tracedCorners);

    double scaleBelowOnePercent =
        percent(s.scaleBelowOne, s.validResidualCorners);

    double scaleAboveOnePercent =
        percent(s.scaleAboveOne, s.validResidualCorners);

    double minClampPercent =
        percent(s.scaleAtMinClamp, s.validResidualCorners);

    double maxClampPercent =
        percent(s.scaleAtMaxClamp, s.validResidualCorners);

    double pruningRateAmongOriginalSubdivisions =
        percent(d.prunedByResidual, originalSubdivide);

    double addedRateAmongOriginalLeaves =
        percent(d.addedByResidual, originalLeaf);

    double originalSubdividePercent =
        percent(originalSubdivide, d.testedCells);

    double correctedSubdividePercent =
        percent(correctedSubdivide, d.testedCells);

    out << std::fixed << std::setprecision(8);

    out << "================================================================================\n";
    out << " RESIDUAL MAIN STATISTICS\n";
    out << "================================================================================\n";

    out << " Threshold:                         " << mCurrentThreshold << "\n";
    out << " Residual correction enabled:       " << (mUseResidualCorrection ? "Yes" : "No") << "\n";
    out << " Prune strength:                    " << mResidualPruneStrength << "\n";
    out << " Refine strength:                   " << mResidualRefineStrength << "\n";
    out << " Rho ignore range min / max:        "
        << std::min(mResidualConfidenceTau0, mResidualConfidenceTau1)
        << " / "
        << std::max(mResidualConfidenceTau0, mResidualConfidenceTau1)
        << "\n";

    out << " Signal confidence epsilon:         "
        << mResidualConfidenceEps << "\n";
    out << " Scale min / max:                   "
        << mResidualCorrectionMinScale << " / "
        << mResidualCorrectionMaxScale << "\n";

    out << "\n";

    out << " Total probes:                      " << mProbes.size() << "\n";
    out << " Leaf probes:                       " << leafCount << "\n";
    out << " Corners:                           " << mCorners.size() << "\n";
    out << " Total memory MB:                   " << mem.totalMB << "\n";
    out << " Build time ms:                     " << mBuildTimeMs << "\n";

    out << "\n";

    out << " Valid residual samples:            "
        << s.validResidualCorners
        << " / " << s.tracedCorners
        << " (" << validResidualPercent << "%)\n";

    out << " Scale < 1:                         "
        << s.scaleBelowOne
        << " (" << scaleBelowOnePercent << "%)\n";

    out << " Scale > 1:                         "
        << s.scaleAboveOne
        << " (" << scaleAboveOnePercent << "%)\n";

    out << " Scale near min limit:                "
        << s.scaleAtMinClamp
        << " (" << minClampPercent << "%)\n";

    out << " Scale near max limit:                "
        << s.scaleAtMaxClamp
        << " (" << maxClampPercent << "%)\n";

    out << " Mean rho:                          "
        << s.residualRatio.mean() << "\n";

    out << " Median rho:                        "
        << percentileForPaper(s.residualRatio.values, 50.0) << "\n";

    out << " P95 rho:                           "
        << percentileForPaper(s.residualRatio.values, 95.0) << "\n";

    out << " Mean scale:                        "
        << s.correctionScale.mean() << "\n";

    out << " Median scale:                      "
        << percentileForPaper(s.correctionScale.values, 50.0) << "\n";

    out << "\n";

    out << " Tested cells:                      " << d.testedCells << "\n";

    out << " Original-Hessian subdivide:        "
        << originalSubdivide
        << " (" << originalSubdividePercent << "%)\n";

    out << " Residual-corrected subdivide:      "
        << correctedSubdivide
        << " (" << correctedSubdividePercent << "%)\n";

    out << " Pruned by residual:                "
        << d.prunedByResidual
        << " (" << pruningRateAmongOriginalSubdivisions
        << "% of original subdivisions)\n";

    out << " Added by residual:                 "
        << d.addedByResidual
        << " (" << addedRateAmongOriginalLeaves
        << "% of original leaves)\n";

    out << " Mean corrected/original error:     "
        << d.correctedOverOriginalRatio.mean() << "\n";

    out << " Median corrected/original error:   "
        << percentileForPaper(d.correctedOverOriginalRatio.values, 50.0) << "\n";

    out << " Prune-side signal rho < rhoIgnoreMin: "
        << s.overEstimateSignal
        << " (" << percent(s.overEstimateSignal, s.validResidualCorners) << "%)\n";

    out << " Refine-side signal rho > rhoIgnoreMax: "
        << s.underEstimateSignal
        << " (" << percent(s.underEstimateSignal, s.validResidualCorners) << "%)\n";

    out << " Ignored dead-zone samples:          "
        << s.nearConsistentSignal
        << " (" << percent(s.nearConsistentSignal, s.validResidualCorners) << "%)\n";

    out << " Zero/ignored correction samples:    "
        << s.lowConfidenceSamples
        << " (" << percent(s.lowConfidenceSamples, s.validResidualCorners) << "%)\n";

    out << " Base-or-stronger correction samples:"
        << s.highConfidenceSamples
        << " (" << percent(s.highConfidenceSamples, s.validResidualCorners) << "%)\n";

    out << " Mean correction gain:               "
        << s.residualConfidence.mean() << "\n";

    out << " Median correction gain:             "
        << percentileForPaper(s.residualConfidence.values, 50.0) << "\n";
    out << "================================================================================\n";
}

void AdaptiveProbeVolume::writeResidualMainStatsLog(
    const std::string& filename
) const
{
    std::ofstream out(filename);
    if (!out)
        return;

    printResidualMainStats(out);
}

void AdaptiveProbeVolume::exportResidualMainStatsCSV(
    const std::string& filename
) const
{
    std::ofstream csv(filename);
    if (!csv)
        return;

    const auto& s = mResidualPaperStats.allSamples;
    const auto& d = mResidualPaperStats.allDecisions;

    auto percent = [](uint64_t n, uint64_t total) -> double
        {
            return total > 0
                ? 100.0 * double(n) / double(total)
                : 0.0;
        };

    uint64_t leafCount = 0;
    for (const Probe& p : mProbes)
    {
        if (p.isLeaf)
            leafCount++;
    }

    MemoryFootprintInfo mem = calculateMemoryFootprint();

    uint64_t originalSubdivide =
        d.bothSubdivide + d.prunedByResidual;

    uint64_t correctedSubdivide =
        d.bothSubdivide + d.addedByResidual;

    uint64_t originalLeaf =
        d.bothLeaf + d.addedByResidual;

    csv << std::fixed << std::setprecision(8);

    csv
        << "threshold,"
        << "residualEnabled,"
        << "pruneStrength,"
        << "refineStrength,"
        << "rhoIgnoreMin,"
        << "rhoIgnoreMax,"
        << "signalConfidenceEps,"
        << "scaleMin,"
        << "scaleMax,"
        << "totalProbes,"
        << "leafProbes,"
        << "corners,"
        << "memoryMB,"
        << "buildTimeMs,"
        << "validResidualSamples,"
        << "totalTracedCorners,"
        << "validResidualPercent,"
        << "scaleBelowOnePercent,"
        << "scaleAboveOnePercent,"
        << "scaleAtMinClampPercent,"
        << "scaleAtMaxClampPercent,"
        << "rhoMean,"
        << "rhoMedian,"
        << "rhoP95,"
        << "scaleMean,"
        << "scaleMedian,"
        << "testedCells,"
        << "originalSubdivide,"
        << "correctedSubdivide,"
        << "originalSubdividePercent,"
        << "correctedSubdividePercent,"
        << "prunedByResidual,"
        << "addedByResidual,"
        << "pruningRateAmongOriginalSubdivisions,"
        << "addedRateAmongOriginalLeaves,"
        << "correctedOverOriginalMean,"
        << "correctedOverOriginalMedian"
        << "\n";

    csv
        << mCurrentThreshold << ","
        << (mUseResidualCorrection ? 1 : 0) << ","
        << mResidualPruneStrength << ","
        << mResidualRefineStrength << ","
        << std::min(mResidualConfidenceTau0, mResidualConfidenceTau1) << ","
        << std::max(mResidualConfidenceTau0, mResidualConfidenceTau1) << ","
        << mResidualConfidenceEps << ","
        << mResidualCorrectionMinScale << ","
        << mResidualCorrectionMaxScale << ","
        << mProbes.size() << ","
        << leafCount << ","
        << mCorners.size() << ","
        << mem.totalMB << ","
        << mBuildTimeMs << ","
        << s.validResidualCorners << ","
        << s.tracedCorners << ","
        << percent(s.validResidualCorners, s.tracedCorners) << ","
        << percent(s.scaleBelowOne, s.validResidualCorners) << ","
        << percent(s.scaleAboveOne, s.validResidualCorners) << ","
        << percent(s.scaleAtMinClamp, s.validResidualCorners) << ","
        << percent(s.scaleAtMaxClamp, s.validResidualCorners) << ","
        << s.residualRatio.mean() << ","
        << percentileForPaper(s.residualRatio.values, 50.0) << ","
        << percentileForPaper(s.residualRatio.values, 95.0) << ","
        << s.correctionScale.mean() << ","
        << percentileForPaper(s.correctionScale.values, 50.0) << ","
        << d.testedCells << ","
        << originalSubdivide << ","
        << correctedSubdivide << ","
        << percent(originalSubdivide, d.testedCells) << ","
        << percent(correctedSubdivide, d.testedCells) << ","
        << d.prunedByResidual << ","
        << d.addedByResidual << ","
        << percent(d.prunedByResidual, originalSubdivide) << ","
        << percent(d.addedByResidual, originalLeaf) << ","
        << d.correctedOverOriginalRatio.mean() << ","
        << percentileForPaper(d.correctedOverOriginalRatio.values, 50.0)
        << "\n";
}

static float quadraticForm3x3(
    const float3x3& H,
    const float3& d
)
{
    float3 Hd;

    Hd.x =
        H[0][0] * d.x +
        H[0][1] * d.y +
        H[0][2] * d.z;

    Hd.y =
        H[1][0] * d.x +
        H[1][1] * d.y +
        H[1][2] * d.z;

    Hd.z =
        H[2][0] * d.x +
        H[2][1] * d.y +
        H[2][2] * d.z;

    return dot(d, Hd);
}

static float3 worldDisplacementToHessianSpace(
    const float3& dWorld
)
{
    // Hessians are computed in your polar/swizzled coordinate order.
    //
    // x_hessian = world Z
    // y_hessian = world X
    // z_hessian = world Y
    return float3(
        dWorld.z,
        dWorld.x,
        dWorld.y
    );
}

static float computeDirectionalHessianResidualPredictionSHSpace(
    const std::vector<float3x3>& hessians,
    const float3& dWorld
)
{
    if (hessians.empty())
        return 0.0f;

    float3 d =
        worldDisplacementToHessianSpace(dWorld);

    const size_t basisCount =
        std::min<size_t>(hessians.size(), 9);

    float sumSq = 0.0f;

    for (size_t basisIdx = 0; basisIdx < basisCount; ++basisIdx)
    {
        // Direction-aware second-order prediction for this SH coefficient:
        //
        //     0.5 * d^T H_lm d
        //
        // We square per coefficient and then take the L2 norm.
        // This avoids cancellation across SH coefficients.
        float q =
            0.5f *
            quadraticForm3x3(
                hessians[basisIdx],
                d
            );

        sumSq += q * q;
    }

    return std::sqrt(sumSq);
}

//float AdaptiveProbeVolume::computeParentPointHessianPrediction(
//    int parentProbeIdx,
//    const float3& position
//) const
//{
//    if (parentProbeIdx < 0 || parentProbeIdx >= int(mProbes.size()))
//        return 0.0f;
//
//    const Probe& p = mProbes[parentProbeIdx];
//
//    float3 pMin = mCorners[p.corners[0]].position;
//    float3 pMax = pMin;
//
//    for (int k = 1; k < 8; ++k)
//    {
//        const float3& cPos = mCorners[p.corners[k]].position;
//
//        pMin.x = std::min(pMin.x, cPos.x);
//        pMin.y = std::min(pMin.y, cPos.y);
//        pMin.z = std::min(pMin.z, cPos.z);
//
//        pMax.x = std::max(pMax.x, cPos.x);
//        pMax.y = std::max(pMax.y, cPos.y);
//        pMax.z = std::max(pMax.z, cPos.z);
//    }
//
//    float3 extent = pMax - pMin;
//
//    const float eps = 1e-8f;
//
//    float tx = extent.x > eps
//        ? (position.x - pMin.x) / extent.x
//        : 0.5f;
//
//    float ty = extent.y > eps
//        ? (position.y - pMin.y) / extent.y
//        : 0.5f;
//
//    float tz = extent.z > eps
//        ? (position.z - pMin.z) / extent.z
//        : 0.5f;
//
//    tx = std::max(0.0f, std::min(1.0f, tx));
//    ty = std::max(0.0f, std::min(1.0f, ty));
//    tz = std::max(0.0f, std::min(1.0f, tz));
//
//    float3 center = 0.5f * (pMin + pMax);
//
//    float predictedPointResidual = 0.0f;
//
//    for (int k = 0; k < 8; ++k)
//    {
//        const Corner& c = mCorners[p.corners[k]];
//        const float3& cPos = c.position;
//
//        float wx = cPos.x >= center.x ? tx : (1.0f - tx);
//        float wy = cPos.y >= center.y ? ty : (1.0f - ty);
//        float wz = cPos.z >= center.z ? tz : (1.0f - tz);
//
//        float w = wx * wy * wz;
//
//        float3 dx = position - cPos;
//        float distSq = dot(dx, dx);
//
//        float cornerPred =
//            0.5f * c.maxLambdaVecL2 * distSq;
//
//        predictedPointResidual += w * cornerPred;
//    }
//
//    return predictedPointResidual;
//}

float AdaptiveProbeVolume::computeParentPointHessianPrediction(
    int parentProbeIdx,
    const float3& position
) const
{
    if (parentProbeIdx < 0 ||
        parentProbeIdx >= int(mProbes.size()))
    {
        return 0.0f;
    }

    const Probe& p =
        mProbes[parentProbeIdx];

    float3 pMin =
        mCorners[p.corners[0]].position;

    float3 pMax =
        pMin;

    for (int k = 1; k < 8; ++k)
    {
        const float3& cPos =
            mCorners[p.corners[k]].position;

        pMin.x = std::min(pMin.x, cPos.x);
        pMin.y = std::min(pMin.y, cPos.y);
        pMin.z = std::min(pMin.z, cPos.z);

        pMax.x = std::max(pMax.x, cPos.x);
        pMax.y = std::max(pMax.y, cPos.y);
        pMax.z = std::max(pMax.z, cPos.z);
    }

    float3 extent =
        pMax - pMin;

    const float eps =
        1e-8f;

    float tx =
        extent.x > eps
        ? (position.x - pMin.x) / extent.x
        : 0.5f;

    float ty =
        extent.y > eps
        ? (position.y - pMin.y) / extent.y
        : 0.5f;

    float tz =
        extent.z > eps
        ? (position.z - pMin.z) / extent.z
        : 0.5f;

    tx = std::max(0.0f, std::min(1.0f, tx));
    ty = std::max(0.0f, std::min(1.0f, ty));
    tz = std::max(0.0f, std::min(1.0f, tz));

    float3 center =
        0.5f * (pMin + pMax);

    // Signed predicted SH residual vector.
    //
    // We accumulate signed quadratic predictions per SH coefficient,
    // then take the L2 norm only after all parent-corner contributions
    // have been blended.
    float predictedResidualSH[9] =
    {
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f
    };

    float worstCaseFallback =
        0.0f;

    bool hasDirectionalHessian =
        false;

    for (int k = 0; k < 8; ++k)
    {
        int cornerIdx =
            p.corners[k];

        if (cornerIdx < 0 ||
            cornerIdx >= int(mCorners.size()))
        {
            continue;
        }

        const Corner& c =
            mCorners[cornerIdx];

        const float3& cPos =
            c.position;

        float wx =
            cPos.x >= center.x
            ? tx
            : (1.0f - tx);

        float wy =
            cPos.y >= center.y
            ? ty
            : (1.0f - ty);

        float wz =
            cPos.z >= center.z
            ? tz
            : (1.0f - tz);

        float w =
            wx * wy * wz;

        float3 dWorld =
            position - cPos;

        float distSq =
            dot(dWorld, dWorld);

        // Old conservative scalar prediction.
        // Used only as fallback if Hessians are unavailable.
        worstCaseFallback +=
            w *
            0.5f *
            c.maxLambdaVecL2 *
            distSq;

        if (c.shHessians.empty())
            continue;

        hasDirectionalHessian = true;

        float3 d =
            worldDisplacementToHessianSpace(dWorld);

        const size_t basisCount =
            std::min<size_t>(c.shHessians.size(), 9);

        for (size_t basisIdx = 0; basisIdx < basisCount; ++basisIdx)
        {
            float q =
                0.5f *
                quadraticForm3x3(
                    c.shHessians[basisIdx],
                    d
                );

            // Important:
            // Do NOT take abs(q) here.
            // Do NOT square here.
            // Keep the sign and blend in SH-vector space.
            predictedResidualSH[basisIdx] +=
                w * q;
        }
    }

    if (!hasDirectionalHessian)
        return worstCaseFallback;

    float sumSq =
        0.0f;

    for (int basisIdx = 0; basisIdx < 9; ++basisIdx)
    {
        sumSq +=
            predictedResidualSH[basisIdx] *
            predictedResidualSH[basisIdx];
    }

    float signedVectorPrediction =
        std::sqrt(sumSq);

    if (!std::isfinite(signedVectorPrediction) ||
        signedVectorPrediction <= 0.0f)
    {
        return worstCaseFallback;
    }

    // Test control.
    //
    // 0.0 = old worst-case denominator
    // 1.0 = signed SH-vector directional denominator
    //
    // Try 1.0 first. If it over-refines too much, try 0.5.
    const float signedDirectionalBlend =
        1.0f;

    float prediction =
        (1.0f - signedDirectionalBlend) *
        worstCaseFallback +
        signedDirectionalBlend *
        signedVectorPrediction;

    return prediction;
}

int AdaptiveProbeVolume::findLeafProbeIndexCPU(float3 posW) const
{
    return traverseOctreeCPU(posW);
}
