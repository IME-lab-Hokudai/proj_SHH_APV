#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

//

// ------------------------------------------------------------------
// CPU PORT OF SHADER LOGIC
// ------------------------------------------------------------------

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

//

// Interpolates between 'a' and 'b' based on weight 't'.
// If a corner is invalid, it is ignored, and the valid corner takes over 100%.
float3 robustLerp(float3 a, float3 b, float t, bool validA, bool validB)
{
    if (validA && validB)
    {
        // Normal case: Both good, use standard lerp
        return lerp(a, b, t);
    }
    else if (validA && !validB)
    {
        // B is broken (Leak source).
        // Force 'a' to extend all the way to 'b'.
        return a;
    }
    else if (!validA && validB)
    {
        // A is broken.
        // Force 'b' to extend all the way to 'a'.
        return b;
    }
    else
    {
        // Both broken. Return black (or handle as failure upstream).
        return float3(0, 0, 0);
    }
}


//

// CPU Port of the "Smart" Weighted Accumulation Shader
// Matches the shader's logic: Sum(Weight * Value) / Sum(Weight)
void AdaptiveProbeVolume::fetchSmartInterpolatedSH(
    int probeIdx,
    float3 pos,
    float3* outCoeffs
)
{
    const Probe& p = mProbes[probeIdx];
    float3 size = p.maxPoint - p.minPoint;
    float3 t = (pos - p.minPoint) / size;

    // Clamp
    t.x = std::max(0.0f, std::min(1.0f, t.x));
    t.y = std::max(0.0f, std::min(1.0f, t.y));
    t.z = std::max(0.0f, std::min(1.0f, t.z));

    // Pre-calculate Trilinear Weights
    float w[8];
    w[0] = (1 - t.x) * (1 - t.y) * (1 - t.z);
    w[1] = (1 - t.x) * (1 - t.y) * (t.z);
    w[2] = (1 - t.x) * (t.y) * (1 - t.z);
    w[3] = (1 - t.x) * (t.y) * (t.z);
    w[4] = (t.x) * (1 - t.y) * (1 - t.z);
    w[5] = (t.x) * (1 - t.y) * (t.z);
    w[6] = (t.x) * (t.y) * (1 - t.z);
    w[7] = (t.x) * (t.y) * (t.z);

    // Initialize Accumulators
    float3 accumSH[9];
    for (int b = 0; b < 9; ++b) accumSH[b] = float3(0, 0, 0);

    float totalWeight = 0.0f;

    // Loop 8 Corners
    for (int k = 0; k < 8; ++k)
    {
        int cIdx = p.corners[k];

        // 1. DATA CHECK: Skip missing corners
        if (cIdx < 0) continue;

        const Corner& c = mCorners[cIdx];

        // 2. ROBUSTNESS CHECK: Skip invalid/black corners
        // Matches Shader: if (dot(c.coeffs[0], c.coeffs[0]) < 1e-6) continue;
        // We use our helper flag if available, or check the data directly.

        if (!c.isValid ) continue;

        // 3. WEIGHT
        // We cannot check "Visibility" (Normal Dot) here because a hanging node
        // is a point in space, not a surface. We assume "Perfect Visibility" (1.0).
        // This relies on the Spatial Weight (Trilinear) only.
        float weight = w[k];

        // 4. ACCUMULATE
        for (int b = 0; b < 9; ++b)
        {
            accumSH[b] += c.shCoeffs[b] * weight;
        }
        totalWeight += weight;
    }

    // 5. RENORMALIZE
    // This "stretches" the valid corners to fill the gap left by invalid ones.
    if (totalWeight > 1e-6f)
    {
        for (int b = 0; b < 9; ++b)
        {
            outCoeffs[b] = accumSH[b] / totalWeight;
        }
    }
    else
    {
        // Fallback: Total failure. Return Black.
        for (int b = 0; b < 9; ++b) outCoeffs[b] = float3(0, 0, 0);
    }
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
//

//

//

void AdaptiveProbeVolume::constrainHangingNodes()
{
    // Face Indices: 0:-X, 1:+X, 2:-Y, 3:+Y, 4:-Z, 5:+Z
    const int FACE_CORNERS[6][4] = {
             {0, 1, 2, 3}, {4, 5, 6, 7},
             {0, 1, 4, 5}, {2, 3, 6, 7},
             {0, 2, 4, 6}, {1, 3, 5, 7}
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
                for (int k = 0; k < 4; ++k)
                {
                    int localCornerIdx = FACE_CORNERS[face][k];
                    int globalCornerIdx = p.corners[localCornerIdx];
                    Corner& c = mCorners[globalCornerIdx];

                    // 1. Allocate Array on Stack
                    float3 smartCoeffs[9];

                    // 2. Fetch Weighted Average (Renormalized)
                    // This now uses the Accumulation Logic, not the 3-Pass Lerp.
                    fetchSmartInterpolatedSH(neighborIdx, c.position, smartCoeffs);

                    // 3. Overwrite Hanging Node
                    if (c.shCoeffs.size() != 9) c.shCoeffs.resize(9);
                    for (int b = 0; b < 9; ++b) c.shCoeffs[b] = smartCoeffs[b];
                }
            }
        }
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
    //const std::vector<float>& distMeans,
    //const std::vector<float>& distMeanSqs
)
{
    // 1. Locate Probe
    if (batchIndex >= mPendingNewCorners.size()) return;

    int cornerIdx = mPendingNewCorners[batchIndex];
    Corner& c = mCorners[cornerIdx];

    // -----------------------------------------------------------
    // PART A: RUNTIME DATA (Store RGB)
    // -----------------------------------------------------------
    c.shCoeffs = coeffs;       // Full RGB for color rendering
    c.shGradients = grads;     // Full RGB for color interpolation
    //c.distMean = distMeans;
    //c.distMeanSq = distMeanSqs;
    // -----------------------------------------------------------
    // PART B: CONSTRUCTION METRICS (Use Luminance)
    // -----------------------------------------------------------
   // if (mUseRelativeError) {
        const float3 kLuma = float3(0.2126f, 0.7152f, 0.0722f);

        // 1. Calculate Norm of LUMINANCE Coefficients (Vector L)
        // ||L|| = sqrt( sum( (c_i . luma)^2 ) )
        float sumSqL = 0.0f;
        for (const auto& coeff : coeffs)
        {
            // Project this band's RGB coefficient to Luminance
            float lum = dot(coeff, kLuma);
            sumSqL += lum * lum;
        }
        c.coeffVecL2 = std::sqrt(sumSqL);

         c.isValid = (c.coeffVecL2 < 0.0001f)?false : true;
    //}

    // 2. Calculate Curvature of LUMINANCE Field (Hessian)
    float sumSquares = 0.0f;
    for (const auto& h : hessians)
    {
        // Compute Eigenvalues of the Luminance Hessian
        float e1, e2, e3;
        computeEigenvalues3x3(h, e1, e2, e3);

        // Max curvature
        float maxLambda = std::max({ std::abs(e1), std::abs(e2), std::abs(e3) });
        sumSquares += maxLambda * maxLambda;
    }

    c.maxLambdaVecL2 = std::sqrt(sumSquares);
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

        float totalError = 0.0f;
        uint countValidCorner = 0;
        for (int k = 0; k < 8; ++k)
        {
            int cIdx = mProbes[probeIdx].corners[k];
            const Corner& c = mCorners[cIdx];
            if (!c.isValid) continue;
            countValidCorner++;
            // E_abs = 0.5 * Lambda * dx^2
            float e = (c.maxLambdaVecL2 * distSq) * 0.5f;

            if (mUseRelativeError)
            {
                float norm = std::max(c.coeffVecL2, 1e-5f);
                totalError += (e / norm);
            }
            else
            {
                totalError += e;
            }
        }
        float avgProbeError = 0.0f;
        //if (countValidCorner == 0) avgProbeError = mCurrentThreshold + 0.1f;//force subdivision if no 
        avgProbeError = totalError / float(countValidCorner);

        // --------------------------------------------------
        // Step 5: Subdivide
        // --------------------------------------------------
        if (avgProbeError > mCurrentThreshold)
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
// CPU Traversal Logic (Matches Shader Logic)
// ------------------------------------------------------------------
int AdaptiveProbeVolume::traverseOctreeCPU(float3 pos) const
{
    if (mProbes.empty()) return -1;

    // 1. Bounds Check: If position is outside the volume, return -1
    const Probe& root = mProbes[0];
    if (pos.x < root.minPoint.x || pos.y < root.minPoint.y || pos.z < root.minPoint.z ||
        pos.x > root.maxPoint.x || pos.y > root.maxPoint.y || pos.z > root.maxPoint.z)
    {
        return -1;
    }

    int currentIdx = 0; // Start at Root

    // Safety loop to prevent infinite recursion
    for (int i = 0; i < mMaxLevel + 5; ++i)
    {
        const Probe& p = mProbes[currentIdx];

        // Found a Leaf
        if (p.children[0] == -1)
        {
            return currentIdx;
        }

        // Determine Octant
        float3 center = (p.minPoint + p.maxPoint) * 0.5f;
        int octant = 0;
        if (pos.x >= center.x) octant |= 4; // Bit 2 = X
        if (pos.y >= center.y) octant |= 2; // Bit 1 = Y
        if (pos.z >= center.z) octant |= 1; // Bit 0 = Z

        int childIdx = p.children[octant];

        // If tree is incomplete, stop here
        if (childIdx == -1) return currentIdx;

        currentIdx = childIdx;
    }
    return currentIdx;
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
    computeNeighbors();
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
        gc.position = c.position;
        gc.isValid = c.isValid ? 1 : 0;
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

    uint leafCounter = 0;
    for (int i = 0; i < mProbes.size(); ++i) {
        if (mProbes[i].isLeaf)
            leafCounter++;
    }

    out << "================================================================================\n";
    out << " ADAPTIVE PROBE VOLUME HIERARCHY (Vertex-Sharing)\n";
    out << "================================================================================\n";
    out << " Probes:    " << leafCounter << "\n";
    out << " Corners:   " << mCorners.size() << "\n";
    out << " Threshold: " << mCurrentThreshold << "\n";
    out << " Max Level: " << mMaxLevel << "\n";
    out << " Metric:    " << (mUseRelativeError ? "Relative (E_rel)" : "Absolute (E_abs)") << "\n";
    out << " Build time in ms: " << mBuildTimeMs << "\n";
    out << "================================================================================\n\n";

    // Recursive Lambda for Tree Traversal
    std::function<void(int, std::string, bool)> printProbe = [&](int probeIdx, std::string prefix, bool isLast)
        {
            const Probe& probe = mProbes[probeIdx];

            // 1. Re-Calculate Geometry
            float3 diag = probe.maxPoint - probe.minPoint;
            float distSq = dot(diag, diag);
            float size = std::sqrt(distSq);

            // 2. Re-Calculate Error (Must match finishBatch logic exactly)
            float totalError = 0.0f;
            float maxLambdaInProbe = 0.0f;

            for (int k = 0; k < 8; ++k)
            {
                int cIdx = probe.corners[k];
                const Corner& c = mCorners[cIdx];

                // Track max curvature in this voxel for display
                maxLambdaInProbe = std::max(maxLambdaInProbe, c.maxLambdaVecL2);

                // Re-compute error
                float e = (c.maxLambdaVecL2 * distSq) * 0.5f;

                if (mUseRelativeError)
                {
                    float norm = std::max(c.coeffVecL2, 1e-5f);
                    totalError += (e / norm);
                }
                else
                {
                    totalError += e;
                }
            }
            float avgError = totalError / 8.0f;

            // 3. Format Output
            out << prefix << (isLast ? "L-- " : "|-- ");
            out << "[Lvl " << probe.level << "] ";

            // Format error: "0.0050"
            std::stringstream ssErr;
            ssErr << std::fixed << std::setprecision(4) << avgError;
            std::string errStr = ssErr.str();

            if (!probe.isLeaf)
            {
                // Case 1: SUBDIVIDED
                out << "SUBDIVIDED ";
                out << "(Err: " << errStr << " > " << mCurrentThreshold << ") ";
                out << "Size: " << std::setprecision(2) << size << " ";
                out << "MaxLambda: " << std::setprecision(3) << maxLambdaInProbe << "\n";

                // Recurse children
                std::vector<int> childrenIdx;
                for (int i = 0; i < 8; ++i)
                    if (probe.children[i] != -1)
                        childrenIdx.push_back(probe.children[i]);

                for (size_t i = 0; i < childrenIdx.size(); ++i)
                {
                    std::string newPrefix = prefix + (isLast ? "    " : "|   ");
                    printProbe(childrenIdx[i], newPrefix, (i == childrenIdx.size() - 1));
                }
            }
            else
            {
                // Case 2: LEAF
                out << "LEAF       ";

                if (probe.level == mMaxLevel)
                {
                    out << "(Max Depth) ";
                }
                else if (avgError > mCurrentThreshold)
                {
                    // This happens if finishBatch hasn't processed this leaf yet, or logic bug
                    out << "PENDING? (Err: " << errStr << " > Thr) ";
                }
                else
                {
                    out << "(Err: " << errStr << " < " << mCurrentThreshold << ") ";
                }
                out << "MaxLambda: " << std::setprecision(3) << maxLambdaInProbe << "\n";
            }
        };

    if (!mProbes.empty())
    {
        printProbe(0, "", true);
    }
    else
    {
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

    // 1. Header & Config (UPDATED VERSION TO V3)
    out << "ADAPTIVE_GRID_V4\n"; //add build time
    out << mCurrentThreshold << " " << mMaxLevel << " " << (mUseRelativeError ? 1 : 0) << "\n";
    out << mBuildTimeMs << "\n";
    // 2. Corners
    out << "NUM_CORNERS " << mCorners.size() << "\n";
    for (const auto& c : mCorners)
    {
        // Basic Data
        out << c.position.x << " " << c.position.y << " " << c.position.z << " ";
        out << c.maxLambdaVecL2 << " " << c.coeffVecL2 << "\n";

        // SH Coeffs (Radiance)
        out << c.shCoeffs.size() << "\n";
        for (const auto& val : c.shCoeffs)
        {
            out << val.x << " " << val.y << " " << val.z << " ";
        }
        out << "\n";

        // SH Gradients
        out << c.shGradients.size() << "\n";
        for (const auto& g : c.shGradients)
        {
            out << g.r.x << " " << g.r.y << " " << g.r.z << " ";
            out << g.g.x << " " << g.g.y << " " << g.g.z << " ";
            out << g.b.x << " " << g.b.y << " " << g.b.z << " ";
        }
        out << "\n";

        // --------------------------------------------------------
        // NEW: Distance Moments (Chebyshev Data)
        // --------------------------------------------------------

        // Mean Distance
        //out << c.distMean.size() << "\n";
        //for (const auto& val : c.distMean) out << val << " ";
        //out << "\n";

        //// Mean Squared Distance
        //out << c.distMeanSq.size() << "\n";
        //for (const auto& val : c.distMeanSq) out << val << " ";
        //out << "\n";
    }

    out << "NUM_PROBES " << mProbes.size() << "\n";
    for (const auto& p : mProbes)
    {
        out << (p.isLeaf ? 1 : 0) << " " << p.level << "\n";

        // Bounds
        out << p.minPoint.x << " " << p.minPoint.y << " " << p.minPoint.z << " ";
        out << p.maxPoint.x << " " << p.maxPoint.y << " " << p.maxPoint.z << "\n";

        // Corners
        for (int i = 0; i < 8; ++i) out << p.corners[i] << " ";
        out << "\n";

        // Children
        for (int i = 0; i < 8; ++i) out << p.children[i] << " ";
        out << "\n";
    }

    out.close();
    logInfo("Successfully saved AdaptiveProbeVolume (V3) to " + filename);
}

void AdaptiveProbeVolume::loadFromFile(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
    {
        logError("Failed to open file for loading: " + filename);
        return;
    }

    // Clear existing
    mProbes.clear();
    mCorners.clear();
    mPendingNewCorners.clear();
    mProbesPendingCheck.clear();
    mCornerLookup.clear(); // Cleared here... but needs to be filled!

    std::string header;
    in >> header;

    if (header != "ADAPTIVE_GRID_V4")
    {
        logError("Invalid file format or version mismatch (Expected V3): " + filename);
        return;
    }

    int useRelErrInt;
    in >> mCurrentThreshold >> mMaxLevel >> useRelErrInt;
    mUseRelativeError = (useRelErrInt != 0);
    in >> mBuildTimeMs;
    // 2. Load Corners
    size_t numCorners;
    std::string tag;
    in >> tag >> numCorners;

    mCorners.resize(numCorners);

    for (size_t i = 0; i < numCorners; ++i)
    {
        Corner& c = mCorners[i];
        in >> c.position.x >> c.position.y >> c.position.z;
        in >> c.maxLambdaVecL2 >> c.coeffVecL2;
        c.isValid = c.coeffVecL2 < 0.001 ? false : true;

        // SH Coeffs
        size_t numCoeffs;
        in >> numCoeffs;
        c.shCoeffs.resize(numCoeffs);
        for (size_t k = 0; k < numCoeffs; ++k)
        {
            in >> c.shCoeffs[k].x >> c.shCoeffs[k].y >> c.shCoeffs[k].z;
        }

        // SH Gradients
        size_t numGrads;
        in >> numGrads;
        c.shGradients.resize(numGrads);
        for (size_t k = 0; k < numGrads; ++k)
        {
            in >> c.shGradients[k].r.x >> c.shGradients[k].r.y >> c.shGradients[k].r.z;
            in >> c.shGradients[k].g.x >> c.shGradients[k].g.y >> c.shGradients[k].g.z;
            in >> c.shGradients[k].b.x >> c.shGradients[k].b.y >> c.shGradients[k].b.z;
        }

        // NEW: Distance Moments
        /*size_t numDist;
        in >> numDist;
        c.distMean.resize(numDist);
        for (size_t k = 0; k < numDist; ++k) in >> c.distMean[k];

        size_t numDistSq;
        in >> numDistSq;
        c.distMeanSq.resize(numDistSq);
        for (size_t k = 0; k < numDistSq; ++k) in >> c.distMeanSq[k];*/
    }

    // 3. Load Probes
    size_t numProbes;
    in >> tag >> numProbes;

    mProbes.resize(numProbes);
    for (size_t i = 0; i < numProbes; ++i)
    {
        Probe& p = mProbes[i];
        int isLeafInt;
        in >> isLeafInt >> p.level;
        p.isLeaf = (isLeafInt != 0);

        in >> p.minPoint.x >> p.minPoint.y >> p.minPoint.z;
        in >> p.maxPoint.x >> p.maxPoint.y >> p.maxPoint.z;

        for (int k = 0; k < 8; ++k) in >> p.corners[k];
        for (int k = 0; k < 8; ++k) in >> p.children[k];
    }
    in.close();
    logInfo("Successfully loaded AdaptiveProbeVolume (V3) from " + filename);
}
