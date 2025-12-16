#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

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

    mCurrentThreshold = errorThreshold;
    mUseRelativeError = useRelativeError;

    auto bounds = pScene->getSceneBounds();

    // 1. Create Root Probe
    Probe root;
    root.level = 0;
    root.minPoint = bounds.minPoint;
    root.maxPoint = bounds.maxPoint;

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
    const std::vector<HessianSHCoeff>& hessians)
{
    // FIX: batchIndex selects the specific corner in the pending list.
    // The input vectors (coeffs, grads) contain ALL bands for THAT corner.

    if (batchIndex >= mPendingNewCorners.size()) return;

    int cornerIdx = mPendingNewCorners[batchIndex];
    Corner& c = mCorners[cornerIdx];

    // 1. Store Full Data
    c.shCoeffs = coeffs;
    c.shGradients = grads;

    // 2. Compute Norm (for Relative Error)
    // We sum the dot product of all bands
    if (mUseRelativeError) {
        float sumSqL = 0.0f;
        for (const auto& coeff : c.shCoeffs) {
            sumSqL += math::dot(coeff, coeff);
        }
        c.coeffNorm = std::sqrt(sumSqL);
    }
    
    // 3. Compute Curvature (Max Eigenvalue of Hessians)
    // We iterate over the bands (hessians vector) to find the max curvature
    float sumSquares = 0.0f;

    for (const auto& h : hessians)
    {
        float r1, r2, r3, g1, g2, g3, b1, b2, b3;
        computeEigenvalues3x3(h.r, r1, r2, r3);
        computeEigenvalues3x3(h.g, g1, g2, g3);
        computeEigenvalues3x3(h.b, b1, b2, b3);

        float maxR = std::max({ std::abs(r1), std::abs(r2), std::abs(r3) });
        float maxG = std::max({ std::abs(g1), std::abs(g2), std::abs(g3) });
        float maxB = std::max({ std::abs(b1), std::abs(b2), std::abs(b3) });

        float maxLambdaBand = std::max({ maxR, maxG, maxB });
        sumSquares += maxLambdaBand * maxLambdaBand;
    }

    // Total curvature magnitude (L2 norm of the curvature vector across bands)
    c.maxLambda = std::sqrt(sumSquares);
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

        for (int k = 0; k < 8; ++k)
        {
            int cIdx = mProbes[probeIdx].corners[k];
            const Corner& c = mCorners[cIdx];

            // E_abs = 0.5 * Lambda * dx^2
            float e = (c.maxLambda * distSq) * 0.5f;

            if (mUseRelativeError)
            {
                float norm = std::max(c.coeffNorm, 1e-5f);
                totalError += (e / norm);
            }
            else
            {
                totalError += e;
            }
        }
        float avgProbeError = totalError / 8.0f;

        // --------------------------------------------------
        // Step 5: Subdivide
        // --------------------------------------------------
        if (avgProbeError > mCurrentThreshold)
        {
            mProbes[probeIdx].isLeaf = false;

            float3 minP = mProbes[probeIdx].minPoint;
            float3 maxP = mProbes[probeIdx].maxPoint;
            float3 center = (minP + maxP) * 0.5f;

            // Helper: Create new corner
            auto addCorner = [&](float3 pos) -> int {
                Corner c; c.position = pos;
                mCorners.push_back(c);
                int idx = (int)mCorners.size() - 1;
                mPendingNewCorners.push_back(idx);
                return idx;
                };

            // A. Generate 19 NEW Corners
            int c_center = addCorner(center);

            int c_fX0 = addCorner({ minP.x, center.y, center.z }); int c_fX1 = addCorner({ maxP.x, center.y, center.z });
            int c_fY0 = addCorner({ center.x, minP.y, center.z }); int c_fY1 = addCorner({ center.x, maxP.y, center.z });
            int c_fZ0 = addCorner({ center.x, center.y, minP.z }); int c_fZ1 = addCorner({ center.x, center.y, maxP.z });

            int c_eX_Y0Z0 = addCorner({ center.x, minP.y, minP.z }); int c_eX_Y1Z0 = addCorner({ center.x, maxP.y, minP.z });
            int c_eX_Y0Z1 = addCorner({ center.x, minP.y, maxP.z }); int c_eX_Y1Z1 = addCorner({ center.x, maxP.y, maxP.z });

            int c_eY_X0Z0 = addCorner({ minP.x, center.y, minP.z }); int c_eY_X1Z0 = addCorner({ maxP.x, center.y, minP.z });
            int c_eY_X0Z1 = addCorner({ minP.x, center.y, maxP.z }); int c_eY_X1Z1 = addCorner({ maxP.x, center.y, maxP.z });

            int c_eZ_X0Y0 = addCorner({ minP.x, minP.y, center.z }); int c_eZ_X1Y0 = addCorner({ maxP.x, minP.y, center.z });
            int c_eZ_X0Y1 = addCorner({ minP.x, maxP.y, center.z }); int c_eZ_X1Y1 = addCorner({ maxP.x, maxP.y, center.z });

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

void AdaptiveProbeVolume::uploadToGPU()
{
    // Stub
}

void AdaptiveProbeVolume::printDebugInfo(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) return;

    out << "================================================================================\n";
    out << " ADAPTIVE PROBE VOLUME HIERARCHY (Vertex-Sharing)\n";
    out << "================================================================================\n";
    out << " Probes:    " << mProbes.size() << "\n";
    out << " Corners:   " << mCorners.size() << "\n";
    out << " Threshold: " << mCurrentThreshold << "\n";
    out << " Max Level: " << mMaxLevel << "\n";
    out << " Metric:    " << (mUseRelativeError ? "Relative (E_rel)" : "Absolute (E_abs)") << "\n";
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
                maxLambdaInProbe = std::max(maxLambdaInProbe, c.maxLambda);

                // Re-compute error
                float e = (c.maxLambda * distSq) * 0.5f;

                if (mUseRelativeError)
                {
                    float norm = std::max(c.coeffNorm, 1e-5f);
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
