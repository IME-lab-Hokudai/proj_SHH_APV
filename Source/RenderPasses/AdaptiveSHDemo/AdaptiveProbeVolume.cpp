#include "AdaptiveProbeVolume.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>


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
        //outCoeffs[band] = float3(
        //    hermite1D(t.x, r[0].x, r_dX[0].x * size.x, r[1].x, r_dX[1].x * size.x),
        //    hermite1D(t.x, r[0].y, r_dX[0].y * size.y, r[1].y, r_dX[1].y * size.x),
        //    hermite1D(t.x, r[0].z, r_dX[0].z * size.z, r[1].z, r_dX[1].z * size.x)
        //);
        outCoeffs[band] = float3(
            hermite1D(t.x, r[0].x, r_dX[0].x * size.x, r[1].x, r_dX[1].x * size.x),
            hermite1D(t.x, r[0].y, r_dX[0].y * size.x, r[1].y, r_dX[1].y * size.x),
            hermite1D(t.x, r[0].z, r_dX[0].z * size.x, r[1].z, r_dX[1].z * size.x)
        );
        // Final Gradients
        // d/dx comes from Hermite deriv of this stage
        float3 gradX = float3(
            hermite1D_deriv(t.x, r[0].x, r_dX[0].x * size.x, r[1].x, r_dX[1].x * size.x),
            hermite1D_deriv(t.x, r[0].y, r_dX[0].y * size.x, r[1].y, r_dX[1].y * size.x),
            hermite1D_deriv(t.x, r[0].z, r_dX[0].z * size.x, r[1].z, r_dX[1].z * size.x)
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

ref<AdaptiveProbeVolume> AdaptiveProbeVolume::create(ref<Device> pDevice)
{
    return ref<AdaptiveProbeVolume>(new AdaptiveProbeVolume(pDevice));
}

AdaptiveProbeVolume::AdaptiveProbeVolume(ref<Device> pDevice) : mpDevice(pDevice) {}

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
    constrainHangingNodesHermite();
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
        //for (int k = 0; k < 6; ++k) gp.coarseNeighbors[k] = p.coarseNeighbors[k];
        //gp.pad3[0] = 0; // Padding maintenance
        //gp.pad3[1] = 0;

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
        //gc.position = c.position;
        //gc.pad = 0.0f;
        //gc.isValid = c.isValid ? 1 : 0;
        // Safety: ensure we don't read out of bounds if corner has fewer bands
        int numBands = std::min((int)c.shCoeffs.size(), 9);

        for (int i = 0; i < 9; ++i)
        {
            if (i < numBands)
            {
                gc.coeffs[i] = float4(c.shCoeffs[i], 0.0f);
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

    mpSeedProbeIndexBuffer = mpDevice->createStructuredBuffer(
        sizeof(int),
        (uint32_t)mSeedProbeIndices.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        mSeedProbeIndices.data()
    );
    mpSeedProbeIndexBuffer->setName("SeedProbeIndexBuffer");

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
