#include <map>
#include <array>
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
    //computeNeighbors();
    //constrainHangingNodes();
    //constrainHangingNodesHermite();
    // 1. Pack the nodes (tree topology), compact layout (GPUProbe in AdaptiveGridData.slang, 36 bytes).
    // Nodes are re-ordered breadth-first from the seed roots so that the 8 children of every node are adjacent:
    // child k = firstChild + k. Bounds are not stored; the shader derives them from the seed cell. The CPU-side
    // mProbes keep their file order (only the GPU copy is re-ordered).
    static_assert(sizeof(GPUProbe) == 36, "GPUProbe must be the 36-byte compact node");

    // Seed as uploaded: the file's seed grid, or one cell covering the root (node 0) for a grid without a seed.
    std::vector<int> seedRoots;
    if (mUseSeedGrid && !mSeedProbeIndices.empty())
    {
        mGpuSeedMinPoint = mSeedMinPoint;
        mGpuSeedCellSize = mSeedCellSize;
        mGpuSeedResolution = mSeedResolution;
        seedRoots = mSeedProbeIndices;
    }
    else
    {
        mGpuSeedMinPoint = mProbes.empty() ? float3(0.f) : mProbes[0].minPoint;
        mGpuSeedCellSize = mProbes.empty() ? float3(1.f) : mProbes[0].maxPoint - mProbes[0].minPoint;
        mGpuSeedResolution = uint3(1);
        seedRoots = { mProbes.empty() ? -1 : 0 };
    }

    // Lattice of the Oat seam handling: the finest node size (= volume extent / 2^maxLevel, as in PlacementComparison).
    mGpuLatticeCellSize = mGpuSeedCellSize;
    for (const Probe& p : mProbes)
    {
        const float3 s = p.maxPoint - p.minPoint;
        if (s.x > 0.f && s.y > 0.f && s.z > 0.f)
            mGpuLatticeCellSize = float3(std::min(mGpuLatticeCellSize.x, s.x), std::min(mGpuLatticeCellSize.y, s.y), std::min(mGpuLatticeCellSize.z, s.z));
    }

    std::vector<int> newIndex(mProbes.size(), -1);
    std::vector<int> order; // order[newIdx] = old index
    order.reserve(mProbes.size());
    for (int r : seedRoots)
        if (r >= 0 && r < (int)mProbes.size() && newIndex[r] < 0) { newIndex[r] = (int)order.size(); order.push_back(r); }

    uint32_t partialNodes = 0;
    for (size_t i = 0; i < order.size(); ++i)
    {
        const Probe& p = mProbes[order[i]];
        if (p.children[0] == -1) continue;
        bool complete = true;
        for (int k = 0; k < 8; ++k)
            if (p.children[k] < 0 || p.children[k] >= (int)mProbes.size() || newIndex[p.children[k]] >= 0) complete = false;
        if (!complete) { ++partialNodes; continue; } // treated as a leaf below
        for (int k = 0; k < 8; ++k) { newIndex[p.children[k]] = (int)order.size(); order.push_back(p.children[k]); }
    }
    if (partialNodes > 0)
        logWarning("[AdaptiveProbeVolume] {} nodes without 8 valid children are uploaded as leaves.", partialNodes);
    if (order.size() != mProbes.size())
        logWarning("[AdaptiveProbeVolume] {} of {} nodes are not reachable from the seed and are not uploaded.",
            mProbes.size() - order.size(), mProbes.size());

    std::vector<GPUProbe> gpuProbes;
    gpuProbes.reserve(order.size());
    for (int oldIdx : order)
    {
        const Probe& p = mProbes[oldIdx];
        GPUProbe gp;
        const int c0 = p.children[0];
        gp.firstChild = (c0 >= 0 && c0 < (int)newIndex.size()) ? newIndex[c0] : -1;
        // A node uploaded as a leaf (partial children) has newIndex[c0] = -1 or a child that is not at c0 + k.
        if (gp.firstChild >= 0)
            for (int k = 0; k < 8; ++k)
                if (newIndex[p.children[k]] != gp.firstChild + k) { gp.firstChild = -1; break; }
        for (int i = 0; i < 8; ++i) gp.corners[i] = p.corners[i];
        gpuProbes.push_back(gp);
    }

    // Oat seam directions for the shader (XAtlasSceneShader.slang, fetchOatSH). For every leaf, a mask of the 26
    // neighbour directions d in {-1, 0, 1}^3 (bit (dx+1)*9 + (dy+1)*3 + (dz+1)) whose same-size neighbour cell is
    // subdivided, i.e. holds leaves finer than this one that touch it. Stored in the leaf as firstChild = -1 - mask
    // (mask 0 = -1, a plain leaf). The shader does the seam lookups only at lattice nodes facing such a direction and
    // uses the plain leaf interpolation everywhere else (e.g. in every finest-level leaf).
    {
        const size_t nodeCount = gpuProbes.size();
        std::vector<float3> nodeMin(nodeCount, float3(0.f)), nodeSize(nodeCount, float3(0.f));
        for (uint32_t z = 0; z < mGpuSeedResolution.z; ++z)
            for (uint32_t y = 0; y < mGpuSeedResolution.y; ++y)
                for (uint32_t x = 0; x < mGpuSeedResolution.x; ++x)
                {
                    const size_t si = (size_t(z) * mGpuSeedResolution.y + y) * mGpuSeedResolution.x + x;
                    const int r = si < seedRoots.size() && seedRoots[si] >= 0 ? newIndex[seedRoots[si]] : -1;
                    if (r < 0) continue;
                    nodeMin[r] = mGpuSeedMinPoint + float3(float(x), float(y), float(z)) * mGpuSeedCellSize;
                    nodeSize[r] = mGpuSeedCellSize;
                }
        for (size_t i = 0; i < nodeCount; ++i) // parents come before their children (breadth-first order)
        {
            const int fc = gpuProbes[i].firstChild;
            if (fc < 0) continue;
            const float3 h = 0.5f * nodeSize[i];
            for (int k = 0; k < 8; ++k)
            {
                nodeMin[fc + k] = nodeMin[i] + float3((k & 4) ? h.x : 0.f, (k & 2) ? h.y : 0.f, (k & 1) ? h.z : 0.f);
                nodeSize[fc + k] = h;
            }
        }
        // True when the node of size ~s containing q (q = centre of a cell of size s) has children.
        auto subdividedAt = [&](const float3& q, float s) -> bool
        {
            const float3 rel = (q - mGpuSeedMinPoint) / mGpuSeedCellSize;
            if (rel.x < 0.f || rel.y < 0.f || rel.z < 0.f || rel.x >= float(mGpuSeedResolution.x) ||
                rel.y >= float(mGpuSeedResolution.y) || rel.z >= float(mGpuSeedResolution.z))
                return false;
            const size_t si = (size_t(rel.z) * mGpuSeedResolution.y + size_t(rel.y)) * mGpuSeedResolution.x + size_t(rel.x);
            int idx = si < seedRoots.size() && seedRoots[si] >= 0 ? newIndex[seedRoots[si]] : -1;
            while (idx >= 0)
            {
                const int fc = gpuProbes[idx].firstChild;
                if (fc < 0) return false;                          // a leaf as large as the cell or larger
                if (nodeSize[idx].x < 1.5f * s) return true;      // the same-size cell itself, subdivided
                const float3 cen = nodeMin[idx] + 0.5f * nodeSize[idx];
                idx = fc + ((q.x >= cen.x) ? 4 : 0) + ((q.y >= cen.y) ? 2 : 0) + ((q.z >= cen.z) ? 1 : 0);
            }
            return false;
        };
        uint32_t seamLeaves = 0;
        for (size_t i = 0; i < nodeCount; ++i)
        {
            if (gpuProbes[i].firstChild >= 0 || nodeSize[i].x <= 0.f) continue;
            const float3 s = nodeSize[i], cen = nodeMin[i] + 0.5f * s;
            int mask = 0;
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz)
                        if ((dx | dy | dz) != 0 && subdividedAt(cen + float3(float(dx), float(dy), float(dz)) * s, s.x))
                            mask |= 1 << ((dx + 1) * 9 + (dy + 1) * 3 + (dz + 1));
            gpuProbes[i].firstChild = -1 - mask;
            seamLeaves += mask != 0;
        }
        logInfo("[AdaptiveProbeVolume] {} leaves touch a finer leaf (Oat seam handling).", seamLeaves);

        // Hanging-node constraint (setHangingConstraint; PrecomputeSHCoefficients' constrainHangingCorners with the
        // seeded tree as uploaded): every corner on a face or edge of a coarser leaf, without being one of its corners,
        // takes that leaf's value there for the runtime interpolation, and the trilinear blend of its corner gradients.
        // Coarsest constraint first. The corners below are packed after this, so the GPU gets the constrained data and
        // the shader needs no seam handling (SH_SEAM_OAT 0).
        mConstrainedCorners = 0;
        if (mHangingConstraintInterp >= 0)
        {
            const int interp = mHangingConstraintInterp;
            const float3 lcell = mGpuLatticeCellSize;
            auto leafAt = [&](const float3& q) -> int
            {
                const float3 rel = (q - mGpuSeedMinPoint) / mGpuSeedCellSize;
                if (rel.x < 0.f || rel.y < 0.f || rel.z < 0.f || rel.x > float(mGpuSeedResolution.x) ||
                    rel.y > float(mGpuSeedResolution.y) || rel.z > float(mGpuSeedResolution.z))
                    return -1;
                const uint32_t cx = std::min(uint32_t(rel.x), mGpuSeedResolution.x - 1);
                const uint32_t cy = std::min(uint32_t(rel.y), mGpuSeedResolution.y - 1);
                const uint32_t cz = std::min(uint32_t(rel.z), mGpuSeedResolution.z - 1);
                const size_t si = (size_t(cz) * mGpuSeedResolution.y + cy) * mGpuSeedResolution.x + cx;
                int idx = si < seedRoots.size() && seedRoots[si] >= 0 ? newIndex[seedRoots[si]] : -1;
                while (idx >= 0 && gpuProbes[idx].firstChild >= 0)
                {
                    const float3 cen = nodeMin[idx] + 0.5f * nodeSize[idx];
                    idx = gpuProbes[idx].firstChild + ((q.x >= cen.x) ? 4 : 0) + ((q.y >= cen.y) ? 2 : 0) + ((q.z >= cen.z) ? 1 : 0);
                }
                return idx;
            };
            auto hermite1 = [&](float t, float v0, float s0, float v1, float s1) -> float
            {
                if (interp == 2)
                {
                    const float d = v1 - v0;
                    if (d == 0.f) { s0 = 0.f; s1 = 0.f; }
                    else
                    {
                        if (s0 * d < 0.f) s0 = 0.f;
                        if (s1 * d < 0.f) s1 = 0.f;
                        const float a = s0 / d, b = s1 / d, r2 = a * a + b * b;
                        if (r2 > 9.f) { const float k = 3.f / std::sqrt(r2); s0 *= k; s1 *= k; }
                    }
                }
                const float t2 = t * t, t3 = t2 * t;
                return (2.f * t3 - 3.f * t2 + 1.f) * v0 + (t3 - 2.f * t2 + t) * s0 + (-2.f * t3 + 3.f * t2) * v1 + (t3 - t2) * s1;
            };
            auto hermite3 = [&](float t, const float3& v0, const float3& s0, const float3& v1, const float3& s1)
            {
                return float3(hermite1(t, v0.x, s0.x, v1.x, s1.x), hermite1(t, v0.y, s0.y, v1.y, s1.y), hermite1(t, v0.z, s0.z, v1.z, s1.z));
            };
            // Value of one band of leaf L at x, as the shader interpolates it (stored gradients: .x = world Z,
            // .y = world X, .z = world Y).
            auto leafValue = [&](int L, const float3& x, size_t band) -> float3
            {
                const float3 s = nodeSize[L];
                float3 t = (x - nodeMin[L]) / s;
                t = float3(std::clamp(t.x, 0.f, 1.f), std::clamp(t.y, 0.f, 1.f), std::clamp(t.z, 0.f, 1.f));
                float3 v[8], gX[8], gY[8], gZ[8];
                for (int k = 0; k < 8; ++k)
                {
                    const Corner& ck = mCorners[gpuProbes[L].corners[k]];
                    v[k] = band < ck.shCoeffs.size() ? ck.shCoeffs[band] : float3(0.f);
                    const GradSHCoeff g = band < ck.shGradients.size() ? ck.shGradients[band] : GradSHCoeff{ float3(0.f), float3(0.f), float3(0.f) };
                    gX[k] = float3(g.r.y, g.g.y, g.b.y);
                    gY[k] = float3(g.r.z, g.g.z, g.b.z);
                    gZ[k] = float3(g.r.x, g.g.x, g.b.x);
                }
                if (interp == 1)
                {
                    float3 r(0.f);
                    for (int k = 0; k < 8; ++k)
                        r += (((k >> 2) & 1) ? t.x : 1.f - t.x) * (((k >> 1) & 1) ? t.y : 1.f - t.y) * ((k & 1) ? t.z : 1.f - t.z) * v[k];
                    return r;
                }
                float3 q[4], qX[4], qY[4];
                for (int i = 0; i < 4; ++i)
                {
                    q[i] = hermite3(t.z, v[2 * i], gZ[2 * i] * s.z, v[2 * i + 1], gZ[2 * i + 1] * s.z);
                    qX[i] = gX[2 * i] + (gX[2 * i + 1] - gX[2 * i]) * t.z;
                    qY[i] = gY[2 * i] + (gY[2 * i + 1] - gY[2 * i]) * t.z;
                }
                float3 r[2], rX[2];
                for (int i = 0; i < 2; ++i)
                {
                    r[i] = hermite3(t.y, q[2 * i], qY[2 * i] * s.y, q[2 * i + 1], qY[2 * i + 1] * s.y);
                    rX[i] = qX[2 * i] + (qX[2 * i + 1] - qX[2 * i]) * t.y;
                }
                return hermite3(t.x, r[0], rX[0] * s.x, r[1], rX[1] * s.x);
            };

            // Lattice coordinates of each corner identify its position; the build stores a few positions twice (separate
            // records with slightly different traced values): those copies take their mean first, and positions, not
            // record indices, decide what a leaf's corner is.
            std::vector<std::array<int, 3>> key(mCorners.size());
            std::map<std::array<int, 3>, std::vector<int>> copies;
            for (size_t i = 0; i < mCorners.size(); ++i)
            {
                const float3 r = (mCorners[i].position - mGpuSeedMinPoint) / lcell;
                key[i] = { int(std::lround(r.x)), int(std::lround(r.y)), int(std::lround(r.z)) };
                copies[key[i]].push_back(int(i));
            }
            for (const auto& [k, ids] : copies)
            {
                if (ids.size() < 2) continue;
                Corner mean = mCorners[ids[0]];
                for (size_t j = 1; j < ids.size(); ++j)
                {
                    const Corner& o = mCorners[ids[j]];
                    for (size_t b = 0; b < mean.shCoeffs.size() && b < o.shCoeffs.size(); ++b) mean.shCoeffs[b] += o.shCoeffs[b];
                    for (size_t b = 0; b < mean.shGradients.size() && b < o.shGradients.size(); ++b)
                    {
                        mean.shGradients[b].r += o.shGradients[b].r;
                        mean.shGradients[b].g += o.shGradients[b].g;
                        mean.shGradients[b].b += o.shGradients[b].b;
                    }
                }
                const float inv = 1.f / float(ids.size());
                for (auto& v : mean.shCoeffs) v = v * inv;
                for (auto& g : mean.shGradients) { g.r = g.r * inv; g.g = g.g * inv; g.b = g.b * inv; }
                for (int id : ids) { mCorners[id].shCoeffs = mean.shCoeffs; mCorners[id].shGradients = mean.shGradients; }
            }

            struct Item { int corner; int leaf; float size; };
            std::vector<Item> items;
            std::vector<char> seen(mCorners.size(), 0);
            for (size_t i = 0; i < nodeCount; ++i)
            {
                if (gpuProbes[i].firstChild >= 0) continue;
                for (int k = 0; k < 8; ++k)
                {
                    const int ci = gpuProbes[i].corners[k];
                    if (ci < 0 || ci >= int(mCorners.size()) || seen[ci]) continue;
                    seen[ci] = 1;
                    const float3 x = mCorners[ci].position;
                    int best = -1;
                    float bestSize = 0.f;
                    for (int d = 0; d < 8; ++d)
                    {
                        const float3 o = x + float3((d & 4) ? 1.f : -1.f, (d & 2) ? 1.f : -1.f, (d & 1) ? 1.f : -1.f) * (1e-3f * lcell);
                        const int nb = leafAt(o);
                        if (nb < 0) continue;
                        bool isCorner = false;
                        for (int j = 0; j < 8; ++j)
                            isCorner = isCorner || (gpuProbes[nb].corners[j] >= 0 && key[gpuProbes[nb].corners[j]] == key[ci]);
                        if (!isCorner && nodeSize[nb].x > bestSize) { best = nb; bestSize = nodeSize[nb].x; }
                    }
                    if (best >= 0) items.push_back({ ci, best, bestSize });
                }
            }
            std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.size > b.size; });
            for (const Item& it : items)
            {
                Corner& c = mCorners[it.corner];
                const size_t bands = c.shCoeffs.size();
                std::vector<float3> v(bands);
                for (size_t b = 0; b < bands; ++b) v[b] = leafValue(it.leaf, c.position, b);
                float3 t = (c.position - nodeMin[it.leaf]) / nodeSize[it.leaf];
                t = float3(std::clamp(t.x, 0.f, 1.f), std::clamp(t.y, 0.f, 1.f), std::clamp(t.z, 0.f, 1.f));
                std::vector<GradSHCoeff> g(bands, GradSHCoeff{ float3(0.f), float3(0.f), float3(0.f) });
                for (int k = 0; k < 8; ++k)
                {
                    const Corner& ck = mCorners[gpuProbes[it.leaf].corners[k]];
                    const float w = (((k >> 2) & 1) ? t.x : 1.f - t.x) * (((k >> 1) & 1) ? t.y : 1.f - t.y) * ((k & 1) ? t.z : 1.f - t.z);
                    for (size_t b = 0; b < bands && b < ck.shGradients.size(); ++b)
                    {
                        g[b].r += w * ck.shGradients[b].r;
                        g[b].g += w * ck.shGradients[b].g;
                        g[b].b += w * ck.shGradients[b].b;
                    }
                }
                c.shCoeffs = v;
                c.shGradients = g;
                ++mConstrainedCorners;
            }
            logInfo("[AdaptiveProbeVolume] Hanging-node constraint ({}): {} of {} corners constrained.",
                interp == 1 ? "trilinear" : interp == 2 ? "limited Hermite" : "Hermite", mConstrainedCorners, mCorners.size());
        }
    }

    std::vector<int> gpuSeedIndices(seedRoots.size());
    for (size_t i = 0; i < seedRoots.size(); ++i)
        gpuSeedIndices[i] = (seedRoots[i] >= 0 && seedRoots[i] < (int)newIndex.size()) ? newIndex[seedRoots[i]] : -1;

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
        (uint32_t)gpuSeedIndices.size(),
        ResourceBindFlags::ShaderResource,
        MemoryType::DeviceLocal,
        gpuSeedIndices.data()
    );
    mpSeedProbeIndexBuffer->setName("SeedProbeIndexBuffer");

    mGpuNodeBytes = uint64_t(gpuProbes.size()) * sizeof(GPUProbe);
    mGpuCornerBytes = uint64_t(gpuCorners.size()) * sizeof(GPUCorner);
    mGpuSeedBytes = uint64_t(gpuSeedIndices.size()) * sizeof(int);
    logInfo("[AdaptiveProbeVolume] Runtime memory: {} probes x {} B = {:.2f} MB, {} nodes x {} B = {:.2f} MB, seed {} B; total {:.2f} MB.",
        gpuCorners.size(), sizeof(GPUCorner), mGpuCornerBytes / 1e6, gpuProbes.size(), sizeof(GPUProbe), mGpuNodeBytes / 1e6,
        mGpuSeedBytes, (mGpuNodeBytes + mGpuCornerBytes + mGpuSeedBytes) / 1e6);

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

    const bool isV5 =
        header == "ADAPTIVE_GRID_V5_SEEDED";

    const bool isV6 =
        header == "ADAPTIVE_GRID_V6_EDGE_ADDED";

    if (!isV5 && !isV6)
    {
        logError(
            "Invalid AdaptiveProbeVolume file format: " +
            filename
        );
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
        int addedByEdgeAwareInt = 0;

        if (isV6)
        {
            in >> isLeafInt
                >> addedByEdgeAwareInt
                >> p.level;
        }
        else
        {
            // Old V5 files contain no edge-added classification.
            in >> isLeafInt
                >> p.level;

            addedByEdgeAwareInt = 0;
        }

        p.isLeaf =
            isLeafInt != 0;

        p.addedByEdgeAware =
            addedByEdgeAwareInt != 0;

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

    logInfo(
        "Successfully loaded AdaptiveProbeVolume from " +
        filename
    );
}
