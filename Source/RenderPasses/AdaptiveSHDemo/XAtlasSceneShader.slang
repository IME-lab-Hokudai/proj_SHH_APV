import Scene.Raster;
#include "AdaptiveGridData.slang"
#include "Utils/Math/MathConstants.slangh"

// ---- SH windowing test switch (see evaluateIrradiance) ----
#ifndef SH_WINDOW
#define SH_WINDOW 0       // 0 = off (original), 1 = Hanning, 2 = Lanczos
#endif
#ifndef SH_WINDOW_WIDTH
#define SH_WINDOW_WIDTH 4  // larger = weaker windowing
#endif

// ---- SH deringing at shading points (AdaptiveSHDemo "SH deringing (Sloan 2017)"; SHNonNegative.slang) ----
// 1: where the interpolated SH has a negative irradiance lobe, scale its bands 1-2 down just enough to remove it, the
// same for R, G and B (the probes themselves are windowed with Sloan's method when the grid is loaded).
#ifndef SH_NONNEG
#define SH_NONNEG 0
#endif
#include "SHNonNegative.slang"


// Shared material, xatlas, and glass rendering for both probe layouts.
SamplerState gLinearSampler;
Texture2D<float4> gOpaqueColor;
Texture2D<float> gOpaqueDepth;
StructuredBuffer<uint> gLiquidMaterials;
cbuffer GlassCB
{
    uint gGlassEnabled;
    uint gHasEnvironment;
    float gSceneScale;
    float gRefractionStrength;
    uint gLiquidStage;
    float gLiquidOpticalDepth;
};

float transmissionAmount(BasicMaterialData material)
{
    return saturate(float(material.specularTransmission) + float(material.diffuseTransmission));
}

// -----------------------------------------------------------------------------
// Bistro xatlas lightmap test resources.
//
// CPU side builds a dense lookup indexed by:
//     instance range offset + SV_PrimitiveID
//
// Each mapped triangle stores:
//     page index
//     xatlas UV at its three original triangle corners
// -----------------------------------------------------------------------------
Texture2D gBistroAtlasPages[10];

StructuredBuffer<uint2>  gBistroInstanceRanges; // x = offset, y = triangleCount
StructuredBuffer<uint>   gBistroPageIndices;    // UINT_MAX = not lightmapped
StructuredBuffer<float4> gBistroUV01;           // uv0.xy, uv1.xy
StructuredBuffer<float2> gBistroUV2;            // uv2.xy

cbuffer BistroLightmapCB
{
    uint gBistroLightmapEnabled;
    uint gBistroPageCount;
    uint gBistroMaxInstanceID;
    uint gBistroPadding;
};

#define MAX_SH_BASIS 9
#if !SH_GRID_UNIFORM
#define OCTREE_MAX_DEPTH 9

StructuredBuffer<GPUCorner> gCornerBuffer;
StructuredBuffer<GPUProbe> gProbeBuffer;
StructuredBuffer<int> gSeedProbeIndices;

cbuffer SeedGridCB
{
    uint gUseSeedGrid;
    float3 gSeedMinPoint;
    float3 gSeedCellSize;
    uint3 gSeedResolution;
    float3 gLatticeCellSize; // finest leaf size of the grid (Oat seam lattice), set by AdaptiveSHDemo
};

#endif

bool getBistroLightmapLookup(
    VSOut vsOut,
    uint triangleIndex,
    float3 barycentrics,
    out uint pageIndex,
    out float2 lightmapUV
)
{
    pageIndex = 0xffffffffu;
    lightmapUV = float2(0.0f);

    if (gBistroLightmapEnabled == 0)
    {
        return false;
    }

    const uint instanceID = vsOut.instanceID.index;

    if (instanceID > gBistroMaxInstanceID)
    {
        return false;
    }

    const uint2 range = gBistroInstanceRanges[instanceID];

    if (triangleIndex >= range.y)
    {
        return false;
    }

    const uint lookupIndex = range.x + triangleIndex;
    pageIndex = gBistroPageIndices[lookupIndex];

    if (pageIndex == 0xffffffffu ||
        pageIndex >= gBistroPageCount)
    {
        return false;
    }

    const float4 uv01 = gBistroUV01[lookupIndex];
    const float2 uv2  = gBistroUV2[lookupIndex];

    const float3 bary = barycentrics;
    const float2 uv0 = uv01.xy;
    const float2 uv1 = uv01.zw;

    lightmapUV =
        bary.x * uv0 +
        bary.y * uv1 +
        bary.z * uv2;

    // XAtlasUVRaster maps xatlas UV.y directly into clip-space Y.
    // The render-target / texture sampling convention therefore requires
    // the vertical flip when sampling the saved EXR.
    lightmapUV.y = 1.0f - lightmapUV.y;

    return true;
}

float3 sampleBistroAtlas(
    uint pageIndex,
    float2 uv
)
{
    uv = saturate(uv);
    if (pageIndex >= gBistroPageCount) return float3(0.f);

    // Constant descriptor indices also handle pixels from different pages
    // within the same wave without requiring non-uniform resource indexing.
    switch (pageIndex)
    {
    case 0: return gBistroAtlasPages[0].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 1: return gBistroAtlasPages[1].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 2: return gBistroAtlasPages[2].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 3: return gBistroAtlasPages[3].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 4: return gBistroAtlasPages[4].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 5: return gBistroAtlasPages[5].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 6: return gBistroAtlasPages[6].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 7: return gBistroAtlasPages[7].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 8: return gBistroAtlasPages[8].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    case 9: return gBistroAtlasPages[9].SampleLevel(gLinearSampler, uv, 0.f).rgb;
    }

    return float3(0.0f);
}


// ------------------------------------------------------------------
// Interpolation mode, as in the placement evaluation (PlacementComparison.cpp). AdaptiveSHDemo sets both from its UI.
//   SH_INTERP:   0 = tricubic Hermite, 1 = trilinear, 2 = Hermite with Fritsch-Carlson slope limiting.
//   SH_SEAM_OAT: adaptive grids only. 1 = exact Oat (GDC 2005) seam handling: in a lattice cell of the finest level
//                that touches a finer leaf, interpolate on that lattice cell, with node values (and gradients) from the
//                smallest leaf touching each node. 0 = interpolate each leaf on its own (seams between levels).
// ------------------------------------------------------------------
#ifndef SH_INTERP
#define SH_INTERP 2
#endif
#ifndef SH_SEAM_OAT
#define SH_SEAM_OAT 1
#endif

// ------------------------------------------------------------------
// 1D Hermite Interpolation, per channel
// t: Normalized distance (0..1)
// v0, v1: Values at endpoints
// s0, s1: Slopes at endpoints (Must be scaled by segment length!)
// With SH_INTERP == 2 the slopes are limited first, exactly as hermite1DLimited in AdaptiveProbeVolume.cpp
// (PrecomputeSHCoefficients): a slope whose sign disagrees with the secant d = v1 - v0 (or any slope when d = 0) is
// set to 0, and both are scaled down when (s0/d)^2 + (s1/d)^2 > 9.
// ------------------------------------------------------------------
float hermite1DScalar(float t, float v0, float s0, float v1, float s1)
{
#if SH_INTERP == 2
    float d = v1 - v0;
    if (d == 0.f)
    {
        s0 = 0.f;
        s1 = 0.f;
    }
    else
    {
        if (s0 * d < 0.f) s0 = 0.f;
        if (s1 * d < 0.f) s1 = 0.f;
        float a = s0 / d, b = s1 / d;
        float r2 = a * a + b * b;
        if (r2 > 9.f)
        {
            float k = 3.f / sqrt(r2);
            s0 *= k;
            s1 *= k;
        }
    }
#endif
    float t2 = t * t;
    float t3 = t2 * t;

    // Standard Cubic Hermite Basis functions
    float h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    float h10 = t3 - 2.0 * t2 + t;
    float h01 = -2.0 * t3 + 3.0 * t2;
    float h11 = t3 - t2;

    return h00 * v0 + h10 * s0 + h01 * v1 + h11 * s1;
}

float3 hermite1D(float t, float3 v0, float3 s0, float3 v1, float3 s1)
{
    return float3(
        hermite1DScalar(t, v0.x, s0.x, v1.x, s1.x),
        hermite1DScalar(t, v0.y, s0.y, v1.y, s1.y),
        hermite1DScalar(t, v0.z, s0.z, v1.z, s1.z));
}

// Progressive tricubic Hermite of one band in one cell: along Falcor Z (corner bit 1), then Y (bit 2), then X (bit 4);
// the transverse gradients are interpolated linearly. v: corner values; gX, gY, gZ: corner gradients in world axes;
// size: cell size (world-space derivatives are scaled to the 0..1 parameter).
float3 hermiteCellBand(float3 t, float3 size, float3 v[8], float3 gX[8], float3 gY[8], float3 gZ[8])
{
    float3 q[4];
    float3 q_dX[4];
    float3 q_dY[4];
    [unroll]
    for (int i = 0; i < 4; ++i)
    {
        int i0 = 2 * i;
        int i1 = 2 * i + 1;
        q[i] = hermite1D(t.z, v[i0], gZ[i0] * size.z, v[i1], gZ[i1] * size.z);
        q_dX[i] = lerp(gX[i0], gX[i1], t.z);
        q_dY[i] = lerp(gY[i0], gY[i1], t.z);
    }

    float3 r[2];
    float3 r_dX[2];
    [unroll]
    for (int i = 0; i < 2; ++i)
    {
        int i0 = 2 * i;
        int i1 = 2 * i + 1;
        r[i] = hermite1D(t.y, q[i0], q_dY[i0] * size.y, q[i1], q_dY[i1] * size.y);
        r_dX[i] = lerp(q_dX[i0], q_dX[i1], t.y);
    }

    return hermite1D(t.x, r[0], r_dX[0] * size.x, r[1], r_dX[1] * size.x);
}

// Trilinear of one band in one cell (corner bits x = 4, y = 2, z = 1).
float3 trilinearCellBand(float3 t, float3 v[8])
{
    float3 z0 = lerp(v[0], v[1], t.z);
    float3 z1 = lerp(v[2], v[3], t.z);
    float3 z2 = lerp(v[4], v[5], t.z);
    float3 z3 = lerp(v[6], v[7], t.z);
    float3 y0 = lerp(z0, z1, t.y);
    float3 y1 = lerp(z2, z3, t.y);
    return lerp(y0, y1, t.x);
}

#if SH_GRID_UNIFORM
#include "UniformGridSampling.slang"
#else
// Stored gradient convention: .x = world Z, .y = world X, .z = world Y.
void cornerGradientsBand(int corner, int band, out float3 gX, out float3 gY, out float3 gZ)
{
    float4 r = gCornerBuffer[corner].gradR[band];
    float4 g = gCornerBuffer[corner].gradG[band];
    float4 b = gCornerBuffer[corner].gradB[band];
    gX = float3(r.y, g.y, b.y);
    gY = float3(r.z, g.z, b.z);
    gZ = float3(r.x, g.x, b.x);
}

// One band of a leaf, interpolated at pos (Hermite, limited Hermite or trilinear, by SH_INTERP).
// nodeMin, size: bounds of the leaf, from traverseOctree (not stored in the node).
float3 leafBand(int leaf, float3 nodeMin, float3 size, float3 pos, int band)
{
    float3 t = saturate((pos - nodeMin) / size);
    int cIdx[8];
    float3 v[8];
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        cIdx[k] = gProbeBuffer[leaf].corners[k];
        v[k] = gCornerBuffer[cIdx[k]].coeffs[band].rgb;
    }
#if SH_INTERP == 1
    return trilinearCellBand(t, v);
#else
    float3 gX[8];
    float3 gY[8];
    float3 gZ[8];
    [unroll]
    for (int k = 0; k < 8; ++k) cornerGradientsBand(cIdx[k], band, gX[k], gY[k], gZ[k]);
    return hermiteCellBand(t, size, v, gX, gY, gZ);
#endif
}

// Gradients of one band of a leaf at pos: trilinear blend of its 8 corner gradients (as leafGradients in the evaluation).
void leafGradientsBand(int leaf, float3 nodeMin, float3 size, float3 pos, int band, out float3 gX, out float3 gY, out float3 gZ)
{
    float3 t = saturate((pos - nodeMin) / size);
    gX = float3(0.f);
    gY = float3(0.f);
    gZ = float3(0.f);
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        float w = (((k >> 2) & 1) != 0 ? t.x : 1.f - t.x) * (((k >> 1) & 1) != 0 ? t.y : 1.f - t.y) * ((k & 1) != 0 ? t.z : 1.f - t.z);
        float3 cx, cy, cz;
        cornerGradientsBand(gProbeBuffer[leaf].corners[k], band, cx, cy, cz);
        gX += w * cx;
        gY += w * cy;
        gZ += w * cz;
    }
}

// One band of a leaf at pos (Hermite or limited Hermite) together with its gradients there (trilinear blend of the
// corner gradients, as leafGradientsBand), loading each corner once. Used for the Oat lattice nodes.
float3 leafBandAndGradients(int leaf, float3 nodeMin, float3 size, float3 pos, int band, out float3 oX, out float3 oY, out float3 oZ)
{
    float3 t = saturate((pos - nodeMin) / size);
    float3 v[8];
    float3 gX[8];
    float3 gY[8];
    float3 gZ[8];
    oX = float3(0.f);
    oY = float3(0.f);
    oZ = float3(0.f);
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        int ci = gProbeBuffer[leaf].corners[k];
        v[k] = gCornerBuffer[ci].coeffs[band].rgb;
        cornerGradientsBand(ci, band, gX[k], gY[k], gZ[k]);
        float w = (((k >> 2) & 1) != 0 ? t.x : 1.f - t.x) * (((k >> 1) & 1) != 0 ? t.y : 1.f - t.y) * ((k & 1) != 0 ? t.z : 1.f - t.z);
        oX += w * gX[k];
        oY += w * gY[k];
        oZ += w * gZ[k];
    }
    return hermiteCellBand(t, size, v, gX, gY, gZ);
}

void fetchLeafSH(int leaf, float3 nodeMin, float3 size, float3 pos, out float3 finalCoeffs[9])
{
    [unroll]
    for (int band = 0; band < 9; ++band) finalCoeffs[band] = leafBand(leaf, nodeMin, size, pos, band);
}

// Seed cell containing pos: its root node, and its bounds (cellMin, cellSize). Every grid has a seed on the GPU
// (uploadToGPU gives a grid without one a single 1 x 1 x 1 seed cell = the root), so gUseSeedGrid is always 1.
// A point on the upper face of the volume counts as inside the last cell.
int findSeedProbe(float3 pos, out float3 cellMin, out float3 cellSize)
{
    cellMin = float3(0.f);
    cellSize = gSeedCellSize;
    if (gUseSeedGrid == 0)
    {
        return -1;
    }

    float3 rel = (pos - gSeedMinPoint) / gSeedCellSize;
    if (any(rel < 0.f) || any(rel > float3(gSeedResolution)))
    {
        return -1;
    }
    int3 cell = min(int3(floor(rel)), int3(gSeedResolution) - 1);

    cellMin = gSeedMinPoint + float3(cell) * gSeedCellSize;
    uint idx = (uint(cell.z) * gSeedResolution.y + uint(cell.y)) * gSeedResolution.x + uint(cell.x);
    return gSeedProbeIndices[idx];
}

// Leaf containing pos, and its bounds. Bounds are not stored per node: they start as the seed cell and are
// halved at every level, moving to the half that contains pos (octant bits x = 4, y = 2, z = 1, as in the build).
int traverseOctree(float3 pos, out float3 nodeMin, out float3 nodeSize)
{
    int currentIdx = findSeedProbe(pos, nodeMin, nodeSize);
    if (currentIdx < 0)
    {
        return -1;
    }

    for (int i = 0; i < OCTREE_MAX_DEPTH; ++i)
    {
        int firstChild = gProbeBuffer[currentIdx].firstChild;
        if (firstChild < 0)
        {
            return currentIdx; // leaf
        }

        float3 halfSize = 0.5 * nodeSize;
        float3 center = nodeMin + halfSize;
        int octant = 0;
        if (pos.x >= center.x) { octant |= 4; nodeMin.x = center.x; }
        if (pos.y >= center.y) { octant |= 2; nodeMin.y = center.y; }
        if (pos.z >= center.z) { octant |= 1; nodeMin.z = center.z; }
        nodeSize = halfSize;

        currentIdx = firstChild + octant;
    }

    return currentIdx;
}

// The finest leaf whose closure contains x: look just off x in each of the 8 diagonal directions (same order and
// tie rule as smallestLeafAt in the evaluation: the first of the finest wins).
int smallestLeafAt(float3 x, out float3 leafMin, out float3 leafSize)
{
    int best = -1;
    float bestSize = 3.402823e38f;
    leafMin = float3(0.f);
    leafSize = float3(0.f);
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        float3 o = x + float3((k & 4) != 0 ? 1.f : -1.f, (k & 2) != 0 ? 1.f : -1.f, (k & 1) != 0 ? 1.f : -1.f) * (1e-3f * gLatticeCellSize);
        float3 mn, sz;
        int nb = traverseOctree(o, mn, sz);
        if (nb >= 0 && sz.x < 0.75f * bestSize)
        {
            best = nb;
            bestSize = sz.x;
            leafMin = mn;
            leafSize = sz;
        }
    }
    return best;
}

#if SH_SEAM_OAT
// The finest leaf strictly finer than maxSize whose closure contains the lattice node x, or -1. Looks just off x in the
// 8 diagonal directions like smallestLeafAt, but skips the directions that stay inside the current leaf (side 0 or
// opposite the node's side on an axis) or lead into a neighbour region without finer leaves (finerMask), since neither
// can hold a finer leaf. Same choice as smallestLeafAt restricted to finer leaves (first of the finest).
int finerLeafAt(float3 x, int3 side, uint finerMask, float maxSize, out float3 leafMin, out float3 leafSize)
{
    int best = -1;
    float bestSize = 0.75f * maxSize;
    leafMin = float3(0.f);
    leafSize = float3(0.f);
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        int3 s = int3((k & 4) != 0 ? 1 : -1, (k & 2) != 0 ? 1 : -1, (k & 1) != 0 ? 1 : -1);
        int3 d = int3(s.x == side.x ? side.x : 0, s.y == side.y ? side.y : 0, s.z == side.z ? side.z : 0);
        if (all(d == int3(0))) continue;
        uint bit = uint((d.x + 1) * 9 + (d.y + 1) * 3 + (d.z + 1));
        if (((finerMask >> bit) & 1u) == 0u) continue;
        float3 mn, sz;
        int nb = traverseOctree(x + float3(s) * (1e-3f * gLatticeCellSize), mn, sz);
        if (nb >= 0 && sz.x < bestSize)
        {
            best = nb;
            bestSize = sz.x;
            leafMin = mn;
            leafSize = sz;
        }
    }
    return best;
}
#endif

#if SH_SEAM_OAT
// Exact Oat seam handling (port of interpAdaptive in PlacementComparison.cpp). The lattice is the finest level of the
// grid, starting at the volume minimum (= gSeedMinPoint). In the lattice cell containing pos (kept inside the leaf),
// each node on the leaf boundary takes the smallest leaf touching it; if any node belongs to a strictly finer leaf, pos
// is interpolated on the lattice cell from those nodes, otherwise the leaf's own interpolation is used. Only nodes that
// face a neighbour direction flagged in the leaf's finer mask are looked up; leaves without one skip all of it. Neighbouring
// lattice cells share their nodes, so faces between leaves are continuous (exactly for trilinear).
void fetchOatSH(int leaf, float3 leafMin, float3 leafSize, float3 pos, out float3 finalCoeffs[9])
{
    // Leaves store which neighbour directions hold a finer leaf (AdaptiveProbeVolume::uploadToGPU): none -> no seam.
    const uint finerMask = uint(-1 - gProbeBuffer[leaf].firstChild);
    if (finerMask == 0u)
    {
        fetchLeafSH(leaf, leafMin, leafSize, pos, finalCoeffs);
        return;
    }
    float3 lcell = gLatticeCellSize;
    float3 volMin = gSeedMinPoint;
    int3 lo = int3(round((leafMin - volMin) / lcell));
    int3 n = max(int3(1), int3(round(leafSize / lcell)));
    int3 ci = clamp(int3(floor((pos - volMin) / lcell)), lo, lo + n - 1);

    int nodeLeaf[8];
    float3 nodeMin[8];
    float3 nodeSize[8];
    float3 nodePos[8];
    bool transition = false;
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        int3 g = ci + int3((k >> 2) & 1, (k >> 1) & 1, k & 1);
        nodePos[k] = volMin + float3(g) * lcell;
        nodeLeaf[k] = leaf;
        nodeMin[k] = leafMin;
        nodeSize[k] = leafSize;
        // Side of the leaf the node lies on, per axis (-1 min face, +1 max face, 0 inside). The neighbour directions
        // the node touches combine each axis' side with 0; only those flagged in finerMask can hold a finer leaf.
        int3 side = int3(g.x == lo.x ? -1 : (g.x == lo.x + n.x ? 1 : 0),
                         g.y == lo.y ? -1 : (g.y == lo.y + n.y ? 1 : 0),
                         g.z == lo.z ? -1 : (g.z == lo.z + n.z ? 1 : 0));
        if (any(side != int3(0)))
        {
            // Only a strictly finer leaf changes the node (an equal-level neighbour shares this leaf's face values),
            // as in the evaluation (PlacementComparison.cpp).
            float3 mn, sz;
            int nb = finerLeafAt(nodePos[k], side, finerMask, leafSize.x, mn, sz);
            if (nb >= 0)
            {
                nodeLeaf[k] = nb;
                nodeMin[k] = mn;
                nodeSize[k] = sz;
                transition = true;
            }
        }
    }

    if (!transition)
    {
        fetchLeafSH(leaf, leafMin, leafSize, pos, finalCoeffs);
        return;
    }

    float3 t = saturate((pos - nodePos[0]) / lcell);
    // Node data, cheapest source first (same values as interpolating each node's leaf at the node):
    //   node in this leaf      -> this leaf's interpolation at the node; its 8 corners are loaded once per band for all
    //                             such nodes (tOwn = node position in the leaf);
    //   node = a corner of its (finer) leaf, the usual case -> that probe corner directly (nodeCorner);
    //   otherwise              -> the node's leaf interpolated at the node.
    float3 tOwn[8];
    int nodeCorner[8];
    [unroll]
    for (int k = 0; k < 8; ++k)
    {
        tOwn[k] = saturate((nodePos[k] - leafMin) / leafSize);
        nodeCorner[k] = -1;
        if (nodeLeaf[k] != leaf)
        {
            float3 u = (nodePos[k] - nodeMin[k]) / nodeSize[k];
            float3 r = round(u);
            if (all(abs(u - r) < 1e-3f) && all(r >= 0.f) && all(r <= 1.f))
                nodeCorner[k] = gProbeBuffer[nodeLeaf[k]].corners[(r.x > 0.5f ? 4 : 0) | (r.y > 0.5f ? 2 : 0) | (r.z > 0.5f ? 1 : 0)];
        }
    }
    int own[8];
    [unroll]
    for (int k = 0; k < 8; ++k) own[k] = gProbeBuffer[leaf].corners[k];
    [unroll]
    for (int band = 0; band < 9; ++band)
    {
        float3 ov[8];
        [unroll]
        for (int k = 0; k < 8; ++k) ov[k] = gCornerBuffer[own[k]].coeffs[band].rgb;
        float3 v[8];
#if SH_INTERP == 1
        [unroll]
        for (int k = 0; k < 8; ++k)
        {
            if (nodeLeaf[k] == leaf) v[k] = trilinearCellBand(tOwn[k], ov);
            else if (nodeCorner[k] >= 0) v[k] = gCornerBuffer[nodeCorner[k]].coeffs[band].rgb;
            else v[k] = leafBand(nodeLeaf[k], nodeMin[k], nodeSize[k], nodePos[k], band);
        }
        finalCoeffs[band] = trilinearCellBand(t, v);
#else
        float3 ogX[8];
        float3 ogY[8];
        float3 ogZ[8];
        [unroll]
        for (int k = 0; k < 8; ++k) cornerGradientsBand(own[k], band, ogX[k], ogY[k], ogZ[k]);
        float3 gX[8];
        float3 gY[8];
        float3 gZ[8];
        [unroll]
        for (int k = 0; k < 8; ++k)
        {
            if (nodeLeaf[k] == leaf)
            {
                v[k] = hermiteCellBand(tOwn[k], leafSize, ov, ogX, ogY, ogZ);
                gX[k] = float3(0.f);
                gY[k] = float3(0.f);
                gZ[k] = float3(0.f);
                [unroll]
                for (int j = 0; j < 8; ++j)
                {
                    float w = (((j >> 2) & 1) != 0 ? tOwn[k].x : 1.f - tOwn[k].x) * (((j >> 1) & 1) != 0 ? tOwn[k].y : 1.f - tOwn[k].y) *
                              ((j & 1) != 0 ? tOwn[k].z : 1.f - tOwn[k].z);
                    gX[k] += w * ogX[j];
                    gY[k] += w * ogY[j];
                    gZ[k] += w * ogZ[j];
                }
            }
            else if (nodeCorner[k] >= 0)
            {
                v[k] = gCornerBuffer[nodeCorner[k]].coeffs[band].rgb;
                cornerGradientsBand(nodeCorner[k], band, gX[k], gY[k], gZ[k]);
            }
            else v[k] = leafBandAndGradients(nodeLeaf[k], nodeMin[k], nodeSize[k], nodePos[k], band, gX[k], gY[k], gZ[k]);
        }
        finalCoeffs[band] = hermiteCellBand(t, lcell, v, gX, gY, gZ);
#endif
    }
}
#endif

bool sampleGridSH(float3 pos, out float3 coeffs[MAX_SH_BASIS])
{
    [unroll]
    for (int band = 0; band < MAX_SH_BASIS; ++band) coeffs[band] = float3(0.f);
    float3 nodeMin, nodeSize;
    int probeIdx = traverseOctree(pos, nodeMin, nodeSize);
    if (probeIdx < 0) return false;
#if SH_SEAM_OAT
    fetchOatSH(probeIdx, nodeMin, nodeSize, pos, coeffs);
#else
    fetchLeafSH(probeIdx, nodeMin, nodeSize, pos, coeffs);
#endif
    return true;
}
#endif

// Band weights of the irradiance (cosine lobe, Sloan 2002), with the optional SH window: used by evaluateIrradiance
// and by the non-negative blend (SHNonNegative.slang), so both see the same irradiance.
float3 shIrradianceBandWeights()
{
    // Cosine lobe coefficients for irradiance (Peter-Pike Sloan 2002)
    float A[3] = { 3.141593f, 2.094395f, 0.785398f };

    // SH windowing (de-ringing test). Scales band l by sigma_l before evaluation.
    // Per-band constants, so this equals windowing every probe's coefficients
    // and gradients: consistent with the Hermite/trilinear interpolation.
    // No rebake needed. Band 0 is untouched (mean irradiance preserved).
    //   SH_WINDOW 0 = off (original), 1 = Hanning, 2 = Lanczos
    //   SH_WINDOW_WIDTH w: larger = weaker. Hanning w=4: (1, 0.854, 0.5); w=3: (1, 0.75, 0.25)
    //                                        Lanczos w=3: (1, 0.827, 0.413)
#if SH_WINDOW != 0
    [unroll]
    for (int l = 1; l < 3; ++l)
    {
        const float t = M_PI * float(l) / float(SH_WINDOW_WIDTH);
#if SH_WINDOW == 1
        const float sigma = 0.5f * (1.f + cos(t));
#else
        const float sigma = sin(t) / t;
#endif
        A[l] *= sigma;
    }
#endif
    return float3(A[0], A[1], A[2]);
}

// Evaluate 9 SH basis functions at direction n
float3 evaluateIrradiance(float3 n, float3 shCoeffs[MAX_SH_BASIS])
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

    float3 Aw = shIrradianceBandWeights();
    float A[3] = { Aw.x, Aw.y, Aw.z };


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
     //return max(result, float3(0.f));
    return result;
}

// SH deringing at a shading point (SH_NONNEG, SHNonNegative.slang): the probes are windowed when the grid is loaded
// (Sloan 2017); where Hermite interpolation still overshoots below zero, bands 1-2 are scaled down just enough to remove
// it, the same for R, G and B.
void applyNonNegativeIrradiance(inout float3 sh[MAX_SH_BASIS])
{
    shDeringShadingPoint(sh, shIrradianceBandWeights());
}

VSOut vsMain(VSIn vIn)
{
    VSOut vsOut = defaultVS(vIn);
    return vsOut;
}

float4 psMain(VSOut vsOut, uint triangleIndex : SV_PrimitiveID, float3 barycentrics : SV_Barycentrics) : SV_TARGET
{
    let materialLod = ExplicitGradientTextureSampler(ddx(vsOut.texC), ddy(vsOut.texC));
    if (alphaTest(vsOut, triangleIndex, materialLod)) discard;

    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);

    BasicMaterialData material = gScene.materials.getBasicMaterialData(vsOut.materialID);
    if (gGlassEnabled != 0 && transmissionAmount(material) > 0.01f) discard;
    let header = gScene.materials.getMaterialHeader(vsOut.materialID);
    SamplerState materialSampler = gScene.materials.getTextureSampler(header.getDefaultTextureSamplerID());
    float3 albedo = gScene.materials.sampleTexture(
        material.texBaseColor, materialSampler, vsOut.texC,
        (float4)material.baseColor, materialLod).rgb;
    float3 emission = gScene.materials.sampleTexture(
        material.texEmissive, materialSampler, vsOut.texC,
        float4(material.emissive, 1.f), materialLod).rgb * material.emissiveFactor;

    uint lightmapPage;
    float2 lightmapUV;
    float3 litColor;
    if (getBistroLightmapLookup(vsOut, triangleIndex, barycentrics, lightmapPage, lightmapUV))
    {
        // Atlas stores E/pi; receiver albedo is applied here.
        litColor = albedo * sampleBistroAtlas(lightmapPage, lightmapUV);
    }
    else
    {
        // Geometry added after the atlas bake uses the selected SH grid.
        float3 finalSH[MAX_SH_BASIS];
        if (sampleGridSH(shadingData.posW, finalSH))
        {
#if SH_NONNEG
            applyNonNegativeIrradiance(finalSH);
#endif
            float3 irradiance = max(
                evaluateIrradiance(shadingData.frame.N, finalSH),
                float3(0.f)
            );
            litColor = albedo * irradiance / M_PI;
        }
        else
        {
            // Preserve the diagnostic color for geometry outside the grid.
            litColor = float3(1.f, 0.f, 1.f);
        }
    }
    return float4(litColor + emission, 1.f);
}

// View dump for the offline image evaluation (AdaptiveSHDemo "Dump view"): per pixel, the data needed to shade it with
// any probe grid offline. Same surface data as psMain: shading position and normal, albedo, emission; pos.w = 1 where
// psMain shades with the SH grid (geometry without a lightmap, e.g. the dynamic object), 0 elsewhere.
struct ViewDumpOut
{
    float4 pos : SV_Target0;
    float4 normal : SV_Target1;
    float4 albedo : SV_Target2;
    float4 emission : SV_Target3;
};

ViewDumpOut psDump(VSOut vsOut, uint triangleIndex : SV_PrimitiveID, float3 barycentrics : SV_Barycentrics)
{
    let materialLod = ExplicitGradientTextureSampler(ddx(vsOut.texC), ddy(vsOut.texC));
    if (alphaTest(vsOut, triangleIndex, materialLod)) discard;

    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);

    BasicMaterialData material = gScene.materials.getBasicMaterialData(vsOut.materialID);
    if (gGlassEnabled != 0 && transmissionAmount(material) > 0.01f) discard;
    let header = gScene.materials.getMaterialHeader(vsOut.materialID);
    SamplerState materialSampler = gScene.materials.getTextureSampler(header.getDefaultTextureSamplerID());
    float3 albedo = gScene.materials.sampleTexture(
        material.texBaseColor, materialSampler, vsOut.texC,
        (float4)material.baseColor, materialLod).rgb;
    float3 emission = gScene.materials.sampleTexture(
        material.texEmissive, materialSampler, vsOut.texC,
        float4(material.emissive, 1.f), materialLod).rgb * material.emissiveFactor;

    uint lightmapPage;
    float2 lightmapUV;
    bool lightmapped = getBistroLightmapLookup(vsOut, triangleIndex, barycentrics, lightmapPage, lightmapUV);

    ViewDumpOut o;
    o.pos = float4(shadingData.posW, lightmapped ? 0.f : 1.f);
    o.normal = float4(shadingData.frame.N, 0.f);
    o.albedo = float4(albedo, 0.f);
    o.emission = float4(emission, 0.f);
    return o;
}

bool projectGlassPoint(float3 p, out float2 uv, out float depth)
{
    float4 clip = mul(gScene.camera.getViewProj(), float4(p, 1.f));
    if (clip.w <= 1e-5f) { uv = 0.f; depth = 0.f; return false; }
    float3 ndc = clip.xyz / clip.w;
    uv = ndc.xy * float2(0.5f, -0.5f) + 0.5f;
    depth = ndc.z;
    return all(uv > 0.f) && all(uv < 1.f) && depth > 0.f && depth < 1.f;
}

float3 glassBackground(float2 uv, float roughness)
{
    uint width, height;
    gOpaqueColor.GetDimensions(width, height);
    float2 texel = 1.f / float2(width, height);
    float2 lo = 0.5f * texel, hi = 1.f - lo;
    float2 radius = texel * (8.f * roughness * roughness);
    float3 color = 2.f * gOpaqueColor.SampleLevel(gLinearSampler, clamp(uv, lo, hi), 0.f).rgb;
    color += gOpaqueColor.SampleLevel(gLinearSampler, clamp(uv + float2(radius.x, 0.f), lo, hi), 0.f).rgb;
    color += gOpaqueColor.SampleLevel(gLinearSampler, clamp(uv - float2(radius.x, 0.f), lo, hi), 0.f).rgb;
    color += gOpaqueColor.SampleLevel(gLinearSampler, clamp(uv + float2(0.f, radius.y), lo, hi), 0.f).rgb;
    color += gOpaqueColor.SampleLevel(gLinearSampler, clamp(uv - float2(0.f, radius.y), lo, hi), 0.f).rgb;
    return color / 6.f;
}

float viewDistanceAt(float2 uv, float depth)
{
    float4 p = mul(gScene.camera.data.invViewProj, float4(uv * float2(2.f, -2.f) + float2(-1.f, 1.f), depth, 1.f));
    return dot(p.xyz / p.w - gScene.camera.getPosition(), normalize(gScene.camera.data.cameraW));
}

bool traceGlassReflection(float3 origin, float3 direction, out float2 hitUV)
{
    uint width, height;
    gOpaqueDepth.GetDimensions(width, height);
    float previous = -1.f;
    float previousT = 0.f;
    hitUV = float2(0.f);
    // Deterministic screen-space marching; no stochastic noise or history.
    for (uint i = 1; i <= 24; ++i)
    {
        float s = float(i) / 24.f;
        float t = gSceneScale * 0.75f * s * s;
        float2 uv; float depth;
        float3 p = origin + direction * t;
        if (!projectGlassPoint(p, uv, depth)) return false;
        uint2 pixel = min(uint2(uv * float2(width, height)), uint2(width - 1, height - 1));
        float sceneDepth = gOpaqueDepth[pixel];
        if (sceneDepth >= 1.f) { previous = -1.f; previousT = t; continue; }
        float difference = viewDistanceAt(uv, depth) - viewDistanceAt(uv, sceneDepth);
        if (previous < 0.f && difference >= 0.f && difference < 2.f * (t - previousT) + 0.002f * gSceneScale)
        {
            // Refine the crossing to reduce visible marching bands.
            float lo = previousT, hi = t;
            for (uint j = 0; j < 4; ++j)
            {
                float mid = 0.5f * (lo + hi);
                float2 q; float z;
                if (!projectGlassPoint(origin + direction * mid, q, z)) break;
                uint2 qp = min(uint2(q * float2(width, height)), uint2(width - 1, height - 1));
                if (z > gOpaqueDepth[qp]) hi = mid; else lo = mid;
            }
            float z;
            return projectGlassPoint(origin + direction * (0.5f * (lo + hi)), hitUV, z);
        }
        previous = difference;
        previousT = t;
    }
    return false;
}

float4 psGlass(VSOut vsOut, uint triangleIndex : SV_PrimitiveID, float3 barycentrics : SV_Barycentrics) : SV_TARGET
{
    BasicMaterialData material = gScene.materials.getBasicMaterialData(vsOut.materialID);
    float transmission = transmissionAmount(material);
    if (transmission <= 0.01f) discard;
    uint materialCount, stride;
    gLiquidMaterials.GetDimensions(materialCount, stride);
    bool isLiquid = vsOut.materialID < materialCount && gLiquidMaterials[vsOut.materialID] != 0;
    if (isLiquid != (gLiquidStage != 0)) discard;
    let lod = ExplicitGradientTextureSampler(ddx(vsOut.texC), ddy(vsOut.texC));
    if (alphaTest(vsOut, triangleIndex, lod)) discard;
    float3 V = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData sd = prepareShadingData(vsOut, triangleIndex, V);
    let mi = gScene.materials.getMaterialInstance(sd, lod, (uint)MaterialInstanceHints::None);
    let properties = mi.getProperties(sd);
    float3 N = properties.guideNormal;
    if (!all(isfinite(N)) || dot(N, N) < 1e-12f) N = vsOut.normalW;
    N = normalize(N);
    if (dot(N, V) < 0.f) N = -N;
    let header = gScene.materials.getMaterialHeader(vsOut.materialID);
    SamplerState glassSampler = gScene.materials.getTextureSampler(header.getDefaultTextureSamplerID());
    float4 spec = gScene.materials.sampleTexture(material.texSpecular, glassSampler, vsOut.texC, (float4)material.specular, lod);
    float roughness = spec.g;
#if MATERIAL_SYSTEM_HAS_SPEC_GLOSS_MATERIALS
    if (material.getShadingModel() == ShadingModel::SpecGloss) roughness = 1.f - spec.a;
#endif
    roughness = saturate(roughness + float(material.diffuseTransmission));
    float3 tint = saturate(gScene.materials.sampleTexture(material.texTransmission, glassSampler,
        vsOut.texC, float4((float3)material.transmission, 1.f), lod).rgb);
    // Wine/beer color is encoded as volume absorption in this scene, not as
    // baseColor or transmission color. Use a stable approximate path length.
    // This also restores absorption in the wine-bottle glass material.
    tint *= exp(-max((float3)material.volumeAbsorption, 0.f) * gLiquidOpticalDepth);
    float ior = max(float(header.getIoR()), 1.f);
    float f0 = (ior - 1.f) / (ior + 1.f);
    f0 *= f0;
    float fresnel = f0 + (1.f - f0) * pow(1.f - saturate(dot(N, V)), 5.f);

    uint width, height;
    gOpaqueColor.GetDimensions(width, height);
    float2 uv = vsOut.posH.xy / float2(width, height);
    float3 incident = -V;
    float3 transmitted = refract(incident, N, 1.f / ior);
    float thickness = max(0.001f * gSceneScale, 0.003f * length(vsOut.posW - gScene.camera.getPosition()));
    float2 refractedUV = uv, a, b;
    float za, zb;
    if (projectGlassPoint(vsOut.posW + transmitted * thickness, a, za) &&
        projectGlassPoint(vsOut.posW + incident * thickness, b, zb))
        refractedUV += gRefractionStrength * (a - b);
    // Keep foreground objects from being pulled through the glass silhouette.
    if (any(refractedUV <= 0.f) || any(refractedUV >= 1.f)) refractedUV = uv;
    uint2 pixel = min(uint2(refractedUV * float2(width, height)), uint2(width - 1, height - 1));
    if (gOpaqueDepth[pixel] < vsOut.posH.z) refractedUV = uv;
    float3 background = glassBackground(refractedUV, roughness);
    float3 reflectedDirection = reflect(incident, N);
    uint page; float2 lmUV;
    float3 incidentLight = float3(0.f);
    if (getBistroLightmapLookup(vsOut, triangleIndex, barycentrics, page, lmUV))
        incidentLight = sampleBistroAtlas(page, lmUV);
    // The baked diffuse incident light is a low-frequency fallback for indoor
    // regions absent from the screen-space reflection buffer.
    float3 reflection = incidentLight;
    if (gHasEnvironment != 0) reflection = gScene.envMap.eval(reflectedDirection, roughness * 6.f);
    float2 reflectedUV;
    if (traceGlassReflection(vsOut.posW, reflectedDirection, reflectedUV))
    {
        float edgeFade = saturate(15.f * min(min(reflectedUV.x, reflectedUV.y), min(1.f - reflectedUV.x, 1.f - reflectedUV.y)));
        reflection = lerp(
            reflection,
            glassBackground(reflectedUV, roughness),
            0.8f * edgeFade
        );
    }
    float3 baseColor = gScene.materials.sampleTexture(material.texBaseColor, glassSampler, vsOut.texC, (float4)material.baseColor, lod).rgb;
    float3 diffuse = baseColor * incidentLight;
    // Appearance controls for the background demo.
    float reflectionStrength = 0.45f;
    float transmissionBrightness = isLiquid ? 0.85f : 0.92f;

    float3 color =
    reflectionStrength * fresnel * reflection
    + (1.f - fresnel) * lerp(
        diffuse,
        transmissionBrightness * tint * background,
        transmission
    );
    //float3 emission = gScene.materials.sampleTexture(material.texEmissive, glassSampler,
    //    vsOut.texC, float4(material.emissive, 1.f), lod).rgb * material.emissiveFactor;
    float3 emission = gScene.materials.sampleTexture(material.texEmissive, glassSampler,
        vsOut.texC, float4(material.emissive, 1.f), lod).rgb;
    // Background is already included, so CPU blending stays disabled.
    return float4(max(color + emission, 0.f), 1.f);
}
