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
// 1D Hermite Interpolation (RGB)
// t: Normalized distance (0..1)
// v0, v1: Values at endpoints
// s0, s1: Slopes at endpoints (Must be scaled by segment length!)
// ------------------------------------------------------------------
float3 hermite1D(float t, float3 v0, float3 s0, float3 v1, float3 s1)
{
    float t2 = t * t;
    float t3 = t2 * t;

    // Standard Cubic Hermite Basis functions
    float h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    float h10 = t3 - 2.0 * t2 + t;
    float h01 = -2.0 * t3 + 3.0 * t2;
    float h11 = t3 - t2;

    return h00 * v0 + h10 * s0 + h01 * v1 + h11 * s1;
}

// ------------------------------------------------------------------
// Standard Trilinear Interpolation (Baseline)
// Traversal Order: Z -> Y -> X (Matches C++ Bit Packing)
// ------------------------------------------------------------------
#if SH_GRID_UNIFORM
#include "UniformGridSampling.slang"
#else
void fetchTrilinearInterpolatedSH(int probeIdx, float3 pos, out float3 finalCoeffs[9])
{
    GPUProbe p = gProbeBuffer[probeIdx];
    float3 size = p.maxPoint - p.minPoint;

    // Normalized position within voxel (0..1)
    float3 t = saturate((pos - p.minPoint) / size);

    // Pre-load indices for the 8 corners
    int cIdx[8];
    [unroll]
    for (int k = 0; k < 8; ++k) cIdx[k] = p.corners[k];

    [unroll]
    for (int band = 0; band < 9; ++band)
    {
        // 1. Fetch Corner Values (No Gradients needed!)
        float3 v[8];
        for (int i = 0; i < 8; ++i)
        {
            v[i] = gCornerBuffer[cIdx[i]].coeffs[band].rgb;
        }

        // 2. Interpolate along Z (Stride 1, Bit 0)
        // Pairs (0,1), (2,3), (4,5), (6,7)
        float3 z0 = lerp(v[0], v[1], t.z);
        float3 z1 = lerp(v[2], v[3], t.z);
        float3 z2 = lerp(v[4], v[5], t.z);
        float3 z3 = lerp(v[6], v[7], t.z);

        // 3. Interpolate along Y (Stride 2, Bit 1)
        // Pairs (z0, z1) -> Low Y
        // Pairs (z2, z3) -> High Y
        float3 y0 = lerp(z0, z1, t.y);
        float3 y1 = lerp(z2, z3, t.y);

        // 4. Interpolate along X (Stride 4, Bit 2)
        finalCoeffs[band] = lerp(y0, y1, t.x);
    }
}

// ------------------------------------------------------------------
// Progressive Tricubic Hermite Interpolation
// ------------------------------------------------------------------
void fetchHermiteInterpolatedSH(int probeIdx, float3 pos, out float3 finalCoeffs[9])
{
    GPUProbe p = gProbeBuffer[probeIdx];
    float3 size = p.maxPoint - p.minPoint;
    float3 t = saturate((pos - p.minPoint) / size);
    // float3 t = (pos - p.minPoint) / size;

    int cIdx[8];
    for (int k = 0; k < 8; ++k) cIdx[k] = p.corners[k];

    for (int band = 0; band < 9; ++band)
    {
        float3 v[8];
        // Gradients mapped to Falcor Axes
        float3 g_FalcorX[8];
        float3 g_FalcorY[8];
        float3 g_FalcorZ[8];

        for (int i = 0; i < 8; ++i)
        {
            GPUCorner c = gCornerBuffer[cIdx[i]];
            v[i] = c.coeffs[band].rgb;

            // ----------------------------------------------------------------
            // FIXED MAPPING:
            // Falcor X = SH Y  --> Use .y (d/dy_sh)
            // Falcor Y = SH Z  --> Use .z (d/dz_sh)
            // Falcor Z = SH X  --> Use .x (d/dx_sh)
            // ----------------------------------------------------------------
            g_FalcorX[i] = float3(c.gradR[band].y, c.gradG[band].y, c.gradB[band].y);
            g_FalcorY[i] = float3(c.gradR[band].z, c.gradG[band].z, c.gradB[band].z);
            g_FalcorZ[i] = float3(c.gradR[band].x, c.gradG[band].x, c.gradB[band].x);
        }

        // ----------------------------------------------------------------
        // PASS 1: Interpolate along Falcor Z (Stride 1)
        // C++: i and i+1 differ by size.z
        // ----------------------------------------------------------------
        float3 q[4];
        float3 q_dX[4];
        float3 q_dY[4];

        for (int i = 0; i < 4; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            // Hermite using Z-Gradient (g_FalcorZ)
            // Note: Multiply by size.z to scale world-derivative to 0..1 space
            q[i] = hermite1D(t.z, v[i0], g_FalcorZ[i0] * size.z, v[i1], g_FalcorZ[i1] * size.z);

            // Linearly interpolate the other gradients
            q_dX[i] = lerp(g_FalcorX[i0], g_FalcorX[i1], t.z);
            q_dY[i] = lerp(g_FalcorY[i0], g_FalcorY[i1], t.z);
        }

        // ----------------------------------------------------------------
        // PASS 2: Interpolate along Falcor Y (Stride 2)
        // C++: i and i+2 differ by size.y
        // ----------------------------------------------------------------
        float3 r[2];
        float3 r_dX[2];

        for (int i = 0; i < 2; ++i)
        {
            int i0 = 2 * i;
            int i1 = 2 * i + 1;

            // Hermite using Y-Gradient (g_FalcorY)
            r[i] = hermite1D(t.y, q[i0], q_dY[i0] * size.y, q[i1], q_dY[i1] * size.y);

            r_dX[i] = lerp(q_dX[i0], q_dX[i1], t.y);
        }

        // ----------------------------------------------------------------
        // PASS 3: Interpolate along Falcor X (Stride 4)
        // C++: i and i+4 differ by size.x
        // ----------------------------------------------------------------
        // Hermite using X-Gradient (g_FalcorX)
        finalCoeffs[band] = hermite1D(t.x, r[0], r_dX[0] * size.x, r[1], r_dX[1] * size.x);
    }
}

int findSeedProbe(float3 pos)
{
    if (gUseSeedGrid == 0)
    {
        return 0;
    }

    float3 rel = (pos - gSeedMinPoint) / gSeedCellSize;
    int3 cell = int3(floor(rel));

    if (cell.x < 0 || cell.y < 0 || cell.z < 0 ||
        cell.x >= int(gSeedResolution.x) ||
        cell.y >= int(gSeedResolution.y) ||
        cell.z >= int(gSeedResolution.z))
    {
        return -1;
    }

    uint idx = (uint(cell.z) * gSeedResolution.y + uint(cell.y)) * gSeedResolution.x + uint(cell.x);
    return gSeedProbeIndices[idx];
}

int traverseOctree(float3 pos)
{
    int currentIdx = findSeedProbe(pos);
    if (currentIdx < 0)
    {
        return -1;
    }

    for (int i = 0; i < OCTREE_MAX_DEPTH; ++i)
    {
        GPUProbe p = gProbeBuffer[currentIdx];

        // Extra safety: position must be inside this probe
        if (pos.x < p.minPoint.x || pos.y < p.minPoint.y || pos.z < p.minPoint.z ||
            pos.x > p.maxPoint.x || pos.y > p.maxPoint.y || pos.z > p.maxPoint.z)
        {
            return -1;
        }

        if (p.children[0] == -1)
        {
            return currentIdx;
        }

        float3 center = 0.5 * (p.minPoint + p.maxPoint);
        int octant = 0;

        if (pos.x >= center.x) octant |= 4;
        if (pos.y >= center.y) octant |= 2;
        if (pos.z >= center.z) octant |= 1;

        int childIdx = p.children[octant];
        if (childIdx == -1)
        {
            return currentIdx;
        }

        currentIdx = childIdx;
    }

    return currentIdx;
}

bool sampleGridSH(float3 pos, out float3 coeffs[MAX_SH_BASIS])
{
    [unroll]
    for (int band = 0; band < MAX_SH_BASIS; ++band) coeffs[band] = float3(0.f);
    int probeIdx = traverseOctree(pos);
    if (probeIdx < 0) return false;
    fetchHermiteInterpolatedSH(probeIdx, pos, coeffs);
    //fetchTrilinearInterpolatedSH(pos, coeffs);
    return true;
}
#endif

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
