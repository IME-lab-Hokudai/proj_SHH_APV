import Scene.Raster;
#include "AdaptiveGridData.slang"
#include "Utils/Math/MathConstants.slangh"

// Lightmaps
SamplerState gLinearSampler;
#define MAX_SH_BASIS 9
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

//Texture2D gFloorLightmap;
//Texture2D gLeftWallLightmap;
//Texture2D gRightWallLightmap;
//Texture2D gRoofLeftLightmap;
//Texture2D gRoofRightLightmap;

//Texture2D gPillar0Lightmap;
//Texture2D gPillar1Lightmap;
//Texture2D gPillar2Lightmap;
//Texture2D gPillar3Lightmap;
//Texture2D gPillar4Lightmap;
//Texture2D gPillar5Lightmap;
//Texture2D gPillar6Lightmap;
//Texture2D gPillar7Lightmap;

// Room shell lightmaps.
//Texture2D gFloorLightmap;
//Texture2D gCeilingLightmap;
//Texture2D gLeftWallLightmap;
//Texture2D gRightWallLightmap;
//Texture2D gBackWallLightmap;
//Texture2D gFrontWallLightmap;

// Box lightmaps.
//Texture2D gTallBoxALightmap;
//Texture2D gWideBoxBLightmap;
//Texture2D gBlockCLightmap;
//Texture2D gLowBoxDLightmap;
//Texture2D gThinSlabELightmap;
//Texture2D gThinSlabFLightmap;
//Texture2D gTallBoxJLightmap;
//Texture2D gShortBoxLLightmap;

//Cornelbox light map
Texture2D gFloorLightmap;
Texture2D gCeilingLightmap;
Texture2D gBackWallLightmap;
Texture2D gLeftWallLightmap;
Texture2D gRightWallLightmap;

// Slab lightmap.
Texture2D gVisibilitySlabLightmap;
cbuffer PerFrameCB
{
    //uint gFloorInstanceID;
    //uint gLeftWallInstanceID;
    //uint gRightWallInstanceID;
    //uint gRoofLeftInstanceID;
    //uint gRoofRightInstanceID;

    //uint gPillar0InstanceID;
    //uint gPillar1InstanceID;
    //uint gPillar2InstanceID;
    //uint gPillar3InstanceID;
    //uint gPillar4InstanceID;
    //uint gPillar5InstanceID;
    //uint gPillar6InstanceID;
    //uint gPillar7InstanceID;

    //float3 gPillar0CenterW;
    //float3 gPillar1CenterW;
    //float3 gPillar2CenterW;
    //float3 gPillar3CenterW;
    //float3 gPillar4CenterW;
    //float3 gPillar5CenterW;
    //float3 gPillar6CenterW;
    //float3 gPillar7CenterW;

    // float3 gPillar0HalfExtentW;
    // float3 gPillar1HalfExtentW;
    // float3 gPillar2HalfExtentW;
    // float3 gPillar3HalfExtentW;
    // float3 gPillar4HalfExtentW;
    // float3 gPillar5HalfExtentW;
    // float3 gPillar6HalfExtentW;
    // float3 gPillar7HalfExtentW;

    //datascene
    //uint gFloorInstanceID;
    //uint gCeilingInstanceID;
    //uint gLeftWallInstanceID;
    //uint gRightWallInstanceID;
    //uint gBackWallInstanceID;
    //uint gFrontWallInstanceID;

    //uint gTallBoxAInstanceID;
    //uint gWideBoxBInstanceID;
    //uint gBlockCInstanceID;
    //uint gLowBoxDInstanceID;
    //uint gThinSlabEInstanceID;
    //uint gThinSlabFInstanceID;
    //uint gTallBoxJInstanceID;
    //uint gShortBoxLInstanceID;

    //float3 gTallBoxACenterW;
    //float3 gWideBoxBCenterW;
    //float3 gBlockCCenterW;
    //float3 gLowBoxDCenterW;
    //float3 gThinSlabECenterW;
    //float3 gThinSlabFCenterW;
    //float3 gTallBoxJCenterW;
    //float3 gShortBoxLCenterW;

    //float3 gTallBoxAHalfExtentW;
    //float3 gWideBoxBHalfExtentW;
    //float3 gBlockCHalfExtentW;
    //float3 gLowBoxDHalfExtentW;
    //float3 gThinSlabEHalfExtentW;
    //float3 gThinSlabFHalfExtentW;
    //float3 gTallBoxJHalfExtentW;
    //float3 gShortBoxLHalfExtentW;

    //cornell box
    uint gFloorInstanceID;
    uint gCeilingInstanceID;
    uint gBackWallInstanceID;
    uint gLeftWallInstanceID;
    uint gRightWallInstanceID;
    uint gVisibilitySlabInstanceID;

    float3 gVisibilitySlabCenterW;
    float3 gVisibilitySlabHalfExtentW;
};

float3 rotateX(float3 v, float a)
{
    float s = sin(a);
    float c = cos(a);
    return float3(v.x, c * v.y - s * v.z, s * v.y + c * v.z);
}

float3 rotateY(float3 v, float a)
{
    float s = sin(a);
    float c = cos(a);
    return float3(c * v.x + s * v.z, v.y, -s * v.x + c * v.z);
}

float3 rotateZ(float3 v, float a)
{
    float s = sin(a);
    float c = cos(a);
    return float3(c * v.x - s * v.y, s * v.x + c * v.y, v.z);
}

float3 worldVectorToBoxLocal(float3 vW, float3 rotationEulerDeg)
{
    float3 r = radians(rotationEulerDeg);

    float3 v = vW;
    v = rotateZ(v, -r.z);
    v = rotateY(v, -r.y);
    v = rotateX(v, -r.x);
    return v;
}

float2 projectCubeFaceUV(
    float3 posW,
    float3 normalW,
    float3 centerW,
    float3 halfExtentW,
    float3 rotationEulerDeg
)
{
    float3 pLocal = worldVectorToBoxLocal(posW - centerW, rotationEulerDeg);
    float3 nLocal = normalize(worldVectorToBoxLocal(normalW, rotationEulerDeg));

    float3 p = pLocal / max(halfExtentW, float3(1e-6f));
    p = clamp(p, -1.0f, 1.0f);

    float3 n = normalize(nLocal);
    float3 an = abs(n);

    float2 uv;
    uint faceIndex = 0;

    if (an.x >= an.y && an.x >= an.z)
    {
        if (n.x > 0.0f)
        {
            uv = float2(-p.z, p.y);
            faceIndex = 0; // +X
        }
        else
        {
            uv = float2(p.z, p.y);
            faceIndex = 1; // -X
        }
    }
    else if (an.z >= an.x && an.z >= an.y)
    {
        if (n.z > 0.0f)
        {
            uv = float2(p.x, p.y);
            faceIndex = 2; // +Z
        }
        else
        {
            uv = float2(-p.x, p.y);
            faceIndex = 3; // -Z
        }
    }
    else
    {
        if (n.y > 0.0f)
        {
            uv = float2(p.x, -p.z);
            faceIndex = 4; // +Y
        }
        else
        {
            uv = float2(p.x, p.z);
            faceIndex = 5; // -Y
        }
    }

    uv = saturate(uv * 0.5f + 0.5f);

    uint tileX = faceIndex % 3;
    uint tileY = faceIndex / 3;

    float2 tileScale = float2(1.0f / 3.0f, 1.0f / 2.0f);
    float2 tileOffset = float2(tileX, tileY) * tileScale;

    return tileOffset + uv * tileScale;
}

float3 shadeWithLightmap(float3 albedo, float3 irradiance)
{
    return albedo * irradiance / 3.14159265f;
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

VSOut vsMain(VSIn vIn)
{
    VSOut vsOut = defaultVS(vIn);
    return vsOut;
}


float4 psMain(VSOut vsOut, uint triangleIndex : SV_PrimitiveID) : SV_TARGET
{
    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);

    float4 baseColorSample = (float4)baseColorMat.baseColor;
    float4 emissiveSample = float4(baseColorMat.emissive, 1.0f);

    float2 quadUV = vsOut.texC;
    quadUV.y = 1.0f - quadUV.y;

    if (baseColorMat.texBaseColor.getMode() == TextureHandle::Mode::Texture) {
        baseColorSample = gScene.materials.sampleTexture(baseColorMat.texBaseColor, gLinearSampler, vsOut.texC, 0.f);
    }

    if (baseColorMat.texEmissive.getMode() == TextureHandle::Mode::Texture) {
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);
    }

    float3 albedo = baseColorSample.rgb;
    float3 litColor = albedo;

    const uint instanceID = vsOut.instanceID.index;
    if (instanceID < 19) {
         //litColor = albedo;
         //discard;
    }

    // Data scene
    //if (instanceID == gFloorInstanceID)
    //{
    //    litColor = shadeWithLightmap(albedo, gFloorLightmap.Sample(gLinearSampler, quadUV).rgb);
    //}
    //else if (instanceID == gCeilingInstanceID)
    //{
    //    litColor = shadeWithLightmap(albedo, gCeilingLightmap.Sample(gLinearSampler, quadUV).rgb);
    //}
    //else if (instanceID == gLeftWallInstanceID)
    //{
    //    litColor = shadeWithLightmap(albedo, gLeftWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    //}
    //else if (instanceID == gRightWallInstanceID)
    //{
    //    litColor = shadeWithLightmap(albedo, gRightWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    //}
    //else if (instanceID == gBackWallInstanceID)
    //{
    //    litColor = shadeWithLightmap(albedo, gBackWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    //}
    //else if (instanceID == gFrontWallInstanceID)
    //{
    //    litColor = shadeWithLightmap(albedo, gFrontWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    //}

    //// Box targets.
    //else if (instanceID == gTallBoxAInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gTallBoxACenterW, gTallBoxAHalfExtentW, float3(0.0f, 0.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gTallBoxALightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gWideBoxBInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gWideBoxBCenterW, gWideBoxBHalfExtentW, float3(0.0f, 0.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gWideBoxBLightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gBlockCInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gBlockCCenterW, gBlockCHalfExtentW, float3(0.0f, 18.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gBlockCLightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gLowBoxDInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gLowBoxDCenterW, gLowBoxDHalfExtentW, float3(0.0f, -22.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gLowBoxDLightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gThinSlabEInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gThinSlabECenterW, gThinSlabEHalfExtentW, float3(0.0f, 0.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gThinSlabELightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gThinSlabFInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gThinSlabFCenterW, gThinSlabFHalfExtentW, float3(0.0f, 12.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gThinSlabFLightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gTallBoxJInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gTallBoxJCenterW, gTallBoxJHalfExtentW, float3(0.0f, -18.0f, 8.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gTallBoxJLightmap.Sample(gLinearSampler, uv).rgb);
    //}
    //else if (instanceID == gShortBoxLInstanceID)
    //{
    //    float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gShortBoxLCenterW, gShortBoxLHalfExtentW, float3(0.0f, 40.0f, 0.0f));
    //    uv.y = 1.0f - uv.y;
    //    litColor = shadeWithLightmap(albedo, gShortBoxLLightmap.Sample(gLinearSampler, uv).rgb);
    //}

    // Cornell box
    if (instanceID == gFloorInstanceID) {
        litColor = shadeWithLightmap(albedo, gFloorLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gCeilingInstanceID) {
        litColor = shadeWithLightmap(albedo, gCeilingLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gBackWallInstanceID) {
        litColor = shadeWithLightmap(albedo, gBackWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gLeftWallInstanceID) {
        litColor = shadeWithLightmap(albedo, gLeftWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gRightWallInstanceID) {
        litColor = shadeWithLightmap(albedo, gRightWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gVisibilitySlabInstanceID) {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gVisibilitySlabCenterW, gVisibilitySlabHalfExtentW, float3(0.0f, 0.0f, 0.0f));
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gVisibilitySlabLightmap.Sample(gLinearSampler, uv).rgb);
    }
    else {
        float3 finalSH[MAX_SH_BASIS];
        // 1. Find the Leaf Voxel
        int probeIdx = traverseOctree(shadingData.posW);
        // 2. Compute interpolated SH Coeffs
        fetchHermiteInterpolatedSH(probeIdx, shadingData.posW, finalSH);
        //fetchTrilinearInterpolatedSH(probeIdx, shadingData.posW, finalSH);

        float3 irradiance = evaluateIrradiance(shadingData.frame.N, finalSH);
        litColor = baseColorSample.rgb * irradiance.rgb / M_PI;
    }

    litColor += emissiveSample.rgb;
    return float4(litColor, 1.0f);
}
