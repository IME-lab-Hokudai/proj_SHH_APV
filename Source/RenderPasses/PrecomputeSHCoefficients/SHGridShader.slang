import Scene.Raster; // Imports defaultVS.
#include "ProbeSamplingData.slang"
#include "Utils/Math/MathConstants.slangh"

SamplerState gLinearSampler;
#define MAX_SH_BASIS 9

StructuredBuffer<UniformGridCorner> gUniformProbes;
cbuffer PerFrameCB
{
    float3 gGridMin;
    float pad0;
    float3 gCellSize;
    float pad1;
    uint3 gProbeCountDim; // Dimensions of VERTICES (e.g., 17, 17, 17)
    uint pad2;
};

// ------------------------------------------------------------------
// Helper: Index Calculation
// ------------------------------------------------------------------
// Converts 3D integer coordinates to flattened 1D buffer index
uint getUniformGridIndex(int3 coord)
{
    // Layout: X changes fastest, then Y, then Z.
    // Index = x + (y * Width) + (z * Width * Height)
    uint strideY = gProbeCountDim.x;
    uint strideZ = gProbeCountDim.x * gProbeCountDim.y;

    return coord.x + (coord.y * strideY) + (coord.z * strideZ);
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
// Progressive Tricubic Hermite Interpolation
// ------------------------------------------------------------------
void fetchHermiteInterpolatedSH(float3 pos, out float3 finalCoeffs[9])
{
    // 1. Calculate Grid Coordinates (Cell Logic)
    // ----------------------------------------------------------
    float3 relPos = pos - gGridMin;

    // Determine which "Box" (Cell) we are in
    int3 cellCoords = int3(floor(relPos / gCellSize));

    // Clamp to valid range (0 to N-2) so we can always look up neighbor +1
    int3 maxCellIndex = int3(gProbeCountDim) - int3(1, 1, 1);
    cellCoords = max(int3(0), min(cellCoords, maxCellIndex - int3(1)));

    // 2. Calculate Local 't' (0..1)
    // ----------------------------------------------------------
    float3 cellOrigin = gGridMin + float3(cellCoords) * gCellSize;
    //float3 t = (pos - cellOrigin) / gCellSize;
    float3 t = saturate((pos - cellOrigin) / gCellSize);
    // 3. Fetch Indices (Using Helper)
    // ----------------------------------------------------------
    // We order these to match the loop bit-logic: Bit 0=Z, Bit 1=Y, Bit 2=X
    int cIdx[8];
    cIdx[0] = getUniformGridIndex(cellCoords + int3(0, 0, 0)); // 000
    cIdx[1] = getUniformGridIndex(cellCoords + int3(0, 0, 1)); // 001
    cIdx[2] = getUniformGridIndex(cellCoords + int3(0, 1, 0)); // 010
    cIdx[3] = getUniformGridIndex(cellCoords + int3(0, 1, 1)); // 011
    cIdx[4] = getUniformGridIndex(cellCoords + int3(1, 0, 0)); // 100
    cIdx[5] = getUniformGridIndex(cellCoords + int3(1, 0, 1)); // 101
    cIdx[6] = getUniformGridIndex(cellCoords + int3(1, 1, 0)); // 110
    cIdx[7] = getUniformGridIndex(cellCoords + int3(1, 1, 1)); // 111
    float3 size = gCellSize;
    for (int band = 0; band < 9; ++band)
    {
        float3 v[8];
        // Gradients mapped to Falcor Axes
        float3 g_FalcorX[8];
        float3 g_FalcorY[8];
        float3 g_FalcorZ[8];

        for (int i = 0; i < 8; ++i)
        {
            UniformGridCorner c = gUniformProbes[cIdx[i]];
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

struct MyVSOut : VSOut
{
    linear float4 irradiance : IRRADIANCE; 
};

MyVSOut vsMain(VSIn vIn)
{
    MyVSOut vsOut = (MyVSOut)defaultVS(vIn);
    return vsOut;
}

//---------------------------------------------------------------------------------------
// Transforms a normal map sample to world space.
//---------------------------------------------------------------------------------------
float3 NormalSampleToWorldSpace(float3 normalMapSample, float3 unitNormalW, float3 tangentW)
{
    // Uncompress each component from [0,1] to [-1,1].
    float3 normalT = 2.0f * normalMapSample - 1.0f;

    // Build orthonormal basis.
    float3 N = unitNormalW;
    float3 T = normalize(tangentW - dot(tangentW, N) * N);
    float3 B = cross(N, T);

    float3x3 TBN = float3x3(T, B, N);

    // Transform from tangent space to world space.
    float3 bumpedNormalW = mul(normalT, TBN);

    return bumpedNormalW;
}

// Evaluate 9 SH basis functions at direction n
float3 evaluateIrradiance(float3 n,  float3 shCoeffs[MAX_SH_BASIS])
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

    //with Condon–Shortley phase(-1) ^ m to match Iwasaki sensei code convention
    // Band 0
    float Y0 = 0.2820947917738781f; // l=0, m=0

     //Band 1
    float Y1 = -0.4886025119029200f * y;  // l=1, m=-1
    float Y2 = 0.4886025119029200f * z;   // l=1, m=0
    float Y3 = -0.4886025119029200f * x; // l=1, m=1

    // Band 2
    float Y4 = 1.0925484305920792f * x * y;       // l=2, m=-2
    float Y5 = -1.0925484305920792f * y * z;      // l=2, m=-1
    float Y6 = 0.9461746957575601f * z * z - 0.3153915652525200f; // l=2, m=0
    float Y7 = -1.0925484305920792f * x * z;                 // l=2, m=1
    float Y8 = 0.546274f * (x * x - y * y);       // l=2, m=2

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
    return max(result, float3(0.f));
}

float4 psMain(MyVSOut vsOut, uint triangleIndex: SV_PrimitiveID) : SV_TARGET
{
    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);
    float4 baseColorSample = (float4)baseColorMat.baseColor;
    float4 emissiveSample = float4(baseColorMat.emissive, 1.0f);
    float4 specularSample = baseColorMat.specular;
    float3 normalSample = normalize(vsOut.normalW);
    float3 normalSampleW = normalSample;
    // 2. Sample textures ONLY if they actually exist in the slot
    if (baseColorMat.texBaseColor.getMode() == TextureHandle::Mode::Texture) {
        baseColorSample = gScene.materials.sampleTexture(baseColorMat.texBaseColor, gLinearSampler, vsOut.texC, 0.f);
    }

    if (baseColorMat.texEmissive.getMode() == TextureHandle::Mode::Texture) {
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);
    }

    if (baseColorMat.texSpecular.getMode() == TextureHandle::Mode::Texture) {
        specularSample = gScene.materials.sampleTexture(baseColorMat.texSpecular, gLinearSampler, vsOut.texC, 0.f);
    }

    if (baseColorMat.texNormalMap.getMode() == TextureHandle::Mode::Texture) {
        float3 texNormal = gScene.materials.sampleTexture(baseColorMat.texNormalMap, gLinearSampler, vsOut.texC, 0.f).rgb;
        normalSampleW = NormalSampleToWorldSpace(texNormal, normalize(vsOut.normalW), vsOut.tangentW.rgb);
    }

  float3 finalSH[MAX_SH_BASIS];

  fetchHermiteInterpolatedSH(vsOut.posW, finalSH);

  float3 irradiance = evaluateIrradiance(normalize(shadingData.frame.N), finalSH);

  float3 diffuseColor = baseColorSample.rgb * irradiance.rgb / M_PI;
  //float3 diffuseColor = irradiance.rgb;
   float3 emissiveColor = emissiveSample.rgb;
   float3 litColor = diffuseColor + emissiveColor;
   //float3 litColor = diffuseColor;
    return float4(litColor,1.0f);
}
