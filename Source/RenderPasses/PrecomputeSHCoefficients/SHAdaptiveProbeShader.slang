import Scene.Raster; // Imports defaultVS.
#include "ProbeSamplingData.slang"
#include "Utils/Math/MathConstants.slangh"

SamplerState gLinearSampler;
#define MAX_SH_BASIS 9
#define OCTREE_MAX_DEPTH 9
StructuredBuffer<GPUCorner> gCornerBuffer;
StructuredBuffer<GPUProbe> gProbeBuffer;

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
void fetchHermiteInterpolatedSH(int probeIdx, float3 pos, out float3 finalCoeffs[9])
{
    GPUProbe p = gProbeBuffer[probeIdx];
    float3 size = p.maxPoint - p.minPoint;
    //float3 t = saturate((pos - p.minPoint) / size);
    float3 t = (pos - p.minPoint) / size;

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
    for (int k = 0; k < 8; ++k) cIdx[k] = p.corners[k];

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
// Taylor-Weighted Interpolation (Annen et al. 2004 / Ward & Heckbert)
// ------------------------------------------------------------------
void fetchTaylorInterpolatedSH(int probeIdx, float3 pos, out float3 finalCoeffs[9])
{
    GPUProbe p = gProbeBuffer[probeIdx];
    float3 size = p.maxPoint - p.minPoint;
    float3 t = saturate((pos - p.minPoint) / size); // Trilinear Weights

    int cIdx[8];
    for (int k = 0; k < 8; ++k) cIdx[k] = p.corners[k];

    for (int band = 0; band < 9; ++band)
    {
        float3 accum = 0.0;

        // Loop through all 8 corners
        // We calculate 8 separate predictions and blend them trilinearly
        for (int i = 0; i < 8; ++i)
        {
            // 1. Get Corner Data
            GPUCorner c = gCornerBuffer[cIdx[i]];
            float3 val = c.coeffs[band].rgb;

            // 2. Get Gradient (Using your Z-Up mapping)
            // Falcor X = SH Y (.y), Falcor Y = SH Z (.z), Falcor Z = SH X (.x)
            float3 grad;
            grad.x = c.gradR[band].y; // d/dx (Falcor)
            grad.y = c.gradR[band].z; // d/dy (Falcor)
            grad.z = c.gradR[band].x; // d/dz (Falcor)
            // Note: Use gradG/gradB for Green/Blue channels in production code

            // 3. Taylor Extrapolation: Pred = L + Gradient * Displacement
            // We need the vector from THIS corner to the pixel 'pos'
            float3 cornerPos = p.minPoint + size * float3(
                (i & 4) ? 1.0 : 0.0, // X bit
                                                          (i & 2) ? 1.0 : 0.0, // Y bit
                                                          (i & 1) ? 1.0 : 0.0  // Z bit
                                                   );
            float3 disp = pos - cornerPos;

            float3 prediction = val + grad * disp; // "First-order Taylor polynomial"

            // 4. Calculate Trilinear Weight for this corner
            // If corner is at 0, weight is (1-t). If at 1, weight is t.
            float wx = (i & 4) ? t.x : (1.0 - t.x);
            float wy = (i & 2) ? t.y : (1.0 - t.y);
            float wz = (i & 1) ? t.z : (1.0 - t.z);
            float weight = wx * wy * wz;

            accum += prediction * weight;
        }
        finalCoeffs[band] = accum;
    }
}

// ------------------------------------------------------------------
// Octree Traversal Logic 
// ------------------------------------------------------------------
// Returns the index of the LEAF probe containing 'pos'
int traverseOctree(float3 pos)
{
    int currentIdx = 0; // Start at Root

    // Safety loop limit (e.g. max depth 16)
    for (int i = 0; i < OCTREE_MAX_DEPTH; ++i)
    {
        GPUProbe p = gProbeBuffer[currentIdx];

        // Check if Leaf: children[0] is -1 (based on C++ init)
        if (p.children[0] == -1)
        {
            return currentIdx;
        }

        // Determine Octant
        // Matches C++: (k & 4)?X, (k & 2)?Y, (k & 1)?Z
        float3 center = (p.minPoint + p.maxPoint) * 0.5;
        int octant = 0;

        if (pos.x >= center.x) octant |= 4; // Bit 2 = X
        if (pos.y >= center.y) octant |= 2; // Bit 1 = Y
        if (pos.z >= center.z) octant |= 1; // Bit 0 = Z

        int childIdx = p.children[octant];

        // If tree is incomplete or we hit a boundary, stay here
        if (childIdx == -1) return currentIdx;

        // Descend
        currentIdx = childIdx;
    }
    return currentIdx;
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
    float4 baseColorSample = float4(0, 0, 0, 1);
    float4 emissiveSample = float4(0, 0, 0, 1);
    float4 specularSample = float4(0, 0, 0, 1);
    float3 normalSample;
    float3 normalSampleW;
   if (baseColorMat.texBaseColor.getMode() == TextureHandle::Mode::Texture) {
        baseColorSample = gScene.materials.sampleTexture(baseColorMat.texBaseColor, gLinearSampler, vsOut.texC, 0.f);
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);
        specularSample = gScene.materials.sampleTexture(baseColorMat.texSpecular, gLinearSampler, vsOut.texC, 0.f);
        normalSample = gScene.materials.sampleTexture(baseColorMat.texNormalMap, gLinearSampler, vsOut.texC, 0.f).rgb;
        normalSampleW = NormalSampleToWorldSpace(normalSample, normalize(vsOut.normalW), vsOut.tangentW.rgb);
   }
  else {
        baseColorSample = (float4)baseColorMat.baseColor;
        emissiveSample = float4(baseColorMat.emissive,1);
        specularSample = baseColorMat.specular;
        normalSample = normalize(vsOut.normalW);
        normalSampleW = normalSample;
  }

  float3 finalSH[MAX_SH_BASIS];
  // 1. Find the Leaf Voxel
   int probeIdx = traverseOctree(shadingData.posW);
   // 2. Compute interpolated SH Coeffs
   fetchHermiteInterpolatedSH(probeIdx, shadingData.posW, finalSH);
   // fetchTaylorInterpolatedSH(probeIdx, shadingData.posW, finalSH);
   // fetchTrilinearInterpolatedSH(probeIdx, shadingData.posW, finalSH);

  // ===============================================================================================
  // STEP 1: INITIAL SETUP & SELF-SAMPLING
  // ===============================================================================================
  // First, we determine which voxel (Probe) the pixel currently resides in.
  //int myIdx = traverseOctree(vsOut.posW);
  //GPUProbe myProbe = gProbeBuffer[myIdx];

  //// Fetch the High-Quality SH data from *this* voxel.
  //// This is our "Base Layer." If we are in the center of the voxel, we use this 100%.
  ////float3 finalSH[9];
  //fetchHermiteInterpolatedSH(myIdx, vsOut.posW, finalSH);

  //// ===============================================================================================
  //// STEP 2: DEFINE BLEND GEOMETRY
  //// ===============================================================================================
  //// We need to know how close this pixel is to the walls of the voxel.
  //float3 myCenter = (myProbe.minPoint + myProbe.maxPoint) * 0.5;
  //float3 mySize = myProbe.maxPoint - myProbe.minPoint;

  //// 'offset' is the vector from center to pixel.
  //float3 offset = vsOut.posW - myCenter;

  //// 'distToWall' tells us how far we are from the boundary (Positive value = Inside).
  //// e.g., if size is 10 and we are at x=4, distToWall is 1.0.
  //float3 distToWall = (mySize * 0.5) - abs(offset);

  //// Define the "Safe Zone" vs. "Blend Zone".
  //// We use 40% (0.4) of the voxel size.
  //// - Inner 20%: Pure High-Res Data.
  //// - Outer 40%: Gradually blend to neighbor.
  //float blendZone = 0.4 * min(min(mySize.x, mySize.y), mySize.z);

  //// ===============================================================================================
  //// STEP 3: INDEPENDENT AXIS BLENDING
  //// ===============================================================================================
  //// Why independent?
  //// If we are in a CORNER, we are close to 2 or 3 walls at once (e.g., Left and Top).
  //// We must blend with *both* the Left Neighbor AND the Top Neighbor to avoid diagonal seams.

  //// -----------------------------------------------------------------------
  //// BLEND PASS A: X-AXIS (Left/Right Neighbors)
  //// -----------------------------------------------------------------------
  //if (distToWall.x < blendZone)
  //{
  //    // 1. Identify which face we are near (-X or +X)
  //    // faceIdx 0 = Left (-X), faceIdx 1 = Right (+X)
  //    int faceIdx = (offset.x < 0) ? 0 : 1;

  //    // 2. Look up the Pre-Computed Neighbor Index
  //    // This array was built in C++. If the neighbor is smaller/same size, this is -1.
  //    // We ONLY blend if the neighbor is Coarser (Larger).
  //    int neighborIdx = myProbe.coarseNeighbors[faceIdx];

  //    if (neighborIdx != -1)
  //    {
  //        // 3. Fetch the Neighbor's Low-Res SH Data
  //        float3 neighborSH[9];
  //        fetchHermiteInterpolatedSH(neighborIdx, vsOut.posW, neighborSH);

  //        // 4. Calculate Blend Weight
  //        // Map distance (BlendZone -> 0) to Weight (0.0 -> 1.0)
  //        // - At blendZone boundary: weight is 0.0 (Pure Self)
  //        // - At the wall: weight is 1.0 (Pure Neighbor)
  //        float t = 1.0 - (distToWall.x / blendZone);
  //        t = smoothstep(0.0, 1.0, t); // Use smooth curve to hide the transition start

  //        // 5. Accumulate
  //        // We mix the neighbor into our current result.
  //        for (int i = 0; i < 9; ++i) finalSH[i] = lerp(finalSH[i], neighborSH[i], t);
  //    }
  //}

  //// -----------------------------------------------------------------------
  //// BLEND PASS B: Y-AXIS (Down/Up Neighbors)
  //// -----------------------------------------------------------------------
  //if (distToWall.y < blendZone)
  //{
  //    // faceIdx 2 = Down (-Y), faceIdx 3 = Up (+Y)
  //    int faceIdx = (offset.y < 0) ? 2 : 3;
  //    int neighborIdx = myProbe.coarseNeighbors[faceIdx];

  //    if (neighborIdx != -1)
  //    {
  //        float3 neighborSH[9];
  //        fetchHermiteInterpolatedSH(neighborIdx, vsOut.posW, neighborSH);

  //        float t = smoothstep(0.0, 1.0, 1.0 - (distToWall.y / blendZone));
  //        for (int i = 0; i < 9; ++i) finalSH[i] = lerp(finalSH[i], neighborSH[i], t);
  //    }
  //}

  //// -----------------------------------------------------------------------
  //// BLEND PASS C: Z-AXIS (Back/Front Neighbors)
  //// -----------------------------------------------------------------------
  //if (distToWall.z < blendZone)
  //{
  //    // faceIdx 4 = Back (-Z), faceIdx 5 = Front (+Z)
  //    int faceIdx = (offset.z < 0) ? 4 : 5;
  //    int neighborIdx = myProbe.coarseNeighbors[faceIdx];

  //    if (neighborIdx != -1)
  //    {
  //        float3 neighborSH[9];
  //        fetchHermiteInterpolatedSH(neighborIdx, vsOut.posW, neighborSH);

  //        float t = smoothstep(0.0, 1.0, 1.0 - (distToWall.z / blendZone));
  //        for (int i = 0; i < 9; ++i) finalSH[i] = lerp(finalSH[i], neighborSH[i], t);
  //    }
  //}

  float3 irradiance = evaluateIrradiance(normalize(shadingData.frame.N), finalSH);

  float3 diffuseColor = baseColorSample.rgb * irradiance.rgb / M_PI;
  //float3 diffuseColor = irradiance.rgb;
   float3 emissiveColor = emissiveSample.rgb;
   //float3 litColor = diffuseColor + emissiveColor;
   float3 litColor = diffuseColor;
    return float4(litColor,1.0f);
}
