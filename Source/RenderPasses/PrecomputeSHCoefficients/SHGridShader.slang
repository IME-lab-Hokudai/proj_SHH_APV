import Scene.Raster; // Imports defaultVS.

SamplerState gLinearSampler;
#define MAX_SH_BASIS 9
// Indexing: probeIndex * numBasis + basisIndex
StructuredBuffer<float4> gSHCoeffs;

cbuffer gProbeGridInfo
{
    int3 resolution;
    int numBasis;
    float3 origin;
    float pad0;
    float3 spacing;
    float pad1;
};

struct MyVSOut : VSOut
{
    linear float4 irradiance : IRRADIANCE; 
};

MyVSOut vsMain(VSIn vIn)
{
    MyVSOut vsOut = (MyVSOut)defaultVS(vIn);
   // float4 shCoeffs[MAX_SH_BASIS];
    //ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    //sampleSHAtPosition(vsOut.posW, shCoeffs);
    //vsOut.irradiance = evaluateIrradiance(normalize(vsOut.normalW), shCoeffs);
    
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

// 3D to 1D index for the probe grid corresponding to cpu code
//int width = res.x;
//int height = res.y;
//int depth = res.z;
//for (int probeZ = 0; probeZ < depth; ++probeZ)
//    for (int probeY = 0; probeY < height; ++probeY)
 //       for (int probeX = 0; probeX < width; ++probeX)
int getProbeIndex(int3 coord)
{
    int width = resolution.x;
    int height = resolution.y;
    int depth = resolution.z;

    int3 clamped = clamp(coord, int3(0, 0, 0), resolution - 1);
    return clamped.x + width * clamped.y + width*height * clamped.z;
}

void sampleSHAtPosition(float3 worldPos, out float4 outSH[MAX_SH_BASIS])
{
    // Convert the world position to the probe grid coordinate space
    // by subtracting the grid origin and dividing by the probe spacing
    float3 gridPos = (worldPos - origin) / spacing;

    // Compute the integer base coordinate of the grid cell containing the position
    //note that this is the lower conner of the voxel (C00)
    int3 base = int3(floor(gridPos));

    // Compute the fractional part within the grid cell for interpolation weights
    //relative position of gridPos in the voxel
    float3 frac = gridPos - float3(base);

    // Initialize output SH coefficients to zero
    [unroll]
    for (int i = 0; i < MAX_SH_BASIS; ++i)
    {
        outSH[i] = float4(0, 0, 0, 0);
    }

    // Perform trilinear interpolation from the 8 surrounding probes
    // Loop over all corners of the cube surrounding the position
    [unroll]
    for (int z = 0; z <= 1; ++z)
    {
        [unroll]
        for (int y = 0; y <= 1; ++y)
        {
            [unroll]
            for (int x = 0; x <= 1; ++x)
            {
                // Calculate interpolation weights for each dimension
                // If coordinate is 0, weight is (1 - frac), else frac
                float3 w = float3(
                x == 0 ? (1.0 - frac.x) : frac.x,
                y == 0 ? (1.0 - frac.y) : frac.y,
                z == 0 ? (1.0 - frac.z) : frac.z
                );

                // Combine weights for 3D interpolation
                float weight = w.x * w.y * w.z;

                // Compute the neighbor probe coordinate
                int3 neighborCoord = base + int3(x, y, z);

                // Convert 3D coordinate to 1D probe index
                int neighborProbeIndex = getProbeIndex(neighborCoord);

                // Accumulate the weighted SH coefficients from the neighbor probe
                [loop]
                for (int i = 0; i < numBasis; ++i)
                {
                    outSH[i] += gSHCoeffs[neighborProbeIndex * numBasis + i] * weight;
                }
            }
        }
    }
}

// Evaluate 9 SH basis functions at direction n
float4 evaluateIrradiance(float3 n,  float4 shCoeffs[MAX_SH_BASIS])
{
    float x = n.x;
    float y = n.y;
    float z = n.z;

    // SH basis function indexing and formulas:
    // Index | l (band) | m (order) | SH function          | Formula (Cartesian coords)
    // -------------------------------------------------------------------------------------
    // 0     | 0        |  0        | Y_0^0                | 0.282095
    // 1     | 1        | -1        | Y_1^-1               | 0.488603 * y       // ~ sqrt(3/(4π)) * y/r * sin(φ)
    // 2     | 1        |  0        | Y_1^0                | 0.488603 * z       // ~ sqrt(3/(4π)) * z/r
    // 3     | 1        |  1        | Y_1^1                | 0.488603 * x       // ~ sqrt(3/(4π)) * x/r * cos(φ)
    // 4     | 2        | -2        | Y_2^-2               | 1.092548 * x * y   // ~ sqrt(15/(4π)) * xy/r² * sin(2φ)
    // 5     | 2        | -1        | Y_2^-1               | 1.092548 * y * z   // ~ sqrt(15/(4π)) * yz/r² * sin(φ)
    // 6     | 2        |  0        | Y_2^0                | 0.315392 * (-x*x - y*y + 2*z*z) // ~ 0.315392 * (3z² - 1)
    // 7     | 2        |  1        | Y_2^1                | 1.092548 * x * z   // ~ sqrt(15/(4π)) * xz/r² * cos(φ)
    // 8     | 2        |  2        | Y_2^2                | 0.546274 * (x*x - y*y)  // ~ 0.546274 * (x² - y²)/r² * cos(2φ)

    // Band 0
    float Y0 = 0.282095;

    // Band 1
    float Y1 = 0.488603 * y;
    float Y2 = 0.488603 * z;
    float Y3 = 0.488603 * x;

    // Band 2
    float Y4 = 1.092548 * x * y;
    float Y5 = 1.092548 * y * z;
    float Y6 = 0.315392 * (3.0 * z * z - 1.0);
    float Y7 = 1.092548 * x * z;
    float Y8 = 0.546274 * (x * x - y * y);
    float A[3] = { 3.141593, 2.094395, 0.785398 }; // Cosine lobe coefficients
    // Dot product: coeffs · basis
float4 result = 
        shCoeffs[0] * Y0 * A[0] +
        shCoeffs[1] * Y1 * A[1] +
        shCoeffs[2] * Y2 * A[1] +
        shCoeffs[3] * Y3 * A[1] +
        shCoeffs[4] * Y4 * A[2] +
        shCoeffs[5] * Y5 * A[2] +
        shCoeffs[6] * Y6 * A[2] +
        shCoeffs[7] * Y7 * A[2] +
        shCoeffs[8] * Y8 * A[2];

    return max(result, float4(0.f)); // Clamp to avoid negative lighting
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

  float4 shCoeffs[MAX_SH_BASIS];
  sampleSHAtPosition(shadingData.posW, shCoeffs);
  vsOut.irradiance = evaluateIrradiance(normalize(shadingData.faceN), shCoeffs);

  //float3 diffuseColor = baseColorSample.rgb * vsOut.irradiance.rgb / 3.14159265;
  float3 diffuseColor = baseColorSample.rgb * vsOut.irradiance.rgb;
  float3 emissiveColor = emissiveSample.rgb;
   ///float3 litColor = diffuseColor + emissiveColor;
   float3 litColor = diffuseColor;
  // float3 litColor = normalSampleW;
  // float3 litColor = shadingData.faceN;
  //float3 litColor = vsOut.irradiance.rgb;
    return float4(litColor,1.0f);
}
