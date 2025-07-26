import Scene.Raster; // Imports defaultVS.

SamplerState gLinearSampler;

cbuffer PerFrameCB
{
    float4 shCoeffs[9];
};

struct MyVSOut : VSOut
{
    linear float4 irradiance : IRRADIANCE; // Choose any free TEXCOORD semantic
};

MyVSOut vsMain(VSIn vIn)
{
    MyVSOut vsOut = (MyVSOut)defaultVS(vIn);
    vsOut.irradiance = evaluateIrradiance(vsOut.normalW);

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
float4 evaluateIrradiance(float3 n)
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

    // Dot product: coeffs · basis
    float4 result = 
        shCoeffs[0] * Y0 +
        shCoeffs[1] * Y1 +
        shCoeffs[2] * Y2 +
        shCoeffs[3] * Y3 +
        shCoeffs[4] * Y4 +
        shCoeffs[5] * Y5 +
        shCoeffs[6] * Y6 +
        shCoeffs[7] * Y7 +
        shCoeffs[8] * Y8;

    return max(result, float4(0.0)); // Clamp to avoid negative lighting
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

    LightData light = gScene.getLight(0);

    // The light vector aims opposite the direction the light rays travel.
    float3 lightPos = mul(light.transMat, float4(0, 0, 0, 1.f)).xyz;
    // float3 toLightDir = normalize(lightPos - vsOut.posW); //this is for area light source
    float3 toLightDir = -light.dirW; //directional light source

    // ambient
    float3 ambientColor = float3(0.01f, 0.01f, 0.01f) * light.intensity;

    // diffuse calculation
    float ndotl = max(dot(toLightDir, normalSampleW), 0.0f);
    //float3 diffuseColor = ndotl * baseColorSample.rgb * light.intensity;
   float3 diffuseColor = baseColorSample.rgb * vsOut.irradiance.rgb;

    // specular calculation
    const float shininess = specularSample.r*256.0f;
    float3 halfVec = normalize(viewDir + toLightDir);
    float ndoth = max(dot(halfVec, normalSampleW), 0.0f);
    float3 specular = pow(ndoth, shininess) * specularSample.rgb * light.intensity;
    // Our spec formula goes outside [0,1] range, but we are
    // doing LDR rendering.  So scale it down a bit.
    specular = specular / (specular + 1.0f);

    // float3 litColor = (diffuseColor + specColor) * light.intensity  + emissiveSample.rgb;
    // float3 litColor =  diffuseColor   + emissiveSample.rgb;
    float3 litColor = diffuseColor;
    return float4(litColor,1.0f);
    // return emissiveSample;
    //return baseColorSample;
}
