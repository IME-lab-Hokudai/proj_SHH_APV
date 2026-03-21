import Scene.Raster;

SamplerState gLinearSampler;

Texture2D gFloorLightmap;
Texture2D gLeftWallLightmap;
Texture2D gRightWallLightmap;
Texture2D gRoofLeftLightmap;
Texture2D gRoofRightLightmap;

Texture2D gPillar0Lightmap;
Texture2D gPillar1Lightmap;
Texture2D gPillar2Lightmap;
Texture2D gPillar3Lightmap;
Texture2D gPillar4Lightmap;
Texture2D gPillar5Lightmap;
Texture2D gPillar6Lightmap;
Texture2D gPillar7Lightmap;

cbuffer PerFrameCB
{
    uint gFloorInstanceID;
    uint gLeftWallInstanceID;
    uint gRightWallInstanceID;
    uint gRoofLeftInstanceID;
    uint gRoofRightInstanceID;

    uint gPillar0InstanceID;
    uint gPillar1InstanceID;
    uint gPillar2InstanceID;
    uint gPillar3InstanceID;
    uint gPillar4InstanceID;
    uint gPillar5InstanceID;
    uint gPillar6InstanceID;
    uint gPillar7InstanceID;
}

Texture2D gLightmap; // The baked lightmap

VSOut vsMain(VSIn vIn)
{
    VSOut vsOut = defaultVS(vIn);
    
    return vsOut;
}

float4  psMain(VSOut vsOut, uint triangleIndex: SV_PrimitiveID) : SV_TARGET
{

    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);
    float4 baseColorSample = (float4)baseColorMat.baseColor;
    float4 emissiveSample = float4(baseColorMat.emissive, 1.0f);

    // FIX: Flip the V coordinate because when we store in exr file it is flipped
    float2 sampledUV = vsOut.texC;
    sampledUV.y = 1.0f - sampledUV.y;

    if (baseColorMat.texBaseColor.getMode() == TextureHandle::Mode::Texture) {
        baseColorSample = gScene.materials.sampleTexture(baseColorMat.texBaseColor, gLinearSampler, vsOut.texC, 0.f);
    }

    if (baseColorMat.texEmissive.getMode() == TextureHandle::Mode::Texture) {
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);
    }
    else {
        baseColorSample = (float4)baseColorMat.baseColor;
        emissiveSample = float4(baseColorMat.emissive, 1);
    }
    // We get the base color from the material system to multiply with our light
    float3 albedo = baseColorSample.rgb;
    float3 litColor = albedo;

    const uint instanceID = vsOut.instanceID.index;

    if (instanceID == gFloorInstanceID)
    {
        float3 irradiance = gFloorLightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gLeftWallInstanceID)
    {
        float3 irradiance = gLeftWallLightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gRightWallInstanceID)
    {
        float3 irradiance = gRightWallLightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gRoofLeftInstanceID)
    {
        float3 irradiance = gRoofLeftLightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gRoofRightInstanceID)
    {
        float3 irradiance = gRoofRightLightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar0InstanceID)
    {
        float3 irradiance = gPillar0Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar1InstanceID)
    {
        float3 irradiance = gPillar1Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar2InstanceID)
    {
        float3 irradiance = gPillar2Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar3InstanceID)
    {
        float3 irradiance = gPillar3Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar4InstanceID)
    {
        float3 irradiance = gPillar4Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar5InstanceID)
    {
        float3 irradiance = gPillar5Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar6InstanceID)
    {
        float3 irradiance = gPillar6Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar7InstanceID)
    {
        float3 irradiance = gPillar7Lightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }

    litColor += emissiveSample.rgb;
    return float4(litColor, 1.0f);
}
