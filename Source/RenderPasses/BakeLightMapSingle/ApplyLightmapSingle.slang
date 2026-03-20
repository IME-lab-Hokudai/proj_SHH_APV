import Scene.Raster;

cbuffer PerFrameCB
{
    uint gReceiverInstanceID;
}

SamplerState gLinearSampler;
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

    if (vsOut.instanceID.index == gReceiverInstanceID)
    {
        float3 irradiance = gLightmap.Sample(gLinearSampler, sampledUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }

    litColor += emissiveSample.rgb;
    return float4(litColor, 1.0f);
}
