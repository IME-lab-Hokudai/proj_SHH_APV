import Scene.Raster;

// Lightmaps
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

    float3 gPillar0CenterW;
    float3 gPillar1CenterW;
    float3 gPillar2CenterW;
    float3 gPillar3CenterW;
    float3 gPillar4CenterW;
    float3 gPillar5CenterW;
    float3 gPillar6CenterW;
    float3 gPillar7CenterW;

    float3 gPillar0HalfExtentW;
    float3 gPillar1HalfExtentW;
    float3 gPillar2HalfExtentW;
    float3 gPillar3HalfExtentW;
    float3 gPillar4HalfExtentW;
    float3 gPillar5HalfExtentW;
    float3 gPillar6HalfExtentW;
    float3 gPillar7HalfExtentW;
};

VSOut vsMain(VSIn vIn)
{
    VSOut vsOut = defaultVS(vIn);
    return vsOut;
}

float2 projectCubeFaceUV(float3 posW, float3 normalW, float3 centerW, float3 halfExtentW)
{
    float3 p = (posW - centerW) / max(halfExtentW, float3(1e-6f));
    p = clamp(p, -1.0f, 1.0f);

    float3 n = normalize(normalW);
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

    float2 tileScale  = float2(1.0f / 3.0f, 1.0f / 2.0f);
    float2 tileOffset = float2(tileX, tileY) * tileScale;

    return tileOffset + uv * tileScale;
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

    if (instanceID == gFloorInstanceID)
    {
        float3 irradiance = gFloorLightmap.Sample(gLinearSampler, quadUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gLeftWallInstanceID)
    {
        float3 irradiance = gLeftWallLightmap.Sample(gLinearSampler, quadUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gRightWallInstanceID)
    {
        float3 irradiance = gRightWallLightmap.Sample(gLinearSampler, quadUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gRoofLeftInstanceID)
    {
        float3 irradiance = gRoofLeftLightmap.Sample(gLinearSampler, quadUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gRoofRightInstanceID)
    {
        float3 irradiance = gRoofRightLightmap.Sample(gLinearSampler, quadUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar0InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar0CenterW, gPillar0HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar0Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar1InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar1CenterW, gPillar1HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar1Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar2InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar2CenterW, gPillar2HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar2Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar3InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar3CenterW, gPillar3HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar3Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar4InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar4CenterW, gPillar4HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar4Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar5InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar5CenterW, gPillar5HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar5Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar6InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar6CenterW, gPillar6HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar6Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }
    else if (instanceID == gPillar7InstanceID)
    {
        float2 pillarUV = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gPillar7CenterW, gPillar7HalfExtentW);
        pillarUV.y = 1.0f - pillarUV.y;
        float3 irradiance = gPillar7Lightmap.Sample(gLinearSampler, pillarUV).rgb;
        litColor = albedo * irradiance / 3.14159265f;
    }

    litColor += emissiveSample.rgb;
    return float4(litColor, 1.0f);
}
