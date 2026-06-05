import Scene.Raster;

SamplerState gLinearSampler;

// Room shell lightmaps.
Texture2D gFloorLightmap;
Texture2D gCeilingLightmap;
Texture2D gLeftWallLightmap;
Texture2D gRightWallLightmap;
Texture2D gBackWallLightmap;
Texture2D gFrontWallLightmap;

// Box lightmaps.
Texture2D gTallBoxALightmap;
Texture2D gWideBoxBLightmap;
Texture2D gBlockCLightmap;
Texture2D gLowBoxDLightmap;
Texture2D gThinSlabELightmap;
Texture2D gThinSlabFLightmap;
Texture2D gTallBoxJLightmap;
Texture2D gShortBoxLLightmap;

cbuffer PerFrameCB
{
    uint gFloorInstanceID;
    uint gCeilingInstanceID;
    uint gLeftWallInstanceID;
    uint gRightWallInstanceID;
    uint gBackWallInstanceID;
    uint gFrontWallInstanceID;

    uint gTallBoxAInstanceID;
    uint gWideBoxBInstanceID;
    uint gBlockCInstanceID;
    uint gLowBoxDInstanceID;
    uint gThinSlabEInstanceID;
    uint gThinSlabFInstanceID;
    uint gTallBoxJInstanceID;
    uint gShortBoxLInstanceID;

    float3 gTallBoxACenterW;
    float3 gWideBoxBCenterW;
    float3 gBlockCCenterW;
    float3 gLowBoxDCenterW;
    float3 gThinSlabECenterW;
    float3 gThinSlabFCenterW;
    float3 gTallBoxJCenterW;
    float3 gShortBoxLCenterW;

    float3 gTallBoxAHalfExtentW;
    float3 gWideBoxBHalfExtentW;
    float3 gBlockCHalfExtentW;
    float3 gLowBoxDHalfExtentW;
    float3 gThinSlabEHalfExtentW;
    float3 gThinSlabFHalfExtentW;
    float3 gTallBoxJHalfExtentW;
    float3 gShortBoxLHalfExtentW;
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

    float2 tileScale = float2(1.0f / 3.0f, 1.0f / 2.0f);
    float2 tileOffset = float2(tileX, tileY) * tileScale;

    return tileOffset + uv * tileScale;
}

float3 shadeWithLightmap(float3 albedo, float3 irradiance)
{
    return albedo * irradiance / 3.14159265f;
}

float4 psMain(VSOut vsOut, uint triangleIndex : SV_PrimitiveID) : SV_TARGET
{
    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);

    float4 baseColorSample = (float4)baseColorMat.baseColor;
    float4 emissiveSample = float4(baseColorMat.emissive, 1.0f);

    if (baseColorMat.texBaseColor.getMode() == TextureHandle::Mode::Texture)
    {
        baseColorSample = gScene.materials.sampleTexture(baseColorMat.texBaseColor, gLinearSampler, vsOut.texC, 0.f);
    }

    if (baseColorMat.texEmissive.getMode() == TextureHandle::Mode::Texture)
    {
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);
    }

    float3 albedo = baseColorSample.rgb;
    float3 litColor = albedo;

    float2 quadUV = vsOut.texC;
    quadUV.y = 1.0f - quadUV.y;

    uint instanceID = vsOut.instanceID.index;

    // Room shell.
    if (instanceID == gFloorInstanceID)
    {
        litColor = shadeWithLightmap(albedo, gFloorLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gCeilingInstanceID)
    {
        litColor = shadeWithLightmap(albedo, gCeilingLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gLeftWallInstanceID)
    {
        litColor = shadeWithLightmap(albedo, gLeftWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gRightWallInstanceID)
    {
        litColor = shadeWithLightmap(albedo, gRightWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gBackWallInstanceID)
    {
        litColor = shadeWithLightmap(albedo, gBackWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }
    else if (instanceID == gFrontWallInstanceID)
    {
        litColor = shadeWithLightmap(albedo, gFrontWallLightmap.Sample(gLinearSampler, quadUV).rgb);
    }

    // Box targets.
    else if (instanceID == gTallBoxAInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gTallBoxACenterW, gTallBoxAHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gTallBoxALightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gWideBoxBInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gWideBoxBCenterW, gWideBoxBHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gWideBoxBLightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gBlockCInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gBlockCCenterW, gBlockCHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gBlockCLightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gLowBoxDInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gLowBoxDCenterW, gLowBoxDHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gLowBoxDLightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gThinSlabEInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gThinSlabECenterW, gThinSlabEHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gThinSlabELightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gThinSlabFInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gThinSlabFCenterW, gThinSlabFHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gThinSlabFLightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gTallBoxJInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gTallBoxJCenterW, gTallBoxJHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gTallBoxJLightmap.Sample(gLinearSampler, uv).rgb);
    }
    else if (instanceID == gShortBoxLInstanceID)
    {
        float2 uv = projectCubeFaceUV(vsOut.posW, vsOut.normalW, gShortBoxLCenterW, gShortBoxLHalfExtentW);
        uv.y = 1.0f - uv.y;
        litColor = shadeWithLightmap(albedo, gShortBoxLLightmap.Sample(gLinearSampler, uv).rgb);
    }

    litColor += emissiveSample.rgb;
    return float4(litColor, 1.0f);
}