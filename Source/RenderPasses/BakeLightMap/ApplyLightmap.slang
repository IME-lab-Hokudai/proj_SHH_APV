import Scene.Raster;

// The baked lightmap (EXR or Result Texture)
Texture2D gLightmap;
SamplerState gLinearSampler;

VSOut vsMain(VSIn vIn)
{
    VSOut vsOut = defaultVS(vIn);
    
    return vsOut;
}

float4 psMain(VSOut vsOut, uint triangleIndex: SV_PrimitiveID) : SV_TARGET
{

    float3 viewDir = normalize(gScene.camera.getPosition() - vsOut.posW);
    ShadingData shadingData = prepareShadingData(vsOut, triangleIndex, viewDir);
    BasicMaterialData baseColorMat = gScene.materials.getBasicMaterialData(vsOut.materialID);
    float4 baseColorSample = float4(0, 0, 0, 1);
    float4 emissiveSample = float4(0, 0, 0, 1);

    if (baseColorMat.texBaseColor.getMode() == TextureHandle::Mode::Texture) {
        baseColorSample = gScene.materials.sampleTexture(baseColorMat.texBaseColor, gLinearSampler, vsOut.texC, 0.f);
        emissiveSample = gScene.materials.sampleTexture(baseColorMat.texEmissive, gLinearSampler, vsOut.texC, 0.f);

    }
    else {
        baseColorSample = (float4)baseColorMat.baseColor;
        emissiveSample = float4(baseColorMat.emissive, 1);
    }
    // We get the base color from the material system to multiply with our light
    float3 albedo = baseColorSample.rgb;

    // 3. Sample the Baked Lightmap
    // We use the UVs that were used during the UVPass rasterization
    float3 irradiance = gLightmap.Sample(gLinearSampler, vsOut.texC).rgb;

    // 4. Final Shading
    // Irradiance represents the total light hitting the surface.
    // For a Lambertian (diffuse) surface: Out = Albedo * Irradiance
    float3 litColor = albedo * irradiance;

    litColor += emissiveSample.rgb;

    return float4(litColor, 1.0f);
}