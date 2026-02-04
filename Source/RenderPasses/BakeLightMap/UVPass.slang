// UVPass.slang
import Scene.Raster;

// vsMain uses the built-in VSIn and VSOut from Raster.slang
VSOut vsMain(VSIn v)
{
    VSOut vout;
    const GeometryInstanceID instanceID = { v.instanceID };

    // Map UVs [0, 1] to NDC [-1, 1] for lightmap unwrapping
    float2 uv = v.texC;
    vout.posH = float4(uv.x * 2.0f - 1.0f, uv.y * 2.0f - 1.0f, 0.0f, 1.0f);

    // Compute standard world space attributes
    float4x4 worldMat = gScene.getWorldMatrix(instanceID);
    vout.posW = mul(worldMat, float4(v.pos, 1.f)).xyz;
    
    // Unpack and transform normal to world space
    float3 normal = v.unpack().normal;
    vout.normalW = mul(gScene.getInverseTransposeWorldMatrix(instanceID), normal);
    //float3x3 invTrans = (float3x3)gScene.getInverseTransposeWorldMatrix(instanceID);
    //vout.normalW = normalize(mul(invTrans, float3(0, 0, 1))); // test with fixed input
    // Pass through IDs for material/instance lookup if needed
    vout.instanceID = instanceID;
    vout.materialID = gScene.getMaterialID(instanceID);
    vout.texC = uv;

    return vout;
}

struct PSOut
{
    float4 posW : SV_Target0;
    float4 normW : SV_Target1;
};

PSOut psMain(VSOut pin) : SV_Target
{
    PSOut pout;
    // .w = 1.0f marks this texel as valid for the extraction compute shader
    pout.posW = float4(pin.posW, 1.0f);
    //pout.posW = float4(1, 0, 0, 1);
    pout.normW = float4(normalize(pin.normalW), 1.0f);
    
    return pout;
}