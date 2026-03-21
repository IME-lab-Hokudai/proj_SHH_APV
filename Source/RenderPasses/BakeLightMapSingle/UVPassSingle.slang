// UVPass.slang
import Scene.Raster;

cbuffer BakeCB
{
    uint gReceiverInstanceID;
    float3 gPillarCenterW;
    float3 gPillarHalfExtentW;
};

struct PSOut
{
    float4 posW : SV_Target0;
    float4 normW : SV_Target1;
};

VSOut vsMain(VSIn v)
{
    VSOut o = defaultVS(v);

    // UV-space rasterization
    float2 uv = v.texC;
    o.posH = float4(uv * 2.f - 1.f, 0.f, 1.f);

    return o;
}

PSOut psMain(VSOut pin)
{
    PSOut o;

    if (pin.instanceID.index != gReceiverInstanceID)
    {
        discard;
    }

    o.posW = float4(pin.posW, 1.f);
    o.normW = float4(normalize(pin.normalW), 1.f);
    return o;
}
