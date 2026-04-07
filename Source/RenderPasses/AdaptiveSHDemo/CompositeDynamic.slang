SamplerState gLinearSampler;

Texture2D gDynamicFiltered;
Texture2D gDynamicMask;

struct VSOut
{
    float4 posH : SV_Position;
    float2 uv   : TEXCOORD0;
};

VSOut vsMain(uint vertexID : SV_VertexID)
{
    VSOut o;

    float2 pos;
    pos.x = (vertexID == 1 || vertexID == 3) ? 1.0f : -1.0f;
    pos.y = (vertexID >= 2) ? -1.0f : 1.0f;

    o.posH = float4(pos, 0.0f, 1.0f);
    o.uv = pos * float2(0.5f, -0.5f) + 0.5f;

    return o;
}

float4 psMain(VSOut i) : SV_Target
{
    float3 dyn  = gDynamicFiltered.Sample(gLinearSampler, i.uv).rgb;
    float  mask = gDynamicMask.Sample(gLinearSampler, i.uv).r;

    return float4(dyn, mask);
}