struct VSIn {
    float3 posW : WORLD_POSITION;
    float3 color : COLOR;
};

struct VSOut {
    float4 posH : SV_POSITION;
    float3 color : COLOR;
};

cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
}

VSOut vsMain(VSIn vIn)
{
    VSOut vOut;
    vOut.posH = mul(viewProjMat, float4(vIn.posW, 1.f));
    vOut.color = vIn.color;
    return vOut;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    return float4(vsOut.color, 1.0);
}
