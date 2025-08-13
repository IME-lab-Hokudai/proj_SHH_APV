struct VSOut
{
    float4 posH : SV_POSITION;
};

cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
}

VSOut vsMain(float4 posW: WORLD_POSITION)
{
    VSOut vOut;
    posW.w = 1.0f;
    vOut.posH = mul(viewProjMat, posW);
    return vOut;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    float3 litColor = float3(0.5, 0.5, 0.5);
    return float4(litColor,1.0f);
}
