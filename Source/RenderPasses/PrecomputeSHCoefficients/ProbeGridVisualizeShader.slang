struct VSOut
{
    float4 posH : SV_POSITION;
};

cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
    float4x4 viewMat;
    float4x4 projMat;
}

VSOut vsMain(float4 posW: WORLD_POSITION)
{
    VSOut vout;
    posW.w = 1.0f;
    //posW.z -= 0.1f;//small epsilon to avoid flickering
    vout.posH = mul(viewProjMat, posW);
    return vout;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    float3 litColor = float3(0.5, 0.5, 0.5);
    return float4(litColor,1.0f);
}
