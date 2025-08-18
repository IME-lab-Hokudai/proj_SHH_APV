import ProbeSamplingData;

StructuredBuffer<float4> gProbeSamplingResults;
StructuredBuffer<ProbeDirSample> gProbeDirSamples;

struct VSIn {
    float3 posW : WORLD_POSITION;
    uint32_t probeDirSampleIndex : DIR_SAMPLE_INDEX;
};

struct VSOut
{
    float4 posH : SV_POSITION;
    uint32_t probeSampleIndex : DIR_SAMPLE_INDEX;
};


cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
    float4x4 viewMat;
    float4x4 projMat;
}

VSOut vsMain(VSIn vIn)
{
    VSOut vOut;
    vOut.posH = mul(viewProjMat, float4(vIn.posW, 1.f));
    vOut.probeSampleIndex = vIn.probeDirSampleIndex;
    return vOut;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    float3 litColor = float3(0.5, 0.5, 0.5);
    float4 probeColor = gProbeSamplingResults[vsOut.probeSampleIndex];
    return float4(probeColor.rgb,1.0f);
}
