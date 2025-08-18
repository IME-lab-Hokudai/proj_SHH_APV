import ProbeSamplingData;

StructuredBuffer<float4> gProbeSamplingResults;
StructuredBuffer<ProbeDirSample> gProbeDirSamples;

struct VSIn {
    float3 posW : WORLD_POSITION;
    uint32_t probeDirSampleIndex : DIR_SAMPLE_INDEX;
    uint32_t probeIdx : PROBE_INDEX;
};

struct VSOut
{
    float4 posH : SV_POSITION;
    uint32_t probeDirSampleIndex : DIR_SAMPLE_INDEX;
    uint32_t probeIdx : PROBE_INDEX;
};

cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
    float4x4 viewMat;
    float4x4 projMat;
    uint32_t sampleCount;
}

VSOut vsMain(VSIn vIn)
{
    VSOut vOut;
    vOut.posH = mul(viewProjMat, float4(vIn.posW, 1.f));
    vOut.probeDirSampleIndex = vIn.probeDirSampleIndex;
    vOut.probeIdx = vIn.probeIdx;
    return vOut;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    uint flatIdx = vsOut.probeIdx * sampleCount + vsOut.probeDirSampleIndex;
    float4 probeColor = gProbeSamplingResults[flatIdx];
    return float4(probeColor.rgb,1.0f);
}
