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
    float4 irradiance: IRRADIANCE;
};

cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
    float4x4 viewMat;
    float4x4 projMat;
    uint32_t numSamplePerProbe;
}

VSOut vsMain(VSIn vIn)
{
    VSOut vOut;
    vOut.posH = mul(viewProjMat, float4(vIn.posW, 1.f));
    vOut.probeDirSampleIndex = vIn.probeDirSampleIndex;
    vOut.probeIdx = vIn.probeIdx;
    uint flatIdx = vIn.probeIdx * numSamplePerProbe + vIn.probeDirSampleIndex;
    //vOut.irradiance = abs(gProbeSamplingResults[flatIdx]);
    vOut.irradiance = abs(gProbeSamplingResults[flatIdx]);
    return vOut;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    // uint flatIdx = vsOut.probeIdx * sampleCount + vsOut.probeDirSampleIndex;
    // float4 probeColor = gProbeSamplingResults[flatIdx];
    return float4(vsOut.irradiance.rgb,1.0f);
}
