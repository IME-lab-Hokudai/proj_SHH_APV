// RGB = accumulated lighting; W = number of traced samples.
StructuredBuffer<float4> gAccumBuffer;
RWTexture2D<float4> gOutput;

cbuffer CB
{
    uint gTotalSamples;
    uint gWidth;
};

[numthreads(16, 16, 1)]
void main(uint3 id : SV_DispatchThreadID)
{
    uint2 res;
    gOutput.GetDimensions(res.x, res.y);
    if (any(id.xy >= res)) return;

    uint flatIndex = id.y * gWidth + id.x;
    float4 accumulated = gAccumBuffer[flatIndex];

    // Preserve coverage for the later blur and dilation passes.
    float3 lighting = accumulated.rgb / float(max(gTotalSamples, 1u));
    gOutput[id.xy] = float4(lighting, accumulated.w > 0.f ? 1.f : 0.f);
}
