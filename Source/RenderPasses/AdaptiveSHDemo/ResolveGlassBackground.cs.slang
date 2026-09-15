// Resolve the existing four-sample scene targets into immutable glass inputs.
Texture2DMS<float4> gColorMS;
Texture2DMS<float> gDepthMS;
RWTexture2D<float4> gColor;
RWTexture2D<float> gDepth;

[numthreads(16, 16, 1)]
void main(uint3 id : SV_DispatchThreadID)
{
    uint width, height, samples;
    gColorMS.GetDimensions(width, height, samples);
    if (id.x >= width || id.y >= height) return;
    float4 color = float4(0.f);
    float depth = 1.f;
    for (uint i = 0; i < samples; ++i)
    {
        color += gColorMS.Load(id.xy, i);
        depth = min(depth, gDepthMS.Load(id.xy, i));
    }
    gColor[id.xy] = color / float(samples);
    gDepth[id.xy] = depth;
}
