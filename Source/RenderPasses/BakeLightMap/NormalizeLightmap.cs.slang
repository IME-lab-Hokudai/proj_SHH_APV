// NormalizeLightmap.cs.slang

// The 1D raw accumulation buffer from the Ray Tracer
StructuredBuffer<float4> gAccumBuffer; 

// The final 2D texture for the RenderGraph output
RWTexture2D<float4> gOutput; 

cbuffer CB {
    uint gTotalSamples;
    uint gWidth;  // We need the width to calculate the 1D index
};

[numthreads(16, 16, 1)]
void main(uint3 id : SV_DispatchThreadID) {
    uint2 res;
    gOutput.GetDimensions(res.x, res.y);
    
    // Bounds check for the 2D texture
    if (any(id.xy >= res)) return;

    // Convert 2D pixel coordinate (x, y) to 1D buffer index
    uint flatIndex = id.y * gWidth + id.x;
    
    // Read raw energy from the 1D buffer
    float4 accumulated = gAccumBuffer[flatIndex];
    
    // 1. Divide by total samples to get the average 
    // 2. Multiply by PI to handle Cosine Hemisphere PDF normalization 
    float3 finalIrradiance = (accumulated.rgb / float(gTotalSamples)) * 3.14159265f;

    // Write the normalized result to the 2D output texture [cite: 114]
    gOutput[id.xy] = float4(finalIrradiance, 1.0f);
}