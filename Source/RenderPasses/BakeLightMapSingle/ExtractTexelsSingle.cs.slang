// ExtractTexels.cs.slang
#include "BakeDataStructures.slang"

Texture2D<float4> gPosW;
Texture2D<float4> gNormW;
RWStructuredBuffer<TexelSample> gTexelSamples;
RWStructuredBuffer<uint> gCounter; // Resets to 0 in C++ before execution

cbuffer CB {
    uint2 gResolution;
};

[numthreads(16, 16, 1)]
void main(uint3 id : SV_DispatchThreadID) {
    if (any(id.xy >= gResolution)) return;

    float4 pos = gPosW[id.xy];
    
    // pos.w = 1.0f was set in UVPass.slang for valid geometry
    if (pos.w > 0.0f) {
        uint bufferIdx;
        // Atomically increment the counter to find the next slot in the 1D buffer
        InterlockedAdd(gCounter[0], 1, bufferIdx);
        
        TexelSample ts;
        ts.posW = pos.xyz;
        ts.normalW = normalize(gNormW[id.xy].xyz);
        ts.texelIndex = id.y * gResolution.x + id.x; // Original pixel address
        
        gTexelSamples[bufferIdx] = ts;
    }
}
