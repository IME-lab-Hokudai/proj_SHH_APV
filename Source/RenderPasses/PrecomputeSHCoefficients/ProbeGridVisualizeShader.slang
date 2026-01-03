struct VSIn {
    float3 posW : WORLD_POSITION;
    float3 color : COLOR;
    uint level : LEVEL_ID;
    uint isLeaf : IS_LEAF; // New Attribute
};

struct VSOut {
    float4 posH : SV_POSITION;
    float3 color : COLOR;
};

cbuffer PerFrameBuffer {
    float4x4 viewProjMat;
    uint32_t numSamplePerProbe;
    uint visibleLevelMask; 
    uint drawLeafOnly; 
}

VSOut vsMain(VSIn vIn)
{
    VSOut vOut;
    
    // 1. Level Mask Check
    uint mask = 1u << vIn.level;
    if ((visibleLevelMask & mask) == 0)
    {
        vOut.posH = float4(-1, -1, -1, -1); // Discard
        vOut.color = float3(0,0,0);
        return vOut;
    }

    // 2. Leaf Only Check
    // If "Draw Leaf Only" is ON, and this node is NOT a leaf, discard it.
    if (drawLeafOnly != 0 && vIn.isLeaf == 0)
    {
        vOut.posH = float4(-1, -1, -1, -1); // Discard
        vOut.color = float3(0,0,0);
        return vOut;
    }

    vOut.posH = mul(viewProjMat, float4(vIn.posW, 1.f));
    vOut.color = vIn.color;
    return vOut;
}

float4 psMain(VSOut vsOut) : SV_TARGET
{
    return float4(vsOut.color, 1.0);
}
