struct VSIn
{
    float3 posW : WORLD_POSITION;
    float3 color : COLOR;

    uint level : LEVEL_ID;
    uint isLeaf : IS_LEAF;
    uint edgeAdded : EDGE_ADDED;
};

struct VSOut
{
    float4 posH : SV_POSITION;
    float3 color : COLOR;

    nointerpolation uint level : LEVEL_ID;
    nointerpolation uint isLeaf : IS_LEAF;
    nointerpolation uint edgeAdded : EDGE_ADDED;
};

cbuffer PerFrameBuffer
{
    float4x4 viewProjMat;

    uint32_t numSamplePerProbe;
    uint visibleLevelMask;
    uint drawLeafOnly;
    uint showNormal;
    uint showEdgeAdded;
};

VSOut vsMain(VSIn vIn)
{
    VSOut vOut;

    // Always initialize all output fields.
    vOut.posH = mul(viewProjMat, float4(vIn.posW, 1.0f));
    vOut.color = vIn.color;

    vOut.level = vIn.level;
    vOut.isLeaf = vIn.isLeaf;
    vOut.edgeAdded = vIn.edgeAdded;

    return vOut;
}

float4 psMain(VSOut vIn) : SV_TARGET
{
    // Hide levels disabled by the UI.
    uint levelBit = 1u << min(vIn.level, 31u);

    if ((visibleLevelMask & levelBit) == 0u)
    {
        discard;
    }

    // Hide non-leaf probes when "Draw Leaf Only" is enabled.
    if (drawLeafOnly != 0u && vIn.isLeaf == 0u)
    {
        discard;
    }

    if (vIn.edgeAdded != 0)
    {
        if (showEdgeAdded == 0)
            discard;
    }
    else
    {
        if (showNormal == 0)
            discard;
    }
    return float4(vIn.color, 1.0f);
}
