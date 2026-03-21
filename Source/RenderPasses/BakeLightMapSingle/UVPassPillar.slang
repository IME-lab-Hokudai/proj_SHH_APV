import Scene.Raster;

cbuffer BakeCB
{
    uint   gReceiverInstanceID;
    float3 gPillarCenterW;
    float3 gPillarHalfExtentW;
};

struct PillarVSOut
{
    float4 posH : SV_Position;
    float3 posW : TEXCOORD0;
    float3 normalW : TEXCOORD1;
    GeometryInstanceID instanceID : INSTANCE_ID;
};

struct PSOut
{
    float4 posW  : SV_Target0;
    float4 normW : SV_Target1;
};

float2 projectCubeFaceUV(float3 posW, float3 normalW)
{
    // Normalize into pillar-local space [-1,1]
    float3 p = (posW - gPillarCenterW) / gPillarHalfExtentW;

    float3 n = normalize(normalW);
    float3 an = abs(n);

    float2 uv;
    uint faceIndex = 0;

    // Face selection
    if (an.x >= an.y && an.x >= an.z)
    {
        if (n.x > 0.0f)
        {
            uv = float2(-p.z, p.y);  // +X
            faceIndex = 0;
        }
        else
        {
            uv = float2(p.z, p.y);   // -X
            faceIndex = 1;
        }
    }
    else if (an.z >= an.x && an.z >= an.y)
    {
        if (n.z > 0.0f)
        {
            uv = float2(p.x, p.y);   // +Z
            faceIndex = 2;
        }
        else
        {
            uv = float2(-p.x, p.y);  // -Z
            faceIndex = 3;
        }
    }
    else
    {
        if (n.y > 0.0f)
        {
            uv = float2(p.x, -p.z);  // +Y
            faceIndex = 4;
        }
        else
        {
            uv = float2(p.x, p.z);   // -Y
            faceIndex = 5;
        }
    }

    // Map [-1,1] → [0,1]
    uv = uv * 0.5f + 0.5f;

    // Pack into 3x2 atlas
    uint tileX = faceIndex % 3;
    uint tileY = faceIndex / 3;

    float2 tileScale  = float2(1.0f / 3.0f, 1.0f / 2.0f);
    float2 tileOffset = float2(tileX, tileY) * tileScale;

    return tileOffset + uv * tileScale;
}

PillarVSOut vsMain(VSIn v)
{
    PillarVSOut o;

    VSOut base = defaultVS(v);

    o.posW = base.posW;
    o.normalW = normalize(base.normalW);
    o.instanceID = base.instanceID;

    float2 atlasUV = projectCubeFaceUV(o.posW, o.normalW);

    // UV → clip space
    o.posH = float4(atlasUV * 2.0f - 1.0f, 0.0f, 1.0f);

    return o;
}

PSOut psMain(PillarVSOut pin)
{
    if (pin.instanceID.index != gReceiverInstanceID)
    {
        discard;
    }

    PSOut o;
    o.posW  = float4(pin.posW, 1.0f);
    o.normW = float4(normalize(pin.normalW), 1.0f);
    return o;
}
