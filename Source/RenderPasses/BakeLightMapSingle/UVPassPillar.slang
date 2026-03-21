import Scene.Raster;

cbuffer BakeCB
{
    uint gReceiverInstanceID;
};

struct PillarVSOut
{
    float4 posH            : SV_Position;
    float3 posW            : POS_W;
    float3 normalW         : NORMAL_W;
    float3 posL            : POS_L;
    float3 normalL : NORMAL_L;
    float2 atlasUV : TEXCOORD5;
    GeometryInstanceID instanceID : INSTANCE_ID;
};

struct PSOut
{
    float4 posW  : SV_Target0;
    float4 normW : SV_Target1;
};

float2 projectCubeFaceUV(float3 posL, float3 nL)
{
    float3 an = abs(nL);

    float2 uv;
    uint faceIndex = 0;

    if (an.x >= an.y && an.x >= an.z)
    {
        if (nL.x > 0.0f)
        {
            uv = float2(-posL.z, posL.y);
            faceIndex = 0; // +X
        }
        else
        {
            uv = float2(posL.z, posL.y);
            faceIndex = 1; // -X
        }
    }
    else if (an.z >= an.x && an.z >= an.y)
    {
        if (nL.z > 0.0f)
        {
            uv = float2(posL.x, posL.y);
            faceIndex = 2; // +Z
        }
        else
        {
            uv = float2(-posL.x, posL.y);
            faceIndex = 3; // -Z
        }
    }
    else
    {
        if (nL.y > 0.0f)
        {
            uv = float2(posL.x, -posL.z);
            faceIndex = 4; // +Y
        }
        else
        {
            uv = float2(posL.x, posL.z);
            faceIndex = 5; // -Y
        }
    }

    // createCube() uses [-0.5, 0.5]
    uv = uv + 0.5;

    uint tileX = faceIndex % 3;
    uint tileY = faceIndex / 3;

    float2 tileScale  = float2(1.0 / 3.0, 1.0 / 2.0);
    float2 tileOffset = float2(tileX, tileY) * tileScale;

    return tileOffset + uv * tileScale;
}

PillarVSOut vsMain(VSIn v)
{
    PillarVSOut o;

    VSOut base = defaultVS(v);

    o.posW = base.posW;
    o.normalW = base.normalW;
    o.instanceID = base.instanceID;

    o.posL = v.pos;
    o.normalL = normalize(v.unpack().normal);

    float2 atlasUV = projectCubeFaceUV(o.posL, o.normalL);
    o.atlasUV = atlasUV;
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
   // o.posW  = float4(pin.posW, 1.0f);
    //o.normW = float4(normalize(pin.normalW), 1.0f);
    o.posW = float4(pin.atlasUV, 0.0, 1.0);
    o.normW = float4(pin.atlasUV, 0.0, 1.0);
    return o;
}
