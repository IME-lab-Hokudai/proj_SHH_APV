import Scene.Raster;

cbuffer BakeCB
{
    uint   gReceiverInstanceID;
    float3 gPillarCenterW;
    float3 gPillarHalfExtentW;
    // Object rotation in degrees, same values as SceneForData.pyscene.
    float3 gPillarRotationEulerDeg;
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

float3 rotateX(float3 v, float a)
{
    float s = sin(a);
    float c = cos(a);
    return float3(v.x, c * v.y - s * v.z, s * v.y + c * v.z);
}

float3 rotateY(float3 v, float a)
{
    float s = sin(a);
    float c = cos(a);
    return float3(c * v.x + s * v.z, v.y, -s * v.x + c * v.z);
}

float3 rotateZ(float3 v, float a)
{
    float s = sin(a);
    float c = cos(a);
    return float3(c * v.x - s * v.y, s * v.x + c * v.y, v.z);
}

float3 worldVectorToBoxLocal(float3 vW)
{
    float3 r = radians(gPillarRotationEulerDeg);

    // Inverse rotation. This assumes the scene transform applies Euler rotation
    // in X -> Y -> Z order. If one object still looks wrong, the rotation order
    // is the first thing to check.
    float3 v = vW;
    v = rotateZ(v, -r.z);
    v = rotateY(v, -r.y);
    v = rotateX(v, -r.x);
    return v;
}

float2 projectCubeFaceUVLocal(float3 posW, float3 normalW)
{
    // Transform world-space point and normal into the box's local orientation.
    float3 pLocal = worldVectorToBoxLocal(posW - gPillarCenterW);
    float3 nLocal = normalize(worldVectorToBoxLocal(normalW));

    // Normalize into local box space [-1, 1].
    float3 p = pLocal / max(gPillarHalfExtentW, float3(1e-6f));
    p = clamp(p, -1.0f, 1.0f);

    float3 n = normalize(nLocal);
    float3 an = abs(n);

    float2 uv;
    uint faceIndex = 0;

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

    uv = uv * 0.5f + 0.5f;

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

    float2 atlasUV = projectCubeFaceUVLocal(o.posW, o.normalW);

    // UV -> clip space.
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
