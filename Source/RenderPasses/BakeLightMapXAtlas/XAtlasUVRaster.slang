import Scene.Raster;
#include "BakeDataStructures.slang"

// Page-local compact triangle data.
// All three buffers use the same local triangle index.
//
// gInstanceIDs[i] -> original Falcor geometry instance ID
// gTriangleIDs[i] -> original triangle ID inside that instance's mesh
// gTriangleUVs[i] -> xatlas UV triplet for this physical atlas page
StructuredBuffer<uint> gInstanceIDs;
StructuredBuffer<uint> gTriangleIDs;
StructuredBuffer<TriangleLightmapUV> gTriangleUVs;

struct BakeVSOut
{
    float3 posW : POSW;
    float3 normalW : NORMAL;
    float4 tangentW : TANGENT;
    float2 texC : TEXCRD;
    float4 posH : SV_POSITION;
    nointerpolation uint instanceID : INSTANCE_ID;
    nointerpolation uint triangleID : TRIANGLE_ID;
};

struct PSOut
{
    float4 posW : SV_Target0;
    float4 normW : SV_Target1;
    float4 filterGuide : SV_Target2;
    uint triangleID : SV_Target3;
};

BakeVSOut vsMain(uint vertexID : SV_VertexID)
{
    BakeVSOut o;

    const uint localTriangleID = vertexID / 3;
    const uint cornerID = vertexID % 3;

    const uint sourceInstanceID =
        gInstanceIDs[localTriangleID];
    o.instanceID = sourceInstanceID;

    const uint sourceTriangleID =
        gTriangleIDs[localTriangleID];
    o.triangleID = sourceTriangleID;

    GeometryInstanceID instanceID =
        { sourceInstanceID };

    const uint3 vertexIndices =
        gScene.getIndices(
            instanceID,
            sourceTriangleID
        );

    const uint vertexIndex =
        cornerID == 0 ? vertexIndices.x :
        cornerID == 1 ? vertexIndices.y :
                        vertexIndices.z;

    const StaticVertexData vertex =
        gScene.getVertex(vertexIndex);

    const float4x4 worldMat =
        gScene.getWorldMatrix(instanceID);

    const float3x3 normalMat =
        gScene.getInverseTransposeWorldMatrix(instanceID);

    o.posW =
        mul(
            worldMat,
            float4(vertex.position, 1.f)
        ).xyz;

    o.normalW =
        mul(
            normalMat,
            vertex.normal
        );

    // Material UVs and tangent handedness are independent of xatlas UVs.
    // Match Scene.Raster.defaultVS() when transforming the source tangent.
    o.tangentW = float4(mul((float3x3)worldMat, vertex.tangent.xyz), vertex.tangent.w);
    o.texC = vertex.texCrd;

    const TriangleLightmapUV triUV =
        gTriangleUVs[localTriangleID];

    const float2 uv =
        cornerID == 0 ? triUV.uv0 :
        cornerID == 1 ? triUV.uv1 :
                        triUV.uv2;

    o.posH =
        float4(
            uv * 2.f - 1.f,
            0.f,
            1.f
        );

    return o;
}

PSOut psMain(BakeVSOut pin)
{
    PSOut o;

    // Preserve smooth interpolation; use the mesh front side, never the camera,
    // to orient the receiver. Falcor accounts for the mesh winding convention.
    GeometryInstanceID instanceID = { pin.instanceID };
    float3 Ng = gScene.getFaceNormalW(instanceID, pin.triangleID);
    float normalLengthSquared = dot(pin.normalW, pin.normalW);
    float3 N = normalLengthSquared > 1e-20f ? pin.normalW * rsqrt(normalLengthSquared) : Ng;
    if (dot(N, Ng) < 0.f) N = -N;

    // Filter material textures over the atlas texel footprint.
    let lod = ExplicitGradientTextureSampler(ddx(pin.texC), ddy(pin.texC));
    VertexData v = {};
    v.posW = pin.posW;
    v.normalW = N;
    v.faceNormalW = Ng;
    v.tangentW = float4(pin.tangentW.xyz, sign(pin.tangentW.w));
    v.texC = pin.texC;

    // This bake stores front-side diffuse lighting. A fixed virtual view along
    // Ng makes material evaluation independent of the camera. Do not enable
    // the view-dependent AdjustShadingNormal hint for a baked receiver.
    ShadingData sd = gScene.materials.prepareShadingData(v, gScene.getMaterialID(instanceID), Ng);
    let mi = gScene.materials.getMaterialInstance(sd, lod, (uint)MaterialInstanceHints::None);
    let properties = mi.getProperties(sd);
    // For StandardMaterial, guideNormal is its final normal-mapped shading
    // normal. Falcor handles texture encoding, sampling and the tangent frame.
    float3 mappedNormal = properties.guideNormal;
    float mappedLengthSquared = dot(mappedNormal, mappedNormal);
    if (all(isfinite(mappedNormal)) && mappedLengthSquared > 1e-20f)
        N = mappedNormal * rsqrt(mappedLengthSquared);

    o.posW = float4(pin.posW, 1.f);

    o.normW = float4(N, 1.f);
    // A compact ID lets extraction recover Ng once, without a full float4 MRT.
    o.triangleID = pin.triangleID;
    // World-space texel footprint and exact instance identity for the blur.
    // Derivatives are evaluated during rasterization, before atlas gaps exist.
    o.filterGuide = float4(length(ddx(pin.posW)), length(ddy(pin.posW)),
        float(pin.instanceID & 0xffffu), float(pin.instanceID >> 16));

    return o;
}
