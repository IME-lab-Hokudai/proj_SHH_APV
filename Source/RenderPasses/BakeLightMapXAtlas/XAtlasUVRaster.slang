import Scene.Raster;
#include "BakeDataStructures.slang"

// Page-local compact triangle data.
// gTriangleIDs maps a page-local triangle to the original triangle ID in
// the Falcor geometry instance. gTriangleUVs uses the same local index.
StructuredBuffer<uint> gTriangleIDs;
StructuredBuffer<TriangleLightmapUV> gTriangleUVs;

cbuffer BakeCB
{
    uint gReceiverInstanceID;
    uint gTriangleCount;
};

struct BakeVSOut
{
    float3 posW : POSW;
    float3 normalW : NORMAL;
    float4 posH : SV_POSITION;
};

struct PSOut
{
    float4 posW : SV_Target0;
    float4 normW : SV_Target1;
};

BakeVSOut vsMain(uint vertexID : SV_VertexID)
{
    BakeVSOut o;

    const uint localTriangleID = vertexID / 3;
    const uint cornerID = vertexID % 3;

    const uint sourceTriangleID =
        gTriangleIDs[localTriangleID];

    GeometryInstanceID instanceID =
        { gReceiverInstanceID };

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
    o.posW = float4(pin.posW, 1.f);
    o.normW = float4(normalize(pin.normalW), 1.f);
    return o;
}
