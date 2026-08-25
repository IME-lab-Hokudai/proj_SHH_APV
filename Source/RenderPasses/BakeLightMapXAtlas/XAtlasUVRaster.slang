import Scene.Raster;

#include "BakeDataStructures.slang"

StructuredBuffer<TriangleLightmapUV> gTriangleUVs;

cbuffer BakeCB
{
    uint gReceiverInstanceID;
    uint gTriangleCount;
};

struct BakeVSOut
{
    float3 posW    : POSW;
    float3 normalW : NORMAL;

    float4 posH : SV_POSITION;
};

struct PSOut
{
    float4 posW  : SV_Target0;
    float4 normW : SV_Target1;
};


// =============================================================
// Vertex shader
//
// We will issue a NON-INDEXED draw with:
//
//     vertexCount = triangleCount * 3
//
// Therefore:
//
//     vertexID / 3 = triangle ID
//     vertexID % 3 = corner ID
//
// Example:
//
//     vertexID 0 -> triangle 0, corner 0
//     vertexID 1 -> triangle 0, corner 1
//     vertexID 2 -> triangle 0, corner 2
//     vertexID 3 -> triangle 1, corner 0
//     ...
// =============================================================

BakeVSOut vsMain(uint vertexID : SV_VertexID)
{
    BakeVSOut o;

    const uint triangleID =
        vertexID / 3;

    const uint cornerID =
        vertexID % 3;

    const GeometryInstanceID instanceID =
        { gReceiverInstanceID };

    // ---------------------------------------------------------
    // Get the ORIGINAL Falcor vertex indices for this triangle.
    //
    // These are global indices into gScene.vertices.
    // ---------------------------------------------------------

    const uint3 vertexIndices =
        gScene.getIndices(
            instanceID,
            triangleID
        );

    uint vertexIndex;

    if (cornerID == 0)
    {
        vertexIndex = vertexIndices.x;
    }
    else if (cornerID == 1)
    {
        vertexIndex = vertexIndices.y;
    }
    else
    {
        vertexIndex = vertexIndices.z;
    }

    // ---------------------------------------------------------
    // Fetch the original Falcor vertex.
    //
    // This gives us its original:
    //     position
    //     normal
    //     tangent
    //     texcoord
    // etc.
    // ---------------------------------------------------------

    const StaticVertexData vertex =
        gScene.getVertex(vertexIndex);

    // ---------------------------------------------------------
    // Transform original geometry into world space.
    // ---------------------------------------------------------

    const float4x4 worldMat =
        gScene.getWorldMatrix(instanceID);

    const float3x3 normalMat =
        gScene.getInverseTransposeWorldMatrix(
            instanceID
        );

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

    // ---------------------------------------------------------
    // Get xatlas UV for this triangle corner.
    // ---------------------------------------------------------

    const TriangleLightmapUV triUV =
        gTriangleUVs[triangleID];

    float2 uv;

    if (cornerID == 0)
    {
        uv = triUV.uv0;
    }
    else if (cornerID == 1)
    {
        uv = triUV.uv1;
    }
    else
    {
        uv = triUV.uv2;
    }

    // ---------------------------------------------------------
    // Rasterize directly in xatlas UV space.
    // ---------------------------------------------------------

    o.posH =
        float4(
            uv * 2.f - 1.f,
            0.f,
            1.f
        );

    return o;
}


// =============================================================
// Pixel shader
// =============================================================

PSOut psMain(BakeVSOut pin)
{
    PSOut o;

    o.posW =
        float4(pin.posW, 1.f);

    o.normW =
        float4(
            normalize(pin.normalW),
            1.f
        );

    return o;
}
