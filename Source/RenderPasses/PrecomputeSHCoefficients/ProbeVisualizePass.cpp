#include "ProbeVisualizePass.h"

namespace
{
    const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/ProbeGridVisualizeShader.slang";
}

ProbeVisualizePass::ProbeVisualizePass(const ref<Device>& pDevice, const ProgramDesc& progDesc, const DefineList& programDefines)
    : BaseGraphicsPass(pDevice, progDesc, programDefines)
{
    // Enable Depth Test
    DepthStencilState::Desc dsDesc;
    dsDesc.setDepthEnabled(true);
    dsDesc.setDepthWriteMask(true);
    ref<DepthStencilState> pDsState = DepthStencilState::create(dsDesc);
    mpState->setDepthStencilState(pDsState);

    // Standard Solid Render (we build the "wireframe" as actual geometry)
    RasterizerState::Desc rasterDesc;
    rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
    rasterDesc.setCullMode(RasterizerState::CullMode::None);
    mpRasterState = RasterizerState::create(rasterDesc);
    mpState->setRasterizerState(mpRasterState);
}

ref<ProbeVisualizePass> ProbeVisualizePass::create(const ref<Device>& pDevice, const DefineList& defines)
{
    ProgramDesc desc;
    desc.addShaderLibrary(kShaderFile);
    desc.vsEntry("vsMain");
    desc.psEntry("psMain");
    return ref<ProbeVisualizePass>(new ProbeVisualizePass(pDevice, desc, defines));
}

void ProbeVisualizePass::execute(RenderContext* pRenderContext, const ref<Fbo>& pFbo, bool autoSetVpSc) const
{
    if (mVertices.empty()) return;

    mpState->setFbo(pFbo, autoSetVpSc);
    pRenderContext->draw(mpState.get(), mpVars.get(), (uint32_t)mVertices.size(), 0);
}

void ProbeVisualizePass::setVolumeData(const std::vector<AdaptiveProbeVolume::Probe>& probes)
{
    mVertices.clear();

    // -------------------------------------------------------------
        // VISUALIZATION COLOR PALETTE
        // -------------------------------------------------------------
    static const float3 kLevelColors[] = {
        float3(0.9f, 0.9f, 0.9f), // Level 0: ⬜ White
        float3(1.0f, 0.0f, 0.0f), // Level 1: 🟥 Red
        float3(1.0f, 0.5f, 0.0f), // Level 2: 🟧 Orange
        float3(1.0f, 1.0f, 0.0f), // Level 3: 🟨 Yellow
        float3(0.0f, 1.0f, 0.0f), // Level 4: 🟩 Green
        float3(0.0f, 1.0f, 1.0f), // Level 5: 🟦 Cyan
        float3(0.0f, 0.0f, 1.0f), // Level 6: 🔷 Blue
        float3(1.0f, 0.0f, 1.0f)  // Level 7: 🟪 Magenta
    };

    for (const auto& probe : probes)
    {
        // Only visualize Leaf nodes to see the final tessellation
        if (probe.isLeaf)
        {
            int lvl = std::min(probe.level, 5);
            generateProbeCube(probe.minPoint, probe.maxPoint, kLevelColors[lvl], mVertices);
        }
    }

    if (mVertices.empty()) return;

    // Create Buffers
    const uint32_t vbSize = (uint32_t)(sizeof(ProbeVertex) * mVertices.size());
    pVertexBuffer = mpDevice->createBuffer(vbSize, ResourceBindFlags::Vertex, MemoryType::Upload, mVertices.data());

    ref<VertexLayout> pLayout = VertexLayout::create();
    ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();
    pBufLayout->addElement("WORLD_POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("COLOR", sizeof(float3), ResourceFormat::RGB32Float, 1, 0);
    pLayout->addBufferLayout(0, pBufLayout);

    Vao::BufferVec buffers{ pVertexBuffer };
    pVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers);
    mpState->setVao(pVao);
}

void ProbeVisualizePass::setCameraData(const float4x4& viewProjMat)
{
    mpVars->getRootVar()["PerFrameBuffer"]["viewProjMat"] = viewProjMat;
}

void ProbeVisualizePass::generateProbeCube(const float3& minP, const float3& maxP, const float3& color, std::vector<ProbeVertex>& outVerts)
{
    // Define the 8 corners
    float3 p0 = float3(minP.x, minP.y, minP.z);
    float3 p1 = float3(maxP.x, minP.y, minP.z);
    float3 p2 = float3(minP.x, maxP.y, minP.z);
    float3 p3 = float3(maxP.x, maxP.y, minP.z);
    float3 p4 = float3(minP.x, minP.y, maxP.z);
    float3 p5 = float3(maxP.x, minP.y, maxP.z);
    float3 p6 = float3(minP.x, maxP.y, maxP.z);
    float3 p7 = float3(maxP.x, maxP.y, maxP.z);

    // 12 Edges (Pairs of indices)
    // Bottom Ring: 0-1, 1-3, 3-2, 2-0? No, standard indexing:
    // 0:min, 7:max.
    // X-axis lines, Y-axis lines, Z-axis lines
    float3 corners[8] = { p0, p1, p2, p3, p4, p5, p6, p7 };

    // Indices map based on standard binary order (x, y, z)
    // 0=000, 1=100, 2=010, 3=110, 4=001, 5=101, 6=011, 7=111
    int edges[12][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7}, // X-axis edges
        {0, 2}, {1, 3}, {4, 6}, {5, 7}, // Y-axis edges
        {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Z-axis edges
    };

    // Thickness of the wireframe lines
    float thickness = 0.001f;

    for (int i = 0; i < 12; ++i)
    {
        float3 start = corners[edges[i][0]];
        float3 end = corners[edges[i][1]];

        // Create a thin box (quads) around the segment
        // Simple approach: Add slight offsets to create volume
        float3 dir = normalize(end - start);
        float3 up = (abs(dir.y) < 0.99f) ? float3(0, 1, 0) : float3(1, 0, 0);
        float3 right = normalize(cross(dir, up)) * thickness;
        up = normalize(cross(right, dir)) * thickness;

        // 4 vertices around start, 4 around end
        float3 v0 = start - right - up;
        float3 v1 = start + right - up;
        float3 v2 = start + right + up;
        float3 v3 = start - right + up;

        float3 v4 = end - right - up;
        float3 v5 = end + right - up;
        float3 v6 = end + right + up;
        float3 v7 = end - right + up;

        // Push Triangles (Front, Back, Top, Bottom, Left, Right)
        // Simplified: Just push 2 triangles for the main face facing camera? 
        // No, let's just push a few random faces to make sure it's visible from any angle.
        // Or simpler: Just 2 triangles to form a "billboard" if we were using geometry shader, 
        // but here we are static. Let's do a full box for the line segment.

        auto pushQuad = [&](float3 a, float3 b, float3 c, float3 d) {
            outVerts.push_back({ a, color }); outVerts.push_back({ b, color }); outVerts.push_back({ c, color });
            outVerts.push_back({ a, color }); outVerts.push_back({ c, color }); outVerts.push_back({ d, color });
            };

        // Top/Bottom/Sides
        pushQuad(v0, v1, v5, v4);
        pushQuad(v1, v2, v6, v5);
        pushQuad(v2, v3, v7, v6);
        pushQuad(v3, v0, v4, v7);
    }
}
