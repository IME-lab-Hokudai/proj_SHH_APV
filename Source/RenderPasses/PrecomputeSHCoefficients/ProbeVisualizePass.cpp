#include "ProbeVisualizePass.h"

namespace
{
    const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/ProbeGridVisualizeShader.slang";
}

ProbeVisualizePass::ProbeVisualizePass(const ref<Device>& pDevice, const ProgramDesc& progDesc, const DefineList& programDefines)
    : BaseGraphicsPass(pDevice, progDesc, programDefines)
{
    DepthStencilState::Desc dsDesc;
    dsDesc.setDepthEnabled(true);
    dsDesc.setDepthWriteMask(true);
    mpState->setDepthStencilState(DepthStencilState::create(dsDesc));

    RasterizerState::Desc rasterDesc;
    rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
    rasterDesc.setCullMode(RasterizerState::CullMode::Back);
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

    auto var = mpVars->getRootVar()["PerFrameBuffer"];
    var["visibleLevelMask"] = mVisibleLevelMask;
    var["drawLeafOnly"] = mDrawLeafOnly ? 1u : 0u; // Send to GPU

    mpState->setFbo(pFbo, autoSetVpSc);
    pRenderContext->draw(mpState.get(), mpVars.get(), (uint32_t)mVertices.size(), 0);
}

void ProbeVisualizePass::setVolumeData(const std::vector<AdaptiveProbeVolume::Probe>& probes)
{
    mVertices.clear();

    static const float3 kLevelColors[] = {
        float3(0.9f, 0.9f, 0.9f), float3(1.0f, 0.0f, 0.0f), float3(1.0f, 0.5f, 0.0f),
        float3(1.0f, 1.0f, 0.0f), float3(0.0f, 1.0f, 0.0f), float3(0.0f, 1.0f, 1.0f),
        float3(0.0f, 0.0f, 1.0f), float3(1.0f, 0.0f, 1.0f)
    };

    for (const auto& probe : probes)
    {
        int lvl = std::min(probe.level, 7);
        // Pass probe.isLeaf to the generator
        generateProbeCube(probe.minPoint, probe.maxPoint, kLevelColors[lvl], probe.level, probe.isLeaf, mVertices);
    }

    if (mVertices.empty()) return;

    const uint32_t vbSize = (uint32_t)(sizeof(ProbeVertex) * mVertices.size());
    pVertexBuffer = mpDevice->createBuffer(vbSize, ResourceBindFlags::Vertex, MemoryType::Upload, mVertices.data());

    ref<VertexLayout> pLayout = VertexLayout::create();
    ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();

    // Layout Offsets (Cummulative):
    // Pos (12) -> Color (12) -> Level (4) -> IsLeaf (4)
    pBufLayout->addElement("WORLD_POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("COLOR", 12, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("LEVEL_ID", 24, ResourceFormat::R32Uint, 1, 0);
    pBufLayout->addElement("IS_LEAF", 28, ResourceFormat::R32Uint, 1, 0); // New Attribute

    pLayout->addBufferLayout(0, pBufLayout);

    Vao::BufferVec buffers{ pVertexBuffer };
    pVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers);
    mpState->setVao(pVao);
}

void ProbeVisualizePass::setCameraData(const float4x4& viewProjMat)
{
    mpVars->getRootVar()["PerFrameBuffer"]["viewProjMat"] = viewProjMat;
}

void ProbeVisualizePass::generateProbeCube(const float3& minP, const float3& maxP, const float3& color, int level, bool isLeaf, std::vector<ProbeVertex>& outVerts)
{
    // ... setup corners and edges ...
    float3 p0(minP.x, minP.y, minP.z); float3 p1(maxP.x, minP.y, minP.z);
    float3 p2(minP.x, maxP.y, minP.z); float3 p3(maxP.x, maxP.y, minP.z);
    float3 p4(minP.x, minP.y, maxP.z); float3 p5(maxP.x, minP.y, maxP.z);
    float3 p6(minP.x, maxP.y, maxP.z); float3 p7(maxP.x, maxP.y, maxP.z);

    float3 corners[8] = { p0, p1, p2, p3, p4, p5, p6, p7 };
    int edges[12][2] = {
        {0,1}, {2,3}, {4,5}, {6,7}, {0,2}, {1,3}, {4,6}, {5,7}, {0,4}, {1,5}, {2,6}, {3,7}
    };

    float thickness = 0.001f;
    uint32_t leafFlag = isLeaf ? 1 : 0; // Convert to uint

    for (int i = 0; i < 12; ++i)
    {
        // ... build box vertices v0-v7 ...
        float3 start = corners[edges[i][0]];
        float3 end = corners[edges[i][1]];
        float3 dir = normalize(end - start);
        float3 up = (abs(dir.y) < 0.99f) ? float3(0, 1, 0) : float3(1, 0, 0);
        float3 right = normalize(cross(dir, up)) * thickness;
        up = normalize(cross(right, dir)) * thickness;

        float3 v0 = start - right - up; float3 v1 = start + right - up;
        float3 v2 = start + right + up; float3 v3 = start - right + up;
        float3 v4 = end - right - up;   float3 v5 = end + right - up;
        float3 v6 = end + right + up;   float3 v7 = end - right + up;

        auto pushQuad = [&](float3 a, float3 b, float3 c, float3 d) {
            // Push vertex with leafFlag
            outVerts.push_back({ a, color, (uint32_t)level, leafFlag });
            outVerts.push_back({ b, color, (uint32_t)level, leafFlag });
            outVerts.push_back({ c, color, (uint32_t)level, leafFlag });
            outVerts.push_back({ a, color, (uint32_t)level, leafFlag });
            outVerts.push_back({ c, color, (uint32_t)level, leafFlag });
            outVerts.push_back({ d, color, (uint32_t)level, leafFlag });
            };

        pushQuad(v0, v1, v5, v4);
        pushQuad(v1, v2, v6, v5);
        pushQuad(v2, v3, v7, v6);
        pushQuad(v3, v0, v4, v7);
    }
}

// Add this function to the end of the file or alongside setVolumeData
void ProbeVisualizePass::setUniformVolumeData(const float3& minPoint, const float3& cellSize, const uint3& cellResolution)
{
    mVertices.clear();

    // Use a distinct color for Uniform Grid (e.g., Cyan) to distinguish it from Adaptive levels
    float3 color = float3(0.0f, 1.0f, 1.0f);
    uint32_t level = 0;
    bool isLeaf = true; // Uniform cells are always leaves

    // Iterate over all cells in the grid
    for (uint32_t z = 0; z < cellResolution.z; ++z)
    {
        for (uint32_t y = 0; y < cellResolution.y; ++y)
        {
            for (uint32_t x = 0; x < cellResolution.x; ++x)
            {
                // Compute bounds for this specific cell
                float3 currentMin = minPoint + float3(x, y, z) * cellSize;
                float3 currentMax = currentMin + cellSize;

                // Generate the wireframe cube
                generateProbeCube(currentMin, currentMax, color, level, isLeaf, mVertices);
            }
        }
    }

    if (mVertices.empty()) return;

    // -----------------------------------------------------------------------
    // Rebuild GPU Resources (Same logic as setVolumeData)
    // -----------------------------------------------------------------------
    const uint32_t vbSize = (uint32_t)(sizeof(ProbeVertex) * mVertices.size());
    pVertexBuffer = mpDevice->createBuffer(vbSize, ResourceBindFlags::Vertex, MemoryType::Upload, mVertices.data());

    ref<VertexLayout> pLayout = VertexLayout::create();
    ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();

    // Layout Offsets (Cummulative):
    // Pos (12) -> Color (12) -> Level (4) -> IsLeaf (4)
    pBufLayout->addElement("WORLD_POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("COLOR", 12, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("LEVEL_ID", 24, ResourceFormat::R32Uint, 1, 0);
    pBufLayout->addElement("IS_LEAF", 28, ResourceFormat::R32Uint, 1, 0);

    pLayout->addBufferLayout(0, pBufLayout);

    Vao::BufferVec buffers{ pVertexBuffer };
    pVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers);
    mpState->setVao(pVao);
}

void ProbeVisualizePass::setDecisionDebugData(
    const std::vector<AdaptiveProbeVolume::DecisionDebugVoxel>& debugVoxels,
    bool showPruned,
    bool showAdded
)
{
    mVertices.clear();
    static const float3 kLevelColors[] =
    {
        float3(0.9f, 0.9f, 0.9f),
        float3(1.0f, 0.0f, 0.0f),
        float3(1.0f, 0.5f, 0.0f),
        float3(1.0f, 1.0f, 0.0f),
        float3(0.0f, 1.0f, 0.0f),
        float3(0.0f, 1.0f, 1.0f),
        float3(0.0f, 0.0f, 1.0f),
        float3(1.0f, 0.0f, 1.0f)
    };

    for (const auto& v : debugVoxels)
    {
        if (v.kind == AdaptiveProbeVolume::DecisionDebugKind::Pruned)
        {
            if (!showPruned) continue;
        }
        else if (v.kind == AdaptiveProbeVolume::DecisionDebugKind::Added)
        {
            if (!showAdded) continue;
        }
        else
        {
            continue;
        }

        uint32_t lvl = std::min<uint32_t>(v.level, 7);

        float3 color = kLevelColors[lvl];

        // Keep pruned voxels reddish if you want
        if (v.kind == AdaptiveProbeVolume::DecisionDebugKind::Pruned)
            color = lerp(color, float3(1.0f, 0.1f, 0.0f), 0.6f);

        float3 size = v.maxPoint - v.minPoint;
        float3 inflate = 0.01f * size;

        generateProbeCube(
            v.minPoint - inflate,
            v.maxPoint + inflate,
            color,
            int(v.level),
            true,
            mVertices
        );
    }

    if (mVertices.empty()) return;

    const uint32_t vbSize =
        uint32_t(sizeof(ProbeVertex) * mVertices.size());

    pVertexBuffer = mpDevice->createBuffer(
        vbSize,
        ResourceBindFlags::Vertex,
        MemoryType::Upload,
        mVertices.data()
    );

    ref<VertexLayout> pLayout = VertexLayout::create();
    ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();

    pBufLayout->addElement("WORLD_POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("COLOR", 12, ResourceFormat::RGB32Float, 1, 0);
    pBufLayout->addElement("LEVEL_ID", 24, ResourceFormat::R32Uint, 1, 0);
    pBufLayout->addElement("IS_LEAF", 28, ResourceFormat::R32Uint, 1, 0);

    pLayout->addBufferLayout(0, pBufLayout);

    Vao::BufferVec buffers{ pVertexBuffer };
    pVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers);
    mpState->setVao(pVao);
}
