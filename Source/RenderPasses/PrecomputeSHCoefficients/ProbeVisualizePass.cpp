/***************************************************************************
 # Copyright (c) 2015-23, NVIDIA CORPORATION. All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions
 # are met:
 #  * Redistributions of source code must retain the above copyright
 #    notice, this list of conditions and the following disclaimer.
 #  * Redistributions in binary form must reproduce the above copyright
 #    notice, this list of conditions and the following disclaimer in the
 #    documentation and/or other materials provided with the distribution.
 #  * Neither the name of NVIDIA CORPORATION nor the names of its
 #    contributors may be used to endorse or promote products derived
 #    from this software without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 # PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **************************************************************************/
#include "ProbeVisualizePass.h"

#include "Scene/TriangleMesh.h"

namespace
{
const char kShaderFile[] = "RenderPasses/PrecomputeSHCoefficients/ProbeGridVisualizeShader.slang";
} // namespace

ProbeVisualizePass::ProbeVisualizePass(const ref<Device>& pDevice, const ProgramDesc& progDesc, const DefineList& programDefines)
    : BaseGraphicsPass(pDevice, progDesc, programDefines)
{
    // default depth stencil state
    DepthStencilState::Desc dsDesc;
    ref<DepthStencilState> pDsState = DepthStencilState::create(dsDesc);
   // auto pDsState = DepthStencilState::create(DepthStencilState::Desc().setDepthEnabled(false));
    mpState->setDepthStencilState(pDsState);

     // rasterizer state
    RasterizerState::Desc rasterDesc;
    rasterDesc.setFillMode(RasterizerState::FillMode::Wireframe);
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
    mpState->setFbo(pFbo, autoSetVpSc);
    pRenderContext->draw(mpState.get(), mpVars.get(), mVertices.size(), 0);
}

void ProbeVisualizePass::setGridData(const ProbeGrid& grid)
{
    mVertices.clear();

    const int3 res = grid.resolution;

    int width = res.x;
    int height = res.y;
    int depth = res.z;

    for (int z = 0; z < depth; ++z)
    {
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                float3 probePos = grid.origin + float3(x * grid.spacing.x, y * grid.spacing.y, z * grid.spacing.z);
                auto tmp = generateProbeCube(probePos, grid.spacing);
                mVertices.insert(mVertices.end(), tmp.begin(), tmp.end());
            }
        }
    }

     const uint32_t vbSize = (uint32_t)(sizeof(ProbeVoxelVertex) * mVertices.size());
     pVertexBuffer = mpDevice->createBuffer(vbSize, ResourceBindFlags::Vertex, MemoryType::Upload, mVertices.data());
     pVertexBuffer->breakStrongReferenceToDevice();

     ref<VertexLayout> pLayout = VertexLayout::create();
     ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();
     pBufLayout->addElement("WORLD_POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
     pLayout->addBufferLayout(0, pBufLayout);

     Vao::BufferVec buffers{pVertexBuffer};
     pVao = Vao::create(Vao::Topology::LineList, pLayout, buffers);
     mpState->setVao(pVao);
}

void ProbeVisualizePass::setCameraData(const float4x4& viewProjMat)
{
    mpVars->getRootVar()["PerFrameBuffer"]["viewProjMat"] = viewProjMat;
}

/*
       v6-------v7
      / |      / |
    v4-------v5  |
    |  |     |   |
    |  v2----|--v3
    | /      | /
    v0-------v1
Index order (triangle strip):
  0, 1, 2, 3,   // front face
  7, 4,         // degenerate
  6, 5,         // right face
  1, 0,         // degenerate
  4, 7, 6, 2,  // top face
  3, 1          // degenerate
  */
std::vector<ProbeVisualizePass::ProbeVoxelVertex> ProbeVisualizePass::generateProbeCube(const float3& center, const float3& spacing)
{
    std::vector<ProbeVoxelVertex> verts;

    // half-size along each axis
    float3 h = spacing * 0.5f;

    // Cube corners
    ProbeVoxelVertex v[8] = {
        {center + float3(-h.x, -h.y, -h.z)}, // 0
        {center + float3(h.x, -h.y, -h.z)},  // 1
        {center + float3(-h.x, h.y, -h.z)},  // 2
        {center + float3(h.x, h.y, -h.z)},   // 3
        {center + float3(-h.x, -h.y, h.z)},  // 4
        {center + float3(h.x, -h.y, h.z)},   // 5
        {center + float3(-h.x, h.y, h.z)},   // 6
        {center + float3(h.x, h.y, h.z)}     // 7
    };

    // Line list order: 12 edges × 2 vertices each
    int edgeIdx[24] = {
        0, 1, 1, 3, 3, 2, 2, 0, // bottom face
        4, 5, 5, 7, 7, 6, 6, 4, // top face
        0, 4, 1, 5, 2, 6, 3, 7  // vertical edges
    };

    verts.reserve(24);
    for (int i = 0; i < 24; i++)
        verts.push_back(v[edgeIdx[i]]);

    return verts;
}
