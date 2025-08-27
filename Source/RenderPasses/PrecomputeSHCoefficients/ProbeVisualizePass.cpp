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
    dsDesc.setDepthEnabled(true);
    dsDesc.setDepthWriteMask(false);
    mpState->setDepthStencilState(pDsState);

     // rasterizer state
    RasterizerState::Desc rasterDesc;
    //rasterDesc.setFillMode(RasterizerState::FillMode::Wireframe);
    rasterDesc.setFillMode(RasterizerState::FillMode::Solid);
    rasterDesc.setCullMode(RasterizerState::CullMode::None);
    rasterDesc.setDepthBias(100000, 1.0f);
    rasterDesc.setDepthClamp(0.0f);
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

void ProbeVisualizePass::setGridData(const ProbeGrid& grid, const std::vector<ProbeDirSample>& dirSamples)
{
    mVertices.clear();

    const int3 res = grid.resolution;

    int width = res.x;
    int height = res.y;
    int depth = res.z;

    // Sphere parameters
    float sphereRadius = 0.15f * std::min({grid.spacing.x, grid.spacing.y, grid.spacing.z});
    //float sphereRadius = 0.3f;
    int segmentsU = 64; // longitude
    int segmentsV = 32;  // latitude
    uint32_t probeIndexCount = 0;
    for (int z = 0; z < depth; ++z)
    {
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                float3 probePos = grid.origin + float3(x * grid.spacing.x, y * grid.spacing.y, z * grid.spacing.z);
                //auto tmp = generateProbeCube(probePos, grid.spacing);

                auto tmp = generateProbeSphere(probePos, sphereRadius, segmentsU, segmentsV, dirSamples, probeIndexCount);
                probeIndexCount ++;
                mVertices.insert(mVertices.end(), tmp.begin(), tmp.end());
            }
        }
    }

     const uint32_t vbSize = (uint32_t)(sizeof(ProbeVertex) * mVertices.size());
     pVertexBuffer = mpDevice->createBuffer(vbSize, ResourceBindFlags::Vertex, MemoryType::Upload, mVertices.data());
     pVertexBuffer->breakStrongReferenceToDevice();

     ref<VertexLayout> pLayout = VertexLayout::create();
     ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();
     pBufLayout->addElement("WORLD_POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
     pBufLayout->addElement("DIR_SAMPLE_INDEX", sizeof(float3), ResourceFormat::R32Uint, 1, 0);
     pBufLayout->addElement("PROBE_INDEX", sizeof(float3) + sizeof(uint32_t), ResourceFormat::R32Uint, 1, 0);
     pLayout->addBufferLayout(0, pBufLayout);

     Vao::BufferVec buffers{pVertexBuffer};
     pVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers);
     mpState->setVao(pVao);
     numSamplePerProbe = (uint32_t)dirSamples.size();
}

void ProbeVisualizePass::setCameraData(const float4x4& viewProjMat, const float4x4& viewMat, const float4x4& projMat)
{
    mpVars->getRootVar()["PerFrameBuffer"]["viewProjMat"] = viewProjMat;
    mpVars->getRootVar()["PerFrameBuffer"]["viewMat"] = viewMat;
    mpVars->getRootVar()["PerFrameBuffer"]["projMat"] = projMat;
}

void ProbeVisualizePass::setProbeSamplingData(ref<Buffer> dirSamples, ref<Buffer> samplingBuffer)
{
    mpVars->getRootVar()["gProbeSamplingResults"] = samplingBuffer;
    mpVars->getRootVar()["gProbeDirSamples"] = dirSamples;
    mpVars->getRootVar()["PerFrameBuffer"]["numSamplePerProbe"] = numSamplePerProbe;
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
std::vector<ProbeVisualizePass::ProbeVertex> ProbeVisualizePass::generateProbeCube(const float3& center, const float3& spacing)
{
    std::vector<ProbeVertex> verts;

    //// half-size along each axis
    //float3 h = spacing * 0.5f;

    //// Cube corners
    //ProbeVertex corners[8] = {
    //    {center + float3(-h.x, -h.y, -h.z)}, // 0
    //    {center + float3(h.x, -h.y, -h.z)},  // 1
    //    {center + float3(-h.x, h.y, -h.z)},  // 2
    //    {center + float3(h.x, h.y, -h.z)},   // 3
    //    {center + float3(-h.x, -h.y, h.z)},  // 4
    //    {center + float3(h.x, -h.y, h.z)},   // 5
    //    {center + float3(-h.x, h.y, h.z)},   // 6
    //    {center + float3(h.x, h.y, h.z)}     // 7
    //};

    // // 12 edges: each defined by a pair of corner indices
    //int edgeIdx[12][2] = {
    //    {0, 1},
    //    {1, 3},
    //    {3, 2},
    //    {2, 0}, // bottom
    //    {4, 5},
    //    {5, 7},
    //    {7, 6},
    //    {6, 4}, // top
    //    {0, 4},
    //    {1, 5},
    //    {2, 6},
    //    {3, 7} // verticals
    //};
    //float thickness = 0.0001f;
    //// For each edge, generate a thin quad along the line
    //for (int e = 0; e < 12; ++e)
    //{
    //    float3 p0 = corners[edgeIdx[e][0]].worldPos;
    //    float3 p1 = corners[edgeIdx[e][1]].worldPos;

    //    // Compute a simple perpendicular offset in screen-aligned direction
    //    // Here we just use a small arbitrary vector for thickness; ideally use camera-facing offset
    //    float3 offset = float3(thickness, thickness, thickness);

    //    // Quad vertices (two triangles)
    //    verts.push_back({p0 - offset});
    //    verts.push_back({p0 + offset});
    //    verts.push_back({p1 - offset});

    //    verts.push_back({p1 - offset});
    //    verts.push_back({p0 + offset});
    //    verts.push_back({p1 + offset});
    //}

    return verts;
}

// New function to generate a UV sphere at a given position
std::vector<ProbeVisualizePass::ProbeVertex> ProbeVisualizePass::generateProbeSphere(
    const float3& center,
    float radius,
    int segmentsU,
    int segmentsV,
    const std::vector<ProbeDirSample>& dirSamples,
    uint32_t probeIndex
)
{
    std::vector<ProbeVertex> verts;

    for (int v = 0; v < segmentsV; ++v)
    {
        float phi0 = float(v) / segmentsV * float(M_PI);
        float phi1 = float(v + 1) / segmentsV * float(M_PI);

        for (int u = 0; u < segmentsU; ++u)
        {
            float theta0 = float(u) / segmentsU * 2.0f * float(M_PI);
            float theta1 = float(u + 1) / segmentsU * 2.0f * float(M_PI);

            // Directions from center to vertex (normalized)
            float3 dir00 = normalize(float3(std::sin(phi0) * std::cos(theta0), std::cos(phi0), std::sin(phi0) * std::sin(theta0)));
            float3 dir01 = normalize(float3(std::sin(phi0) * std::cos(theta1), std::cos(phi0), std::sin(phi0) * std::sin(theta1)));
            float3 dir10 = normalize(float3(std::sin(phi1) * std::cos(theta0), std::cos(phi1), std::sin(phi1) * std::sin(theta0)));
            float3 dir11 = normalize(float3(std::sin(phi1) * std::cos(theta1), std::cos(phi1), std::sin(phi1) * std::sin(theta1)));

            // Find closest sample index for each direction
            uint32_t idx00 = findClosestDirIndex(dir00, dirSamples);
            uint32_t idx01 = findClosestDirIndex(dir01, dirSamples);
            uint32_t idx10 = findClosestDirIndex(dir10, dirSamples);
            uint32_t idx11 = findClosestDirIndex(dir11, dirSamples);

            // Four points of the quad
            float3 p00 = center + radius * dir00;
            float3 p01 = center + radius * dir01;
            float3 p10 = center + radius * dir10;
            float3 p11 = center + radius * dir11;

            // Two triangles per quad
            verts.push_back({p00, idx00, probeIndex});
            verts.push_back({p10, idx10, probeIndex});
            verts.push_back({p11, idx11, probeIndex});

            verts.push_back({p00, idx00, probeIndex});
            verts.push_back({p11, idx11, probeIndex});
            verts.push_back({p01, idx01, probeIndex});
        }
    }
    return verts;
}

// Helper: find the closest direction index
uint32_t ProbeVisualizePass::findClosestDirIndex(const float3& dir, const std::vector<ProbeDirSample>& dirSamples)
{
    uint32_t bestIdx = 0;
    float bestDot = -1.0f;
    for (int i = 0; i < dirSamples.size(); ++i)
    {
        float d = dot(dir, normalize(dirSamples[i].dir));
        if (d > bestDot)
        {
            bestDot = d;
            bestIdx = i;
        }
    }
    return bestIdx;
}
