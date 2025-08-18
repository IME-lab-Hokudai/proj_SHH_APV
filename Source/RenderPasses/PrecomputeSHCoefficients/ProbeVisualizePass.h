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
#pragma once
#include "envMap_SH.h"
#include "Falcor.h"
#include "Core/Pass/BaseGraphicsPass.h"
#include "RenderGraph/RenderPass.h"

using namespace Falcor;

class ProbeVisualizePass : public BaseGraphicsPass
{
struct ProbeVertex
{
    float3 worldPos;
    uint32_t probeSampleIndex;
};

public:

    static ref<ProbeVisualizePass> create(
        const ref<Device>& pDevice,
        const DefineList& defines = DefineList()
    );
    virtual void execute(RenderContext* pRenderContext, const ref<Fbo>& pFbo, bool autoSetVpSc = true) const;

    void setGridData(const ProbeGrid& grid, const std::vector<ProbeDirSample>& dirSamples);

    void setCameraData(const float4x4& viewProjMat, const float4x4& viewMat, const float4x4& projMat);

    void setProbeSamplingData(ref<Buffer> dirSamples,  ref<Buffer> samplingBuffer);

protected:
    ProbeVisualizePass(const ref<Device>& pDevice, const ProgramDesc& progDesc, const DefineList& programDefines);

private:
    std::vector<ProbeVertex> generateProbeCube(const float3& center, const float3& spacing);
    std::vector<ProbeVertex> generateProbeSphere(const float3& center, float radius, int segmentsU, int segmentsV, const std::vector<ProbeDirSample>& dirSamples);
    ref<Program> mpProgram;
    ref<RasterizerState> mpRasterState;
    ref<Buffer> pVertexBuffer;
    ref<Vao> pVao;

    std::vector<ProbeVertex> mVertices;

    uint32_t findClosestDirIndex(const float3& dir, const std::vector<ProbeDirSample>& dirSamples);
};
