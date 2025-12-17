#pragma once
#include "Falcor.h"
#include "Core/Pass/BaseGraphicsPass.h"
#include "AdaptiveProbeVolume.h" // Include your volume header

using namespace Falcor;

class ProbeVisualizePass : public BaseGraphicsPass
{
public:
    // Simple vertex for debug rendering
    struct ProbeVertex
    {
        float3 worldPos;
        float3 color;
    };

    static ref<ProbeVisualizePass> create(
        const ref<Device>& pDevice,
        const DefineList& defines = DefineList()
    );

    virtual void execute(RenderContext* pRenderContext, const ref<Fbo>& pFbo, bool autoSetVpSc = true) const;

    // Updated: Takes the linear list of probes from your AdaptiveProbeVolume
    void setVolumeData(const std::vector<AdaptiveProbeVolume::Probe>& probes);

    void setCameraData(const float4x4& viewProjMat);

protected:
    ProbeVisualizePass(const ref<Device>& pDevice, const ProgramDesc& progDesc, const DefineList& programDefines);

private:
    // Generates a wireframe box (thick lines) for a probe
    void generateProbeCube(const float3& minPoint, const float3& maxPoint, const float3& color, std::vector<ProbeVertex>& outVerts);

    ref<RasterizerState> mpRasterState;
    ref<Buffer> pVertexBuffer;
    ref<Vao> pVao;

    std::vector<ProbeVertex> mVertices;
};
