#pragma once
#include "Falcor.h"
#include "Core/Pass/BaseGraphicsPass.h"
#include "AdaptiveProbeVolume.h"

using namespace Falcor;

class ProbeVisualizePass : public BaseGraphicsPass
{
public:
    struct ProbeVertex
    {
        float3 worldPos;
        float3 color;
        uint32_t level;
        uint32_t isLeaf; // New member
        uint32_t addedByEdgeAware;
    };

    static ref<ProbeVisualizePass> create(const ref<Device>& pDevice, const DefineList& defines = DefineList());
    virtual void execute(RenderContext* pRenderContext, const ref<Fbo>& pFbo, bool autoSetVpSc = true) const;

    void setVolumeData(const std::vector<AdaptiveProbeVolume::Probe>& probes);
    // New function to visualize uniform grid
    void setUniformVolumeData(const float3& minPoint, const float3& cellSize, const uint3& cellResolution);
    void setCameraData(const float4x4& viewProjMat);

    void toggleLevel(int level, bool visible)
    {
        if (visible) mVisibleLevelMask |= (1u << level);
        else         mVisibleLevelMask &= ~(1u << level);
    }

    // New Setter
    void setDrawLeafOnly(bool enable) { mDrawLeafOnly = enable; }

    void setShowEdgeAdded(bool show)
    {
        mShowEdgeAdded = show;
    }
protected:
    ProbeVisualizePass(const ref<Device>& pDevice, const ProgramDesc& progDesc, const DefineList& programDefines);

private:
    void generateProbeCube(
        const float3& minP,
        const float3& maxP,
        const float3& color,
        int level,
        bool isLeaf,
        bool addedByEdgeAware,
        std::vector<ProbeVertex>& outVerts
    );

    ref<RasterizerState> mpRasterState;
    ref<Buffer> pVertexBuffer;
    ref<Vao> pVao;

    std::vector<ProbeVertex> mVertices;
    uint32_t mVisibleLevelMask = 0xFFFFFFFF;
    bool mDrawLeafOnly = false; // New State
    bool mShowEdgeAdded = true;
};
