add_plugin(AdaptiveSHDemo)

target_sources(AdaptiveSHDemo PRIVATE
    AdaptiveSHDemo.cpp
    AdaptiveSHDemo.h
	AdaptiveGridShader.slang
	DynamicPass.slang
	DynamicFilter.cs.slang
	CompositeDynamic.slang
	AdaptiveGridData.slang
	AdaptiveProbeVolume.h
	AdaptiveProbeVolume.cpp
    UniformProbeVolume.h
    UniformProbeVolume.cpp
	ProbeVisualizePass.h
	ProbeVisualizePass.cpp
	UniformGridShader.slang
    AdaptiveGridShaderXAtlasTest.slang
)

target_copy_shaders(AdaptiveSHDemo RenderPasses/AdaptiveSHDemo)
target_source_group(AdaptiveSHDemo "RenderPasses")
