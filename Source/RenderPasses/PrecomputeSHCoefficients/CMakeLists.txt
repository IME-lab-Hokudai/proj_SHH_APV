add_plugin(PrecomputeSHCoefficients)

# Probe tracing imports Falcor's original PathTracer shader modules.
add_dependencies(PrecomputeSHCoefficients PathTracer)

target_sources(PrecomputeSHCoefficients PRIVATE
    CommonDefines.h
    PrecomputeSHCoefficients.cpp
    PrecomputeSHCoefficients.h
	sphericalHarmonics.h
	sphericalHarmonics.cpp
    envMap_SH.h
	envMap_SH.cpp
	AdaptiveProbeVolume.h
	AdaptiveProbeVolume.cpp
    UniformProbeVolume.h
    UniformProbeVolume.cpp
	ProbeVisualizePass.h
	ProbeVisualizePass.cpp
    SHShader.slang
	EnvMapShader.slang
	SHGridShader.slang
	ProbeGridVisualizeShader.slang
	ProbeSamplingData.slang
	ProbeSampling.rt.slang
    ProbePathTracer.slang
    SHAdaptiveProbeShader.slang
)

target_copy_shaders(PrecomputeSHCoefficients RenderPasses/PrecomputeSHCoefficients)
target_source_group(PrecomputeSHCoefficients "RenderPasses")
