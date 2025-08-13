add_plugin(PrecomputeSHCoefficients)

target_sources(PrecomputeSHCoefficients PRIVATE
    CommonDefines.h
    PrecomputeSHCoefficients.cpp
    PrecomputeSHCoefficients.h
	sphericalHarmonics.h
	sphericalHarmonics.cpp
    envMap_SH.h
	envMap_SH.cpp
    SHShader.slang
	EnvMapShader.slang
	SHGridShader.slang
	ProbeVisualizePass.h
	ProbeVisualizePass.cpp
	ProbeGridVisualizeShader.slang
)

target_copy_shaders(PrecomputeSHCoefficients RenderPasses/PrecomputeSHCoefficients)
target_source_group(PrecomputeSHCoefficients "RenderPasses")
