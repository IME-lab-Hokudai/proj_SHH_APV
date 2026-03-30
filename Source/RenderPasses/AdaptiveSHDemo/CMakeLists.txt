add_plugin(AdaptiveSHDemo)

target_sources(AdaptiveSHDemo PRIVATE
    AdaptiveSHDemo.cpp
    AdaptiveSHDemo.h
	StaticPass.slang
	DynamicPass.slang
)

target_copy_shaders(AdaptiveSHDemo RenderPasses/AdaptiveSHDemo)
target_source_group(AdaptiveSHDemo "RenderPasses")
