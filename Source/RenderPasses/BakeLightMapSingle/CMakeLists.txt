add_plugin(BakeLightMapSingle)

target_sources(BakeLightMapSingle PRIVATE
    BakeLightMapSingle.cpp
    BakeLightMapSingle.h
	UVPassSingle.slang
    BakeDataStructures.slang
	ExtractTexelsSingle.cs.slang
	LightmapBakingSingle.rt.slang
	NormalizeLightmapSingle.cs.slang
	ApplyLightmapSingle.slang
)

target_copy_shaders(BakeLightMapSingle RenderPasses/BakeLightMapSingle)
target_source_group(BakeLightMapSingle "RenderPasses")
