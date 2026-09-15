#include "BakeLightMapXAtlas.h"
#include "Utils/Logger.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace
{
const char kSeamShader[] = "RenderPasses/BakeLightMapXAtlas/AtlasSeams.cs.slang";
constexpr uint32_t kInvalid = std::numeric_limits<uint32_t>::max();
constexpr uint32_t kStitchIterations = 32;
constexpr uint32_t kTransferBatchSize = 1u << 20;
constexpr float kEndpointNormalCosine = 0.9999f;
constexpr uint32_t kCorrectionRadius = 6;

struct Edge
{
    uint32_t triangle, c0, c1;
    bool forward;
};
struct EdgeBucket
{
    Edge a{}, b{};
    uint32_t count = 0;
};
struct EdgePair { Edge a, b; bool smooth; };

float2 getUV(const TriangleLightmapUV& uv, uint32_t corner)
{
    return corner == 0 ? uv.uv0 : corner == 1 ? uv.uv1 : uv.uv2;
}

bool sameUV(float2 a, float2 b)
{
    // Connected xatlas corners come from the same output vertex; no tolerance
    // that might accidentally join close but separate charts is needed.
    return a.x == b.x && a.y == b.y;
}

bool aligned(float3 a, float3 b, float threshold)
{
    const float aa = dot(a, a), bb = dot(b, b);
    return std::isfinite(aa) && std::isfinite(bb) && aa > 1e-20f && bb > 1e-20f &&
        dot(a, b) >= threshold * std::sqrt(aa * bb);
}
}

void BakeLightMapXAtlas::prepareAtlasSeams(const std::vector<AtlasPageData>& pages)
{
    mSeamConstraints.clear();
    mSeamTexels.clear();
    mSeamPages.clear();
    mSeamPages.resize(pages.size());
    mSeamTriangleCharts.clear();

    struct TriangleRef { uint32_t page = kInvalid; const AtlasPageTriangleData* triangle = nullptr; };
    std::map<uint32_t, std::vector<TriangleRef>> instances;
    std::vector<uint32_t> instanceIDs;
    for (const auto& page : pages)
        for (const auto& tri : page.triangles)
        {
            auto& refs = instances[tri.instanceID];
            if (refs.size() <= tri.triangleID) refs.resize(size_t(tri.triangleID) + 1);
            refs[tri.triangleID] = {page.pageIndex, &tri};
        }
    if (instances.empty()) return;
    logInfo("Preparing smooth atlas seam connectivity for {} instances (no charting or packing).", instances.size());
    for (const auto& entry : instances) instanceIDs.push_back(entry.first);
    // Also needed when filtering an existing bake without running xatlas.
    bool missingGeometry = false;
    for (uint32_t instanceID : instanceIDs)
        missingGeometry |= mMeshGeometryCache.find(mpScene->getGeometryInstance(instanceID).geometryID) == mMeshGeometryCache.end();
    if (missingGeometry) buildMeshGeometryCache(instanceIDs);

    auto normalPass = ComputePass::create(mpDevice, kSeamShader, "readNormals", mpScene->getSceneDefines());
    std::unordered_map<uint32_t, std::vector<EdgePair>> meshEdges;
    for (uint32_t instanceID : instanceIDs)
    {
        const uint32_t meshID = mpScene->getGeometryInstance(instanceID).geometryID;
        if (meshEdges.find(meshID) != meshEdges.end()) continue;
        auto& geometry = mMeshGeometryCache.at(meshID);
        auto& pairs = meshEdges[meshID];
        if (geometry.positions.empty()) continue;
        if (geometry.normals.empty())
        {
            const uint32_t total = uint32_t(geometry.positions.size());
            auto buffer = mpDevice->createStructuredBuffer(sizeof(float3), std::min(total, kTransferBatchSize));
            auto var = normalPass->getRootVar();
            mpScene->bindShaderData(var["gScene"]);
            var["gNormals"] = buffer;
            geometry.normals.resize(geometry.positions.size());
            for (uint32_t offset = 0; offset < total; offset += kTransferBatchSize)
            {
                const uint32_t count = std::min(kTransferBatchSize, total - offset);
                var["CB"]["gCount"] = count;
                var["CB"]["gVertexOffset"] = mpScene->getMesh(MeshID(meshID)).vbOffset + offset;
                normalPass->execute(mpDevice->getRenderContext(), count, 1);
                buffer->getBlob(geometry.normals.data() + offset, 0, size_t(count) * sizeof(float3));
            }
            var["gNormals"] = ref<Buffer>();
        }

        // Weld exact coincident positions only, within this mesh. UV/tangent
        // splits duplicate vertices; nearby disconnected surfaces are not welded.
        std::map<std::array<float, 3>, uint32_t> welded;
        std::vector<uint32_t> vertexIDs(geometry.positions.size());
        for (uint32_t i = 0; i < geometry.positions.size(); ++i)
        {
            const float3 p = geometry.positions[i];
            FALCOR_CHECK(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z), "Non-finite seam geometry.");
            auto result = welded.emplace(std::array<float, 3>{p.x, p.y, p.z}, uint32_t(welded.size()));
            vertexIDs[i] = result.first->second;
        }
        std::unordered_map<uint64_t, EdgeBucket> edges;
        edges.reserve(geometry.triangles.size() * 2);
        for (uint32_t t = 0; t < geometry.triangles.size(); ++t)
        {
            const uint3 tri = geometry.triangles[t];
            for (uint32_t c = 0; c < 3; ++c)
            {
                uint32_t c0 = c, c1 = (c + 1) % 3;
                uint32_t a = vertexIDs[tri[c0]], b = vertexIDs[tri[c1]];
                if (a == b) continue;
                const bool forward = a < b;
                if (!forward) { std::swap(a, b); std::swap(c0, c1); }
                auto& bucket = edges[(uint64_t(a) << 32) | b];
                Edge edge{t, c0, c1, forward};
                if (bucket.count == 0) bucket.a = edge;
                else if (bucket.count == 1) bucket.b = edge;
                ++bucket.count;
            }
        }
        for (const auto& entry : edges)
        {
            const auto& e = entry.second;
            // Ignore open boundaries, non-manifold edges and coincident faces
            // with incompatible winding rather than guessing their adjacency.
            if (e.count != 2 || e.a.forward == e.b.forward || e.a.triangle == e.b.triangle) continue;
            const uint3 a = geometry.triangles[e.a.triangle], b = geometry.triangles[e.b.triangle];
            const bool smooth = aligned(geometry.normals[a[e.a.c0]], geometry.normals[b[e.b.c0]], kEndpointNormalCosine) &&
                aligned(geometry.normals[a[e.a.c1]], geometry.normals[b[e.b.c1]], kEndpointNormalCosine);
            pairs.push_back({e.a, e.b, smooth});
        }
    }

    std::vector<std::unordered_map<uint32_t, uint32_t>> pixelNodes(pages.size());
    uint32_t nextChart = 0;
    uint64_t seamEdges = 0;
    for (const auto& instance : instances)
    {
        const uint32_t instanceID = instance.first;
        const auto& refs = instance.second;
        const auto& pairs = meshEdges.at(mpScene->getGeometryInstance(instanceID).geometryID);
        std::vector<uint32_t> parent(refs.size());
        std::iota(parent.begin(), parent.end(), 0u);
        auto root = [&](uint32_t i) {
            while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
            return i;
        };
        auto present = [&](const EdgePair& e) {
            return e.a.triangle < refs.size() && e.b.triangle < refs.size() &&
                refs[e.a.triangle].triangle && refs[e.b.triangle].triangle;
        };
        auto continuous = [&](const EdgePair& e) {
            const auto& a = refs[e.a.triangle]; const auto& b = refs[e.b.triangle];
            return a.page == b.page && sameUV(getUV(a.triangle->uv, e.a.c0), getUV(b.triangle->uv, e.b.c0)) &&
                sameUV(getUV(a.triangle->uv, e.a.c1), getUV(b.triangle->uv, e.b.c1));
        };
        // A region is both UV-continuous and smooth. Keeping hard edges out of
        // the union also prevents bilinear seam footprints from modifying their
        // covered texels when an xatlas chart spans a hard shading boundary.
        for (const auto& e : pairs)
            if (e.smooth && present(e) && continuous(e)) parent[root(e.b.triangle)] = root(e.a.triangle);
        auto& charts = mSeamTriangleCharts[instanceID];
        charts.resize(refs.size(), kInvalid);
        std::unordered_map<uint32_t, uint32_t> chartIDs;
        for (uint32_t t = 0; t < refs.size(); ++t)
            if (refs[t].triangle)
            {
                auto entry = chartIDs.emplace(root(t), nextChart);
                if (entry.second) ++nextChart;
                charts[t] = entry.first->second;
            }

        auto pixelPosition = [&](uint32_t page, float2 uv) {
            return float2(uv.x * pages[page].width, (1.f - uv.y) * pages[page].height);
        };
        auto stencil = [&](const TriangleRef& ref, float2 position, SeamStencil& out) {
            const auto& page = pages[ref.page];
            const float2 texel = position - float2(0.5f);
            const int x = int(std::floor(texel.x)), y = int(std::floor(texel.y));
            if (x < 0 || y < 0 || x + 1 >= int(page.width) || y + 1 >= int(page.height)) return false;
            const float fx = texel.x - x, fy = texel.y - y;
            out.weights = {(1.f-fx)*(1.f-fy), fx*(1.f-fy), (1.f-fx)*fy, fx*fy};
            out.chartID = charts[ref.triangle->triangleID];
            for (uint32_t i = 0; i < 4; ++i)
            {
                const uint2 pixel(uint32_t(x) + (i & 1u), uint32_t(y) + (i >> 1));
                const uint32_t index = pixel.y * page.width + pixel.x;
                auto entry = pixelNodes[ref.page].emplace(index, uint32_t(mSeamTexels.size()));
                if (entry.second)
                {
                    SeamTexel empty{};
                    empty.lighting.w = -1.f; // Detect pages that never produced lighting.
                    mSeamTexels.push_back(empty);
                    mSeamPages[ref.page].pixels.push_back(pixel);
                    mSeamPages[ref.page].nodes.push_back(entry.first->second);
                }
                out.nodes[i] = entry.first->second;
            }
            return true;
        };
        for (const auto& e : pairs)
        {
            if (!e.smooth || !present(e) || continuous(e)) continue;
            const auto& a = refs[e.a.triangle]; const auto& b = refs[e.b.triangle];
            const float2 a0 = pixelPosition(a.page, getUV(a.triangle->uv, e.a.c0));
            const float2 a1 = pixelPosition(a.page, getUV(a.triangle->uv, e.a.c1));
            const float2 b0 = pixelPosition(b.page, getUV(b.triangle->uv, e.b.c0));
            const float2 b1 = pixelPosition(b.page, getUV(b.triangle->uv, e.b.c1));
            // Half-texel sampling follows the finer side; endpoints are excluded
            // to avoid forcing a hard corner shared with an unrelated third face.
            const uint32_t count = std::max(1u, uint32_t(std::ceil(2.f * std::max(length(a1-a0), length(b1-b0)))));
            for (uint32_t i = 0; i < count; ++i)
            {
                const float t = (float(i) + 0.5f) / float(count);
                SeamConstraint constraint;
                if (stencil(a, a0 + t * (a1-a0), constraint.a) && stencil(b, b0 + t * (b1-b0), constraint.b))
                    mSeamConstraints.push_back(constraint);
            }
            ++seamEdges;
        }
    }
    if (!mpSeamGatherPass) mpSeamGatherPass = ComputePass::create(mpDevice, kSeamShader, "gather", mpScene->getSceneDefines());
    logInfo("Atlas seams: {} smooth split edges, {} bilinear constraints, {} unique texels across {} pages.",
        seamEdges, mSeamConstraints.size(), mSeamTexels.size(), pages.size());
}

void BakeLightMapXAtlas::gatherAtlasSeams(RenderContext* context, const AtlasPageData& page)
{
    const auto& samples = mSeamPages[page.pageIndex];
    if (samples.pixels.empty()) return;
    const uint32_t total = uint32_t(samples.pixels.size());
    const uint32_t capacity = std::min(total, kTransferBatchSize);
    auto pixels = mpDevice->createStructuredBuffer(sizeof(uint2), capacity, ResourceBindFlags::ShaderResource);
    static_assert(sizeof(SeamTexel) == 40);
    auto output = mpDevice->createStructuredBuffer(sizeof(SeamTexel), capacity);
    auto var = mpSeamGatherPass->getRootVar();
    var["gPixels"] = pixels;
    var["gReadback"] = output;
    var["gInput"] = mpResultTex;
    var["gNormW"] = mpUVFbo->getColorTexture(1);
    var["gFilterGuide"] = mpUVFbo->getColorTexture(2);
    var["gTriangleID"] = mpUVFbo->getColorTexture(3);
    std::vector<SeamTexel> values(capacity);
    for (uint32_t offset = 0; offset < total; offset += kTransferBatchSize)
    {
        const uint32_t count = std::min(kTransferBatchSize, total - offset);
        pixels->setBlob(samples.pixels.data() + offset, 0, size_t(count) * sizeof(uint2));
        var["CB"]["gCount"] = count;
        mpSeamGatherPass->execute(context, count, 1);
        output->getBlob(values.data(), 0, size_t(count) * sizeof(SeamTexel));
        for (uint32_t i = 0; i < count; ++i) mSeamTexels[samples.nodes[offset + i]] = values[i];
    }
    var["gPixels"] = ref<Buffer>(); var["gReadback"] = ref<Buffer>();
    var["gInput"] = ref<Texture>(); var["gNormW"] = ref<Texture>();
    var["gFilterGuide"] = ref<Texture>(); var["gTriangleID"] = ref<Texture>();
}

void BakeLightMapXAtlas::stitchAndSaveAtlasSeams(RenderContext* context, std::vector<AtlasPageData>& pages)
{
    struct ActiveConstraint
    {
        const SeamConstraint* edge;
        std::array<float, 4> moveA, moveB;
    };
    auto movableStencil = [&](const SeamStencil& s, std::array<float, 4>& weights) {
        weights.fill(0.f);
        float movableWeight = 0.f;
        for (uint32_t i = 0; i < 4; ++i)
        {
            if (s.weights[i] <= 1e-6f) continue;
            const auto& sample = mSeamTexels[s.nodes[i]];
            const float3 color = sample.lighting.xyz();
            if (sample.lighting.w < 0.f || !std::isfinite(color.x) || !std::isfinite(color.y) || !std::isfinite(color.z)) return false;
            if (sample.normalCoverage.w > 0.f)
            {
                const auto it = mSeamTriangleCharts.find(sample.owner.x);
                if (it == mSeamTriangleCharts.end() || sample.owner.y >= it->second.size() || it->second[sample.owner.y] != s.chartID)
                    continue; // Keep this tap fixed instead of rejecting the entire seam.
            }
            // Mesh endpoint normals already established smooth connectivity.
            // Comparing normals at different interior texel positions incorrectly
            // rejects curved seams. Gutters may participate even on thin charts.
            weights[i] = s.weights[i];
            movableWeight += weights[i];
        }
        return movableWeight > 1e-4f;
    };
    std::vector<ActiveConstraint> active;
    for (const auto& constraint : mSeamConstraints)
    {
        ActiveConstraint candidate{&constraint, {}, {}};
        if (movableStencil(constraint.a, candidate.moveA) && movableStencil(constraint.b, candidate.moveB))
            active.push_back(candidate);
    }
    if (active.empty())
    {
        logInfo("Atlas seams: no valid smooth seam constraints to apply.");
        mSeamConstraints.clear(); mSeamPages.clear(); mSeamTexels.clear(); mSeamTriangleCharts.clear();
        return;
    }

    // Simultaneous relaxed projections enforce L_A(edge) = L_B(edge) for the
    // actual four-tap bilinear footprints. Both covered and gutter texels are
    // variables. Starting from the filtered bake retains the original detail;
    // only texels participating in seam constraints receive corrections.
    std::vector<float3> values(mSeamTexels.size()), delta(values.size());
    std::vector<uint32_t> uses(values.size(), 0);
    for (size_t i = 0; i < values.size(); ++i)
    {
        const float3 color = mSeamTexels[i].lighting.xyz();
        // Invalid samples are rejected above. Sanitize unused zero-weight taps
        // as well so that 0 * NaN cannot contaminate an otherwise valid stencil.
        values[i] = std::isfinite(color.x) && std::isfinite(color.y) && std::isfinite(color.z) ? color : float3(0.f);
    }
    for (const auto& entry : active)
        for (uint32_t i = 0; i < 4; ++i)
        {
            if (entry.moveA[i] > 0.f) ++uses[entry.edge->a.nodes[i]];
            if (entry.moveB[i] > 0.f) ++uses[entry.edge->b.nodes[i]];
        }
    for (uint32_t iteration = 0; iteration < kStitchIterations; ++iteration)
    {
        std::fill(delta.begin(), delta.end(), float3(0.f));
        for (const auto& entry : active)
        {
            const auto* c = entry.edge;
            float3 error(0.f);
            float denominator = 0.f;
            for (uint32_t i = 0; i < 4; ++i)
            {
                error += c->a.weights[i] * values[c->a.nodes[i]] - c->b.weights[i] * values[c->b.nodes[i]];
                denominator += entry.moveA[i] * entry.moveA[i] + entry.moveB[i] * entry.moveB[i];
            }
            const float3 correction = 0.5f * error / std::max(denominator, 0.125f);
            for (uint32_t i = 0; i < 4; ++i)
            {
                delta[c->a.nodes[i]] -= entry.moveA[i] * correction;
                delta[c->b.nodes[i]] += entry.moveB[i] * correction;
            }
        }
        for (size_t i = 0; i < values.size(); ++i)
            if (uses[i])
            {
                const float3 next = values[i] + delta[i] / float(uses[i]);
                values[i] = float3(std::max(0.f, next.x), std::max(0.f, next.y), std::max(0.f, next.z));
            }
    }

    if (!mpSeamSpreadPass) mpSeamSpreadPass = ComputePass::create(mpDevice, kSeamShader, "spreadCorrection", mpScene->getSceneDefines());
    if (!mpSeamApplyPass) mpSeamApplyPass = ComputePass::create(mpDevice, kSeamShader, "apply", mpScene->getSceneDefines());
    uint32_t maxInstance = 0;
    for (const auto& entry : mSeamTriangleCharts) maxInstance = std::max(maxInstance, entry.first);
    std::vector<uint2> ranges(size_t(maxInstance) + 1, uint2(0));
    std::vector<uint32_t> regions;
    for (const auto& entry : mSeamTriangleCharts)
    {
        ranges[entry.first] = uint2(uint32_t(regions.size()), uint32_t(entry.second.size()));
        regions.insert(regions.end(), entry.second.begin(), entry.second.end());
    }
    auto rangeBuffer = mpDevice->createStructuredBuffer(sizeof(uint2), uint32_t(ranges.size()), ResourceBindFlags::ShaderResource);
    rangeBuffer->setBlob(ranges.data(), 0, ranges.size() * sizeof(uint2));
    auto regionBuffer = mpDevice->createStructuredBuffer(sizeof(uint32_t), uint32_t(regions.size()), ResourceBindFlags::ShaderResource);
    regionBuffer->setBlob(regions.data(), 0, regions.size() * sizeof(uint32_t));
    for (auto& page : pages)
    {
        const auto& samples = mSeamPages[page.pageIndex];
        std::vector<uint2> pixels;
        std::vector<float4> colors;
        for (size_t i = 0; i < samples.nodes.size(); ++i)
            if (uses[samples.nodes[i]])
            {
                pixels.push_back(samples.pixels[i]);
                colors.push_back(float4(values[samples.nodes[i]] - mSeamTexels[samples.nodes[i]].lighting.xyz(), 1.f));
            }
        if (pixels.empty()) continue;
        auto input = Texture::createFromFile(mpDevice, page.outputPath, false, false, ResourceBindFlags::ShaderResource);
        FALCOR_CHECK(input && input->getWidth() == page.width && input->getHeight() == page.height,
            "Cannot load atlas page '{}' for seam stitching.", page.outputPath);
        // Reuse the bake's output allocation; only one disk page is resident.
        if (!mpResultTex || mpResultTex->getWidth() != page.width || mpResultTex->getHeight() != page.height)
            mpResultTex = mpDevice->createTexture2D(page.width, page.height, ResourceFormat::RGBA32Float, 1, 1, nullptr,
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess);
        createAtlasPageGpuBuffers(page);
        rasterizeAtlasPage(context, page);
        // Reuse the old blur texture as sparse correction seeds.
        if (!mpFilteredTex || mpFilteredTex->getWidth() != page.width || mpFilteredTex->getHeight() != page.height)
            mpFilteredTex = mpDevice->createTexture2D(page.width, page.height, ResourceFormat::RGBA32Float, 1, 1, nullptr,
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess);
        context->clearUAV(mpFilteredTex->getUAV().get(), float4(0.f));
        const uint32_t total = uint32_t(pixels.size());
        const uint32_t capacity = std::min(total, kTransferBatchSize);
        auto addresses = mpDevice->createStructuredBuffer(sizeof(uint2), capacity, ResourceBindFlags::ShaderResource);
        auto updates = mpDevice->createStructuredBuffer(sizeof(float4), capacity, ResourceBindFlags::ShaderResource);
        auto apply = mpSeamApplyPass->getRootVar();
        apply["gPixels"] = addresses; apply["gValues"] = updates; apply["gOutput"] = mpFilteredTex;
        for (uint32_t offset = 0; offset < total; offset += kTransferBatchSize)
        {
            const uint32_t count = std::min(kTransferBatchSize, total - offset);
            addresses->setBlob(pixels.data() + offset, 0, size_t(count) * sizeof(uint2));
            updates->setBlob(colors.data() + offset, 0, size_t(count) * sizeof(float4));
            apply["CB"]["gCount"] = count;
            mpSeamApplyPass->execute(context, count, 1);
        }
        auto spread = mpSeamSpreadPass->getRootVar();
        spread["CB"]["gRadius"] = kCorrectionRadius;
        spread["gInput"] = input;
        spread["gSeeds"] = mpFilteredTex;
        spread["gOutput"] = mpResultTex;
        spread["gNormW"] = mpUVFbo->getColorTexture(1);
        spread["gPosW"] = mpUVFbo->getColorTexture(0);
        spread["gFilterGuide"] = mpUVFbo->getColorTexture(2);
        spread["gTriangleID"] = mpUVFbo->getColorTexture(3);
        spread["gChartRanges"] = rangeBuffer;
        spread["gChartIDs"] = regionBuffer;
        mpSeamSpreadPass->execute(context, page.width, page.height);
        // Regenerate outer gutters from the corrected surface. Explicit edge
        // gutter seeds retain coverage and therefore keep the seam constraint.
        if (!mpRawTex || mpRawTex->getWidth() != page.width || mpRawTex->getHeight() != page.height)
            mpRawTex = mpDevice->createTexture2D(page.width, page.height, ResourceFormat::RGBA32Float, 1, 1, nullptr,
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess);
        auto dilate = mpDilatePass->getRootVar();
        dilate["gInput"] = mpResultTex; dilate["gOutput"] = mpRawTex;
        // gPadding retains the configured value from filterAndSaveAtlasPage().
        mpDilatePass->execute(context, page.width, page.height);
        mpRawTex->captureToFile(0, 0, page.outputPath, Bitmap::FileFormat::ExrFile, Bitmap::ExportFlags::Uncompressed, false);
        for (const char* name : {"gInput", "gSeeds", "gOutput", "gNormW", "gPosW", "gFilterGuide", "gTriangleID"}) spread[name] = ref<Texture>();
        spread["gChartRanges"] = ref<Buffer>(); spread["gChartIDs"] = ref<Buffer>();
        dilate["gInput"] = ref<Texture>(); dilate["gOutput"] = ref<Texture>();
        apply["gPixels"] = ref<Buffer>(); apply["gValues"] = ref<Buffer>(); apply["gOutput"] = ref<Texture>();
        releaseAtlasPageGpuBuffers(page);
        input = nullptr; addresses = nullptr; updates = nullptr;
        mpDevice->wait();
        logInfo("Seam-stitched atlas page {}: {} boundary/gutter texels updated.", page.pageIndex, pixels.size());
    }
    logInfo("Atlas seam stitching complete: {} / {} constraints accepted, {} iterations.",
        active.size(), mSeamConstraints.size(), kStitchIterations);
    mSeamConstraints.clear(); mSeamPages.clear(); mSeamTexels.clear(); mSeamTriangleCharts.clear();
}
