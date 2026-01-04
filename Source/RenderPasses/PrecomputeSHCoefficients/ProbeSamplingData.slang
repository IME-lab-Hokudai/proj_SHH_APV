/***************************************************************************
 # Copyright (c) 2015-21, NVIDIA CORPORATION. All rights reserved.
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
#include "Utils/HostDeviceShared.slangh"

BEGIN_NAMESPACE_FALCOR

/** This is a host/device structure that describes a probe sampling dir.
 */
struct ProbeDirSample {
    float3 dir;
    float dOmega;
};

// -----------------------------------------------
// Probe sampling data shared between CPU and GPU
// -----------------------------------------------
struct ProbeSampleData
{
    float4 Li; // incoming radiance L_i
    float4 s;  // Sample position s_i  Polar coord
    float4 n;  // Sample normal n_i Polar coord
    float hitT; // Ray hit distance t_i to check ray miss/hit
};

struct GPUCorner
{
    // We pack 9 bands. For each band, we have Value + 3 Gradients.
    // float4 allows efficient alignment.
    // .w components can be padding or used for other data.

    float4 coeffs[9]; // .rgb = SH Value
    float4 gradR[9];  // .xyz = Gradient of Red Channel
    float4 gradG[9];  // .xyz = Gradient of Green Channel
    float4 gradB[9];  // .xyz = Gradient of Blue Channel
};

struct GPUProbe
{
    float3 minPoint;
    int pad1;
    float3 maxPoint;
    int pad2;
    int children[8]; // Indices into GPUProbe buffer (-1 if leaf)
    int corners[8];  // Indices into GPUCorner buffer
    // Order: -X, +X, -Y, +Y, -Z, +Z
    int coarseNeighbors[6];
    int pad3[2]; // Keep 16-byte alignment if needed
};

END_NAMESPACE_FALCOR
