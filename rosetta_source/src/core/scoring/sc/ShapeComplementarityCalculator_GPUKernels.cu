// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file     core/scoring/sc/ShapeComplementarityCalculator_CPUKernels.cu
/// @brief    GPU Kernel code. This file needs to be compiled with NVIDIA's nvcc
/// @detailed Lawrence & Coleman shape complementarity calculator (based on CCP4's sc)
/// @author   Luki Goldschmidt <luki@mbi.ucla.edu>

/// This code was ported from the original Fortran code found in CCP4:
/// Sc (Version 2.0): A program for determining Shape Complementarity
/// Copyright Michael Lawrence, Biomolecular Research Institute
/// 343 Royal Parade Parkville Victoria Australia
///
/// This version contains support for GPU-acceleration by CUDA-capable devices,
/// which provides a 10-25x speed up over the CPU-only code using a regular desktop
/// video card with 4 processors (32 cores). Define CUDA_CPU and compile with

#include <cuda_runtime_api.h>

#define MIN(a,b) ((a) < (b) ? (a): (b))

//////////////////////////////////////////////////////////////////////
// Collision checking GPU kernel for TrimPeripheralBand

__global__ void _cuda_TrimPeripheralBand_kernel(
	float3 *dAccDotCoords,
	uint nAcc,
	float3 *dBurDotCoords,
	char *dDotColl,
	float r2)
{
	register int i, j, l;
	register float3 dot1;
	__shared__ char sColl[1024];
	__shared__ float3 sCoords[1024];

	sColl[threadIdx.x] = 0;
 	dot1 = dBurDotCoords[blockIdx.x*blockDim.x + threadIdx.x];

	for(i = 0; i < nAcc; i += blockDim.x) {
		__syncthreads();
		sCoords[threadIdx.x] = dAccDotCoords[i + threadIdx.x];
		__syncthreads();

		l = MIN(nAcc - i, blockDim.x);
		for(j = 0; j < l; j++) {
			register float3 dot2 = sCoords[j];
			dot2.x -= dot1.x;
			dot2.y -= dot1.y;
			dot2.z -= dot1.z;
			sColl[threadIdx.x] |= (dot2.x*dot2.x + dot2.y*dot2.y + dot2.z*dot2.z) <= r2;
		}
	}
	dDotColl[blockIdx.x*blockDim.x + threadIdx.x] = sColl[threadIdx.x];
}

//////////////////////////////////////////////////////////////////////
// Finding closest dot neighbor GPU kernel

__global__ void _cuda_FindClosestNeighbor_kernel(
	float3 *dMyDotCoords,
	float3 *dTheirDotCoords,
	uint nTheirDots,
	uint *dNeighbors)
{
	register int i, j, l;
	register float3 dot1;
	__shared__ uint sNeighbors[512];
	__shared__ float3 sCoords[512];
	float distmin = 99999.0, d2;

 	dot1 = dMyDotCoords[blockIdx.x*blockDim.x + threadIdx.x];

	for(i = 0; i < nTheirDots; i += blockDim.x) {
		__syncthreads();
		sCoords[threadIdx.x] = dTheirDotCoords[i + threadIdx.x];
		__syncthreads();

		l = MIN(nTheirDots - i, blockDim.x);
		for(j = 0; j < l; j++) {
			register float3 dot2 = sCoords[j];
			dot2.x -= dot1.x;
			dot2.y -= dot1.y;
			dot2.z -= dot1.z;
			d2 = dot2.x*dot2.x + dot2.y*dot2.y + dot2.z*dot2.z;
			if(d2 <= distmin) {
				distmin = d2;
				sNeighbors[threadIdx.x] = i+j;
			}
		}
	}
	dNeighbors[blockIdx.x*blockDim.x + threadIdx.x] = sNeighbors[threadIdx.x];
}

//////////////////////////////////////////////////////////////////////
// Stubs called from CPU code

void _cuda_sccalc_TrimPeripheralBand(int x, int y, float3 *dAccDotCoords, uint nAcc, float3 *dBurDotCoords, char *dDotColl, float r2)
{
	_cuda_TrimPeripheralBand_kernel<<<x, y>>>(dAccDotCoords, nAcc, dBurDotCoords, dDotColl, r2);
}

void _cuda_sccalc_FindClosestNeighbor(int x, int y, float3 *dMyDotCoords, float3 *dTheirDotCoords, uint nTheirDotsCoords, uint *dNeighbors)
{
	_cuda_FindClosestNeighbor_kernel<<<x, y>>>(dMyDotCoords, dTheirDotCoords, nTheirDotsCoords, dNeighbors);
}

