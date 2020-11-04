
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <cfloat>
#include <cstdint>

#include "math_constants.h"

// Open Mesh
#include "MeshHeader.h"
// Size Function
#include "../BackgroundGrid/SizeGrid.h"

/* Keep NQ = 8 for two dimensional meshes! This value was chosen because it gives optimal
performance considering a warp-size of 32 because NQ = 8 results in 8 * 8 = 64 nodes
which is double the warp size. Each vertex is computed using one warp where each warp
computes two grid nodes.
Another implementation used 2 warps for one grid but it was slower as syncthreads is
too expensive.
*/
// Size of Quality Mesh
constexpr int NQ = 8;
// number of refinement steps within DMO
constexpr int  DMO_DEPTH = 3;
// double the maximal number of allowed vertices on the one-ring neighborhood
constexpr int  MAX_ONE_RING_SIZE = 64;


// Error output
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		//fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		fprintf(stderr, "GPUassert: %s. Line %d\n", cudaGetErrorString(code), line);
		if (abort) exit(code);
	}
}


typedef union {
	float floats[2];                 // floats[0] = lowest
	int32_t ints[2];                     // ints[1] = lowIdx
	uint64_t ulong;    // for atomic update
} my_atomics;

__device__ uint64_t my_atomicArgMax(uint64_t* address, float val, int32_t idx)
{
	my_atomics loc, newValue;
	loc.floats[0] = val;
	loc.ints[1] = idx;
	newValue.ulong = *address;
	while (newValue.floats[0] < val)
		newValue.ulong = atomicCAS(address, newValue.ulong, loc.ulong);
		
	return newValue.ulong;
}

struct Vertex {
	int oneRingID;
	int n_oneRing;
	int id;				// own vertex id
};

__device__ __forceinline__ float conditionQuad( const real_cu p[4][2] ) {
	real_cu e[4][2];
	real_cu e_length_squared[4];

	for( int i = 0; i < 4; ++i ) {
		int j = ( i + 1 ) % 4;
		e[i][0] = p[j][0] - p[i][0];
		e[i][1] = p[j][1] - p[i][1];

		e_length_squared[i] = e[i][0] * e[i][0] + e[i][1] * e[i][1];
	}
	
	float det0 = e[0][1] * e[3][0] - e[0][0] * e[3][1];
	float det1 = e[1][1] * e[0][0] - e[1][0] * e[0][1];
	float det3 = e[3][1] * e[2][0] - e[3][0] * e[2][1];

	float det = fminf( det0, det1 );
	det = fminf( det, det3 );

	if( det < 0 )
		return det;

	float c0 = 2 * det0 / ( e_length_squared[0] + e_length_squared[3] );
	float c1 = 2 * det1 / ( e_length_squared[1] + e_length_squared[0] );
	float c3 = 2 * det3 / ( e_length_squared[3] + e_length_squared[2] );

	float c = fminf( c0, c1 );
	c = fminf( c, c3 );

	return c;
}

__device__ __forceinline__ float computeConditionQuality( const int n_oneRing, const real_cu oneRing[MAX_ONE_RING_SIZE], const real_cu p[2] ) {
	float q = FLT_MAX;
	for( int k = 0; k < n_oneRing - 1; k += 2 ) {
		real_cu v[4][2] = { { p[0], p[1] },{ oneRing[2 * k], oneRing[2 * k + 1] },{ oneRing[2 * ( k + 1 )], oneRing[2 * ( k + 1 ) + 1] },{ oneRing[2 * ( k + 2 )], oneRing[2 * ( k + 2 ) + 1] } };
		q = fminf( q, conditionQuad( v ) );
	}
	return q;
}

__device__ __forceinline__ float computeLaplaceConditionQuality( const int n_oneRing, const real_cu oneRing[MAX_ONE_RING_SIZE], const real_cu p[2] ) {

	float q = FLT_MAX;

	// compute laplace point
	real_cu lp[2] = { 0,0 };
	for( int k = 0; k < n_oneRing - 1; ++k ) {
		lp[0] += oneRing[2 * k];
		lp[1] += oneRing[2 * k + 1];
	}
	lp[0] /= ( n_oneRing - 1 );
	lp[1] /= ( n_oneRing - 1 );
	lp[0] = p[0] - lp[0];
	lp[1] = p[1] - lp[1];

	for( int k = 0; k < n_oneRing - 1; k += 2 ) {
		real_cu v[4][2] = { { p[0], p[1] },{ oneRing[2 * k], oneRing[2 * k + 1] },{ oneRing[2 * ( k + 1 )], oneRing[2 * ( k + 1 ) + 1] },{ oneRing[2 * ( k + 2 )], oneRing[2 * ( k + 2 ) + 1] } };
		q = fminf( q, conditionQuad( v ) );
	}
	
	if( q < 0.5 )
		return q;
	else
		return 0.5f + 1.f / ( lp[0] * lp[0] + lp[1] * lp[1] + 1 );
}

template<int type = 0>
__device__ __forceinline__ float quality(const int n_oneRing, const real_cu oneRing[MAX_ONE_RING_SIZE], const real_cu p[2], const int q_crit) {
	if constexpr( type == 0 ) {
		return computeConditionQuality( n_oneRing, oneRing, p );
	} 
	if constexpr( type == 1 ) {
		return computeLaplaceConditionQuality( n_oneRing, oneRing, p );
	}
	return -1;
}

template<int type = 0>
__global__ void optimizeHierarchical(int* coloredVertexIDs, const int cOff, const Vertex* vertices, real_cu* vertexPos, int* oneRingVec, const real_cu affineFactor, const real_cu grid_scale) {
	const int i1 = threadIdx.x / NQ;
	const int j1 = threadIdx.x % NQ;

	const int i2 = (threadIdx.x + NQ * NQ / 2) / NQ;
	const int j2 = (threadIdx.x + NQ * NQ / 2) % NQ;

	const Vertex& v = vertices[coloredVertexIDs[cOff + blockIdx.x]];

	float q = -FLT_MAX;

	__shared__ real_cu xPos, yPos;
	__shared__ real_cu maxDistx;
	__shared__ real_cu maxDisty;

	__shared__ my_atomics argMaxVal;
	argMaxVal.floats[0] = -FLT_MAX;
	argMaxVal.ints[1] = NQ*NQ;

	__shared__ real_cu oneRing[MAX_ONE_RING_SIZE];

	// min/max search + loading oneRing
	if (threadIdx.x == 0) {
		maxDistx = -FLT_MAX;
		maxDisty = -FLT_MAX;

		for (int k = 0; k < v.n_oneRing - 1; ++k) {
			real_cu oneRingX = vertexPos[2 * oneRingVec[v.oneRingID + k]];
			real_cu oneRingY = vertexPos[2 * oneRingVec[v.oneRingID + k] + 1];
			oneRing[2 * k] = oneRingX;
			oneRing[2 * k + 1] = oneRingY;

			real_cu xDist = abs(vertexPos[2 * v.id] - oneRingX);
			real_cu yDist = abs(vertexPos[2 * v.id + 1] - oneRingY);

			maxDistx = fmaxf(maxDistx, xDist);
			maxDisty = fmaxf(maxDisty, yDist);
		}
		
		oneRing[2 * v.n_oneRing - 2] = vertexPos[2 * oneRingVec[v.oneRingID + v.n_oneRing - 1]];
		oneRing[2 * v.n_oneRing - 1] = vertexPos[2 * oneRingVec[v.oneRingID + v.n_oneRing - 1] + 1];

		xPos = vertexPos[2 * v.id];
		yPos = vertexPos[2 * v.id + 1];
	}
	__syncwarp();

	// start depth iteration
	real_cu depth_scale = grid_scale;
	real_cu argMax = 0;
	for (int depth = 0; depth < DMO_DEPTH; ++depth) {

		real_cu xMax, xMin, yMax, yMin;
		xMax = xPos + depth_scale * maxDistx;
		xMin = xPos - depth_scale * maxDistx;
		yMax = yPos + depth_scale * maxDisty;
		yMin = yPos - depth_scale * maxDisty;


		real_cu pos_i1 = affineFactor * (i1 * xMin + (NQ - 1 - i1) * xMax);
		real_cu pos_j1 = affineFactor * (j1 * yMin + (NQ - 1 - j1) * yMax);
		real_cu pos_i2 = affineFactor * (i2 * xMin + (NQ - 1 - i2) * xMax);
		real_cu pos_j2 = affineFactor * (j2 * yMin + (NQ - 1 - j2) * yMax);

		real_cu p1[2] = { pos_i1, pos_j1 };
		float q1 = quality<type>( v.n_oneRing, oneRing, p1, 8 );
		real_cu p2[2] = { pos_i2, pos_j2 };
		float q2 = quality<type>( v.n_oneRing, oneRing, p2, 8 );

		if (q1 > q2) {
			q = q1;
			argMax = 1;
		}
		else {
			q = q2;
			argMax = 2;
		}
		__syncwarp();
		my_atomicArgMax( (uint64_t *)&( argMaxVal.ulong ), q, i1 * NQ + j1 );

		real_cu pCurrent[2] = { xPos, yPos };
		float qOld = quality(v.n_oneRing, oneRing, pCurrent, 8);
		if (i1 * NQ + j1 == argMaxVal.ints[1] && qOld < q) {
			if (argMax == 1) {
				xPos = pos_i1;
				yPos = pos_j1;
			}
			else {
				xPos = pos_i2;
				yPos = pos_j2;
			}
		}
		
		//depth dependent scaling factor
		depth_scale = depth_scale * (2.f / (NQ - 1));
	}

	// set new position if it is better than the old one
	real_cu pOld[2] = { vertexPos[2 * v.id] , vertexPos[2 * v.id + 1] };
	float qOld = quality( v.n_oneRing, oneRing, pOld, 8 );
	if (i1 * NQ + j1 == argMaxVal.ints[1] && qOld < q) {
		vertexPos[2 * v.id] = xPos;
		vertexPos[2 * v.id + 1] = yPos;
	}
}

struct UniformGrid
{
	int nx, ny;
	real_cu hx, hy, xMin, yMin, xMax, yMax;
};

inline void copyOpenMeshData( PolyMesh& mesh, real_cu* vertexPos, Vertex* vertices, int* oneRingVec) {

	int interior_counter = 0;
	int oneRing_counter = 0;
	for( auto vh : mesh.vertices() ) {
		auto p = mesh.point( vh );

		vertexPos[2 * vh.idx()] = p[0];
		vertexPos[2 * vh.idx() + 1] = p[1];

		if( !vh.is_boundary() ) {
			// fill vertex struct

			Vertex& v = vertices[interior_counter];
			v.id = vh.idx();

			v.n_oneRing = vh.valence() * 2 + 1;

			v.oneRingID = oneRing_counter;

			LOG_ASSERT( v.n_oneRing <= MAX_ONE_RING_SIZE / 2 );

			auto heh = vh.out();
			auto heh_init = heh;

			do {
				oneRingVec[oneRing_counter++] = heh.to().idx();
				heh = heh.next();
				oneRingVec[oneRing_counter++] = heh.to().idx();
				heh = heh.next().next().opp();
			} while( heh != heh_init );

			oneRingVec[oneRing_counter] = heh.to().idx();
			++oneRing_counter;
			++interior_counter;
		}
	}
}

inline void createColoring( PolyMesh& mesh, const int n_free_vertices, int** coloredVertexIDs, std::vector<int>& colorOffset) {

	// create coloring scheme
	std::vector<int>colorScheme(mesh.n_vertices(), -1);
	int colorSchemeIt = 0;

	// set boundarys to a value that can be ignored
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		if (mesh.is_boundary(*v_it)) {
			colorScheme[v_it->idx()] = -2;
		}
	}

	while (std::find(colorScheme.begin(), colorScheme.end(), -1) != colorScheme.end()) {
		for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
	
			if (colorScheme[v_it->idx()] != -1) { continue; }		// vertex is already colored
	
			bool neighborIsCurrent = false;
			for (auto voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); ++voh_it) {
				PolyMesh::VertexHandle vh1 = mesh.to_vertex_handle(*voh_it);
				PolyMesh::VertexHandle vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(*voh_it));
				if (colorScheme[vh1.idx()] == colorSchemeIt || colorScheme[vh2.idx()] == colorSchemeIt) {
					neighborIsCurrent = true;
					break;
				}
			}
			if (neighborIsCurrent) { continue; }			// a neighboring vertex is already in this color
	
			colorScheme[v_it->idx()] = colorSchemeIt;
		}
		++colorSchemeIt;
	}

	int n_colors = *(std::max_element(colorScheme.begin(), colorScheme.end())) + 1;

	std::vector<int> n_color_vecs(n_colors, 0);
	for (int i = 0; i < colorScheme.size(); ++i) {
		if (colorScheme[i] > -1)
			++n_color_vecs[colorScheme[i]];
	}

	*coloredVertexIDs = new int[n_free_vertices];

	colorOffset = std::vector<int>(n_colors + 1, 0);
	for (int i = 1; i < n_colors; ++i) {
		colorOffset[i] = colorOffset[i - 1] + n_color_vecs[i - 1];
	}
	colorOffset[n_colors] = n_free_vertices;		// mark the end of the colored-vertices vector

													// add vertex ids
	std::vector<int>colorCounter(n_colors, 0);
	int interior_counter = 0;
	for (int i = 0; i < colorScheme.size(); ++i) {
		if (colorScheme[i] < 0) { continue; }
		(*coloredVertexIDs)[colorOffset[colorScheme[i]] + colorCounter[colorScheme[i]]++] = interior_counter++;
	}
}

template<int type = 0>
void discreteMeshOptimization( PolyMesh& mesh, const float grid_scale = 0.5f, int n_iter = 100) {
	
	int n_free_vertices = 0;
	int oneRingVecLength = 0;
#pragma omp parallel for reduction(+:n_free_vertices,oneRingVecLength)
	for (int i = 0; i < mesh.n_vertices(); ++i) {
		PolyMesh::VertexHandle vh = mesh.vertex_handle(i);
		if (mesh.is_boundary(vh)) { continue; }
		++n_free_vertices;

		for( auto voh : mesh.voh_range( vh ) ) {
			oneRingVecLength += 2;
		}
		++oneRingVecLength;		// additional count s.th. last element is again the first element
	}

	// convert OpenMesh to a basic structure
	real_cu* vertexPos = new real_cu[2 * mesh.n_vertices()];
	Vertex* vertices = new Vertex[n_free_vertices];
	int* oneRingVec = new int[oneRingVecLength];

	real_cu* vertexPos_d;
	Vertex* vertices_d;
	int* oneRingVec_d;
	int* coloredVertexIDs_d;

	int* coloredVertexIDs;
	std::vector<int> colorOffset;


#pragma omp parallel sections num_threads(2)
	{
#pragma omp section
		{
			gpuErrchk(cudaMalloc((void**)&vertexPos_d, 2 * mesh.n_vertices() * sizeof(real_cu)));
			gpuErrchk(cudaMalloc((void**)&vertices_d, n_free_vertices * sizeof(Vertex)));
			gpuErrchk(cudaMalloc((void**)&oneRingVec_d, oneRingVecLength * sizeof(int)));
			gpuErrchk(cudaMalloc((void**)&coloredVertexIDs_d, n_free_vertices * sizeof(int)));

			createColoring(mesh, n_free_vertices, &coloredVertexIDs, colorOffset);

			gpuErrchk(cudaMemcpyAsync(coloredVertexIDs_d, coloredVertexIDs, n_free_vertices * sizeof(int), cudaMemcpyHostToDevice));
		}
#pragma omp section 
		{
			copyOpenMeshData(mesh, vertexPos, vertices, oneRingVec);
		}
	}

	gpuErrchk(cudaMemcpyAsync(vertexPos_d, vertexPos, 2 * mesh.n_vertices() * sizeof(real_cu), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyAsync(vertices_d, vertices, n_free_vertices * sizeof(Vertex), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyAsync(oneRingVec_d, oneRingVec, oneRingVecLength * sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk( cudaDeviceSynchronize() );
	
	const int n_colors = (int)colorOffset.size() - 1;


	constexpr real_cu affineFactor = 1.f / (real_cu)(NQ - 1);

	for (int i = 0; i < n_iter; ++i) {
		for (int cid = 0; cid < n_colors; ++cid) {
			const int nBlocks = colorOffset[cid + 1] - colorOffset[cid];
			const int nThreads = NQ * NQ / 2;
			//std::cout << "i = " << i << "  |  color = " << cid << std::endl;
			optimizeHierarchical << <nBlocks, nThreads >> >(coloredVertexIDs_d, colorOffset[cid], vertices_d, vertexPos_d, oneRingVec_d, affineFactor, grid_scale);
			gpuErrchk( cudaDeviceSynchronize() );
		}
	}
	
	gpuErrchk( cudaDeviceSynchronize() );
	cudaMemcpy(vertexPos, vertexPos_d, 2 * mesh.n_vertices() * sizeof(real_cu), cudaMemcpyDeviceToHost);

	cudaFree(vertexPos_d);
	cudaFree(vertices_d);
	cudaFree(oneRingVec_d);
	cudaFree(coloredVertexIDs_d);

	delete[] vertices;
	delete[] oneRingVec;
	delete[] coloredVertexIDs;
	
	// write vertex positions back to mesh
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
		int id = v_it->idx();
		TriMesh::Point p = { vertexPos[2 * id], vertexPos[2 * id + 1], 0.f };
		mesh.set_point(*v_it, p);
	}

	delete[] vertexPos;
}


template void discreteMeshOptimization<0>( PolyMesh&, const float, int );
template void discreteMeshOptimization<1>( PolyMesh&, const float, int );