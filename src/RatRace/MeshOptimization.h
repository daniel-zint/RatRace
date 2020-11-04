#pragma once

#include "MeshHeader.h"
#include "PolySurfMesh.h"

namespace MeshOptimization
{
	void getOneRing( PolyMesh& mesh, const TriMesh::VertexHandle& vh, std::vector<TriMesh::VertexHandle>& oneRingVh );

	real jacobianQuadMetric( TriMesh::Point points[4] );

	real minQuality( PolyMesh& mesh, const TriMesh::VertexHandle& vh, const std::vector<TriMesh::VertexHandle>& oneRingVh, const TriMesh::Point& p );

	void optimizeHierarchicalCPU( PolyMesh& mesh, const TriMesh::VertexHandle& vh );

	void discreteMeshOptimizationCPU( PolyMesh& mesh, const int n_iter );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	void laplaceSmoothing( PolySurfMesh& psm, const int n_iter );


	PolyMesh removeInnerValence2Vertices( PolyMesh& mesh );

	/*Do not use!!! Not working properly!!!*/
	PolyMesh removeBoundaryValence2Vertices( PolyMesh& mesh, real qMin = 0.02 );
}