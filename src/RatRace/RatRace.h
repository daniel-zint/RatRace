#pragma once

#include <glog/logging.h>
#include <experimental/filesystem>

#include "MeshHeader.h"
#include "DualHybridGraph.h"
#include "PolySurfMesh.h"
#include "../BackgroundGrid/SizeGrid.h"

class RatRace
{
	PolySurfMesh mesh_;
	real featureAngle_ = 30;

	OpenMesh::VPropHandleT<bool> vRemesh;
	OpenMesh::VPropHandleT<bool> vRemeshNew;
	OpenMesh::VPropHandleT<int> vValence;

public:
	void initBlossomQuad( TriMesh& tmesh );
	void initTriangleMerging( TriMesh& tmesh );
	void initHybrid( TriMesh& tmesh );
	void run();
	void postProcessing(const int& nQuads, SizeGrid* sizegrid = nullptr, const std::experimental::filesystem::path& outputPath = std::experimental::filesystem::path() );
	auto& mesh() { return mesh_; }

private:
	int swapValence();
	int removeInteriorValence2Vertices();
	int valence2SplitBoundary();
	int lowQualSplitBoundary();
	int removeBoundaryValence2Vertices();
	int collapse3x3y( SizeGrid* sizegrid = nullptr );
	int valenceSplit( int nQuadsMax = INT_MAX, SizeGrid* sizegrid = nullptr );

	void initRemeshProps();
	void swapRemeshProps();
	void deleteRemeshProps();

	int boundaryValence( const OpenMesh::SmartVertexHandle& vh );
};