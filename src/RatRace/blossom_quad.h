#pragma once

#include "graph_classes.h"
#include <stdio.h>
#include <tchar.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <glog/logging.h>

#include "MeshHeader.h"

#include "../Blossom5\PerfectMatching.h"

constexpr int EXTERNAL_EDGE_COST = 10000;

class blossom_quad {
	PolyMesh workMesh;
	PolyMesh blossomMesh_step1;
	std::vector<PolyMesh::VertexHandle> blossomMesh_step1_vertex_vec;

public:
	PolyMesh do_blossom_algo(PolyMesh loadedMesh);

	/// <summary>
	/// Perform Blossom-Quad on boundaries
	/// </summary>
	/// <param name="mesh">triangle mesh</param>
	/// <returns>hybrid mesh with quads on the boundary almost everywhere</returns>
	PolyMesh quadrangleBoundaries( TriMesh& mesh );

private:
	
	/// <summary>
	/// connect two triangles to a quad
	/// </summary>
	/// <param name="heh">smart halfedge handle of the edge that should be removed</param>
	/// <returns></returns>
	std::vector<PolyMesh::VertexHandle> connect_triangles( const OpenMesh::SmartHalfedgeHandle& heh );
};