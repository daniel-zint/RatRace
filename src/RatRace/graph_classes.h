#pragma once
#include <stdio.h>
#include <vector>
#include "MeshHeader.h"

/// <summary>
/// Edge within the connectivity graph
/// </summary>
struct graph_edge {
	
	PolyMesh::FaceHandle face0;
	PolyMesh::FaceHandle face1;

	std::vector<PolyMesh::VertexHandle> quad_vertices;
	PolyMesh::VertexHandle connecting_vertex;

	bool to_be_used = false;
	int cost;
};

/// <summary>
/// Connectivity graph for Blossom-Quad
/// </summary>
struct con_graph {
	
	std::vector<graph_edge> edges;
	std::vector<graph_edge> external_edges;
	
	void print_graph() {
		for( int i = 0; i < edges.size(); i++ ) {
			std::cout << "edge between face " << edges[i].face0.idx() << " and face " << edges[i].face1.idx() << std::endl;
		}
		for( int j = 0; j < external_edges.size(); j++ ) {
			std::cout << "external edge between face " << external_edges[j].face0.idx() << " and face " << external_edges[j].face1.idx() << " over vertex " << external_edges[j].connecting_vertex.idx() << std::endl;
		}
	}
};

