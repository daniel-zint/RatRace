#pragma once

#include "MeshHeader.h"
#include <glog/logging.h>

class PolySurfMesh : public PolyMesh
{
	
public:
	PolySurfMesh& operator=( PolyMesh& m ) {
		clear();
		for( auto vh : m.vertices() ) {
			add_vertex( m.point( vh ) );
		}
		for( auto fh : m.faces() ) {
			std::vector<VertexHandle>vertices;
			for( auto fv : m.fv_range( fh ) ) {
				vertices.push_back( fv );
			}
			add_face( vertices );
		}
		return *this;
	}

	PolySurfMesh& operator=( TriMesh& m ) {
		clear();
		for( auto vh : m.vertices() ) {
			add_vertex( m.point( vh ) );
		}
		for( auto fh : m.faces() ) {
			std::vector<VertexHandle>vertices;
			for( auto fv : m.fv_range( fh ) ) {
				vertices.push_back( fv );
			}
			add_face( vertices );
		}
		return *this;
	}
	
	void edgeSwapNext( const EdgeHandle& eh );

	void edgeSwapPrev( const EdgeHandle& eh );

	/*heh is the halfedge pointing in the direction of the split. The split vertex is heh.from*/
	OpenMesh::SmartVertexHandle vertexSplit( const HalfedgeHandle& heh );
	/*heh is the splitted halfedge. The edge is splitted in three segments*/
	void edge3Split( const HalfedgeHandle& heh );

	/*collapse the face to which heh belongs by moving heh.from to the diagonal vertex*/
	void diagonalCollapse( const HalfedgeHandle& heh );

	void diagonalCollapseDirect( const HalfedgeHandle& heh );

	/*collapse halfedge for quad mesh. First perform vertex rotate around vh_from and then diagonal collapse */
	void halfedgeCollapse( const HalfedgeHandle& heh );

	void vertexRotate( const VertexHandle& vh );

};


PolyMesh mergeDoubleVertices( PolyMesh& pm );