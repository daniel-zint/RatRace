#include "DualHybridGraph.h"

#include <iostream>
#include <queue>

#include "QualityMetrics.h"


DualHybridGraph::DualHybridGraph( PolySurfMesh& mesh ) : mesh_( mesh ), isTriangleVec_(mesh.n_faces(), false) {
	halfedges.reserve( mesh_.n_halfedges() * 2 );
	nodes.reserve( mesh_.n_faces() );
	// get nodes
	for ( auto fh : mesh_.faces() ) {
		auto& n = nodes.add();
	}

	for ( auto heh : mesh_.halfedges() ) {
		auto n0 = mesh_.face_handle( heh ).idx();
		auto n1 = mesh_.opposite_face_handle( heh ).idx();
		auto v = mesh_.to_vertex_handle( heh );
		
		auto& graphHalfedge = halfedges.add();
		graphHalfedge.from() = n0;
		graphHalfedge.vertex() = v;

		graphHalfedge.to() = n1;
		graphHalfedge.oppositeIdx() = mesh_.opposite_halfedge_handle( heh ).idx();
	}
	
	// get next and prev halfedge
	for ( auto fh : mesh_.faces() ) {
		for ( auto fhh : mesh_.fh_range( fh ) ) {
			auto& gheh = halfedges[fhh.idx()];

			auto heh = mesh_.opposite_halfedge_handle( mesh_.next_halfedge_handle( fhh ) );
			gheh.prevIdx() = heh.idx();
		}
	}

	for ( auto& gheh : halfedges ) {
		if ( gheh.prevIdx() == -1 ) continue;
		auto& ghehPrev = halfedges[gheh.prevIdx()];
		ghehPrev.nextIdx() = gheh.idx();
	}

	// assign each node a halfedge
	for ( auto fh : mesh_.faces() ) {
		nodes[fh.idx()].halfedgeIdx() = mesh_.halfedge_handle( fh ).idx();
	}

	for( auto& n : nodes ) {
		if( isTriangle( n.idx() ) ) {
			isTriangleVec_[n.idx()] = true;
		}
	}
}

std::vector<DualHybridGraph::NodeIdxT> DualHybridGraph::neighborIdxs( const NodeIdxT & idx ) {

	LOG_ASSERT( idx >= 0 );
	LOG_ASSERT( idx < nodes.size() );

	std::vector<NodeIdxT> neighs;

	auto gheh = halfedges[nodes[idx].halfedgeIdx()];
	auto ghehInit = gheh;
	do {
		neighs.push_back( gheh.to() );
		gheh = halfedges[gheh.prevIdx()];
		gheh = halfedges[gheh.oppositeIdx()];
	} while ( gheh.idx() != ghehInit.idx() );


	return neighs;
}

std::vector<DualHybridGraph::NodeIdxT> DualHybridGraph::innerNeighborIdxs( const NodeIdxT & idx ) {

	LOG_ASSERT( idx >= 0 );
	LOG_ASSERT( idx < nodes.size() );

	std::vector<NodeIdxT> neighs;

	auto gheh = halfedges[nodes[idx].halfedgeIdx()];
	auto ghehInit = gheh;
	do {
		if(gheh.to() != -1 )
			neighs.push_back( gheh.to() );
		gheh = halfedges[gheh.prevIdx()];
		gheh = halfedges[gheh.oppositeIdx()];
	} while ( gheh.idx() != ghehInit.idx() );


	return neighs;
}

auto DualHybridGraph::outgoingHalfedges( const NodeIdxT & idx ) {
	return OutgoingHalfedgeRange(idx, *this);
}

DualHybridGraph::HalfedgeIdxT DualHybridGraph::findHalfedge( const NodeIdxT& nodeFrom, const NodeIdxT& nodeTo ) {

	LOG_ASSERT( nodeFrom >= 0 );
	LOG_ASSERT( nodeTo >= 0 );
	LOG_ASSERT( nodeFrom < nodes.size() );
	LOG_ASSERT( nodeTo < nodes.size() );

	DualHybridGraph::HalfedgeIdxT heIdx = -1;
	for( auto o : outgoingHalfedges( nodeFrom ) ) {
		if( o.to() == nodeTo ) {
			heIdx = o.idx();
			break;
		}
	}

	return heIdx;
}

PolyMesh DualHybridGraph::getMesh() {
	PolyMesh m;

	for ( auto vh : mesh_.vertices() ) {
		m.add_vertex( mesh_.point( vh ) );
	}

	for ( auto n : nodes ) {
		if ( !n.isValid() ) continue;
		std::vector<PolyMesh::VertexHandle> vhs;
		for( auto oh : outgoingHalfedges( n.idx() ) ) {
			vhs.push_back( oh.vertex() );
		}
		m.add_face( vhs );
	}

	return m;
}

void DualHybridGraph::printMesh( const std::experimental::filesystem::path & filename ) {
	namespace fs = std::experimental::filesystem;

	LOG_ASSERT( fs::exists( filename.parent_path() ) );

	auto m = getMesh();

	OpenMesh::IO::write_mesh( m, filename.string() );
}

DualHybridGraph::NodeIdxT DualHybridGraph::halfedgeCollapse( const HalfedgeIdxT & idx ) {

	LOG_ASSERT( idx >= 0 );
	LOG_ASSERT( idx < halfedges.size() );

	// collect multi-edges
	std::vector<HalfedgeIdxT> multiIdxs;
	auto idxIter = idx;
	while( halfedges[idxIter].nextIdx() == halfedges[idxIter].prevIdx() ) {
		idxIter = halfedges[idxIter].nextIdx();
		idxIter = halfedges[idxIter].oppositeIdx();
	}
	auto idxReal = idxIter;
	multiIdxs.push_back( idxReal );
	idxIter = halfedges[idxReal].oppositeIdx();
	multiIdxs.push_back( idxIter );
	while( halfedges[idxIter].nextIdx() == halfedges[idxIter].prevIdx() ) {
		idxIter = halfedges[idxIter].nextIdx();
		multiIdxs.push_back( idxIter );
		idxIter = halfedges[idxIter].oppositeIdx();
		multiIdxs.push_back( idxIter );
	}
	auto idxOppReal = idxIter;

	auto& h = halfedges[idxReal];
	auto& hOpp = halfedges[idxOppReal];
	auto& nodeFrom = nodes[h.from()];
	auto& nodeTo = nodes[h.to()];

	// merge multi-edges
	if( multiIdxs.size() > 2 ) {
		for( auto i = 1; i < multiIdxs.size() - 1; i += 2 ) {
			auto vh = halfedges[multiIdxs[i]].vertex();
			nodeTo.vertices().push_back( vh );
		}

		for( auto i = 1; i < multiIdxs.size() - 1; ++i ) {
			halfedges[multiIdxs[i]].isValid() = false;
		}

		h.oppositeIdx() = hOpp.idx();
		hOpp.oppositeIdx() = h.idx();

		// make sure that the vertex holds a valid halfedge idx
		nodeFrom.halfedgeIdx() = h.idx();
		nodeTo.halfedgeIdx() = hOpp.idx();
	}


	auto& hNext = halfedges[h.nextIdx()];
	auto& hPrev = halfedges[h.prevIdx()];

	auto& hOppNext = halfedges[hOpp.nextIdx()];
	auto& hOppPrev = halfedges[hOpp.prevIdx()];
	

	for( auto& oh : outgoingHalfedges( nodeFrom.idx() ) ) {
		oh.from() = nodeTo.idx();
		auto& ih = halfedges[oh.oppositeIdx()];
		ih.to() = nodeTo.idx();
	}

	hNext.prevIdx() = h.prevIdx();
	hPrev.nextIdx() = h.nextIdx();
	
	hOppNext.prevIdx() = hOpp.prevIdx();
	hOppPrev.nextIdx() = hOpp.nextIdx();

	h.isValid() = false;
	hOpp.isValid() = false;

	nodeFrom.isValid() = false;

	// give nodeTo a valid halfedge (it currently has hOpp as halfedge which is now invalid)
	nodeTo.halfedgeIdx() = hNext.idx();

	return nodeFrom.idx();
}

void DualHybridGraph::split3to4( const HalfedgeIdxT & idx, const NodeIdxT& invalidNodeIdx ) {

	LOG_ASSERT( idx >= 0 );
	LOG_ASSERT( idx < halfedges.size() );

	// collect multi-edges
	std::vector<HalfedgeIdxT> multiIdxs;
	auto idxIter = idx;
	while ( halfedges[idxIter].nextIdx() == halfedges[idxIter].prevIdx() ) {
		idxIter = halfedges[idxIter].nextIdx();
		idxIter = halfedges[idxIter].oppositeIdx();
	}
	auto idxReal = idxIter;
	multiIdxs.push_back( idxReal );
	idxIter = halfedges[idxReal].oppositeIdx();
	multiIdxs.push_back( idxIter );
	while ( halfedges[idxIter].nextIdx() == halfedges[idxIter].prevIdx() ) {
		idxIter = halfedges[idxIter].nextIdx();
		multiIdxs.push_back( idxIter );
		idxIter = halfedges[idxIter].oppositeIdx();
		multiIdxs.push_back( idxIter );
	}
	auto idxOppReal = idxIter;
	LOG_ASSERT( multiIdxs.size() <= 4 );

	const auto h = halfedges[idxReal];
	LOG_ASSERT( h.isValid() );

	const auto nodeTri = nodes[h.from()];
	LOG_ASSERT( nodeTri.vertices().size() < 2 ) << "Cannot handle nodes with more than one inner vertex";

	idxIter = idxReal;

	std::vector<HalfedgeIdxT> ohs;	// nodeFrom outgoing halfedges
	std::vector<HalfedgeIdxT> ihs;	// nodeFrom ingoing halfedges		// Not used so far
	do {
		ohs.push_back( halfedges[idxIter].idx() );
		ihs.push_back( halfedges[idxIter].oppositeIdx() );
		idxIter = halfedges[idxIter].prevIdx();
		idxIter = halfedges[idxIter].oppositeIdx();
	} while ( idxIter != idxReal );

	std::vector<PolyMesh::VertexHandle> vhs;
	vhs.reserve( ohs.size() );
	for ( auto& oh : ohs ) {
		vhs.push_back( halfedges[oh].vertex() );
	}

	// add halfedges and node
	auto hNewIdx = halfedges.add().idx();
	auto hNewOppIdx = halfedges.add().idx();
	//DualHybridGraph::Node& nodeQuad = nodes.add();
	auto& nodeQuad = nodes[invalidNodeIdx];
	nodeQuad = Node( invalidNodeIdx );
	auto& hNew = halfedges[hNewIdx];
	auto& hNewOpp = halfedges[hNewOppIdx];
	// set from/to and opposite
	hNew.from() = nodeQuad.idx();
	hNew.to() = nodeTri.idx();
	hNew.oppositeIdx() = hNewOpp.idx();
	hNewOpp.from() = nodeTri.idx();
	hNewOpp.to() = nodeQuad.idx();
	hNewOpp.oppositeIdx() = hNew.idx();
	// update halfedges in nodes
	nodeQuad.halfedgeIdx() = hNew.idx();
	nodes[nodeTri.idx()].halfedgeIdx() = hNewOpp.idx();


	std::vector<std::vector<PolyMesh::VertexHandle>> quadCases, triCases;
	if ( nodeTri.vertices().size() == 1 ) {
		// valence 2 special case
		LOG_ASSERT( ohs.size() == 3 ) << "Cannot perform 3/4 split on this node";
		for ( auto i = 0; i < multiIdxs.size(); i += 2 ) {
			auto b = halfedges[multiIdxs[i]].vertex();
			auto a = halfedges[multiIdxs[i + 1]].vertex();
			auto c = nodeTri.vertices()[0];
			triCases.push_back( { a,b,c } );
		}
	} else {
		// standard case
		LOG_ASSERT( ohs.size() == 5 ) << "Cannot perform 3/4 split on this node";
		auto a = halfedges[multiIdxs[1]].vertex();
		auto b = halfedges[multiIdxs[0]].vertex();
		auto c = halfedges[halfedges[halfedges[multiIdxs[0]].prevIdx()].oppositeIdx()].vertex();
		triCases.push_back( { a,b,c } );

		if ( multiIdxs.size() == 4 ) {
			auto a = halfedges[multiIdxs[3]].vertex();
			auto b = halfedges[multiIdxs[2]].vertex();
			auto c = halfedges[halfedges[halfedges[multiIdxs[3]].nextIdx()].oppositeIdx()].vertex();
			triCases.push_back( { a,b,c } );
		} else {
			auto c = halfedges[halfedges[halfedges[multiIdxs[1]].nextIdx()].oppositeIdx()].vertex();
			triCases.push_back( { a,b,c } );
		}
	}

	// get halfedges for the corresponding points
	std::vector<std::vector<HalfedgeIdxT>> quadOh;
	for ( const auto& triVh : triCases ) {
		std::vector<HalfedgeIdxT> quad;
		for ( auto& oh : ohs ) {
			auto b = halfedges[oh].vertex();
			auto a = halfedges[halfedges[oh].oppositeIdx()].vertex();
			if ( std::find( triVh.begin(), triVh.end(), a ) == triVh.end() || std::find( triVh.begin(), triVh.end(), b ) == triVh.end() ) {
				quad.push_back( oh );
			}
		}
		quadOh.push_back( quad );
	}

	// get quadVhs
	for ( const auto& oh : quadOh ) {
		std::vector<PolyMesh::VertexHandle> quadVh;
		quadVh.push_back( halfedges[halfedges[oh[0]].oppositeIdx()].vertex() );
		quadVh.push_back( halfedges[oh[0]].vertex() );

		for ( auto i = 1; i < oh.size(); ++i ) {
			quadVh.push_back( halfedges[oh[i]].vertex() );
		}

		if ( nodeTri.vertices().size() == 1 ) {
			quadVh.push_back( nodeTri.vertices()[0] );
		}

		quadCases.push_back( quadVh );
	}

	// get best quadCase
	int bestCase = -1;
	real bestQuality = -FLT_MAX;
	for ( auto i = 0; i < quadCases.size(); ++i ) {
		auto& quad = quadCases[i];
		LOG_ASSERT( quad.size() == 4 );

		std::vector<PolyMesh::Point> p;
		p.reserve( 4 );
		for ( auto& vh : quad ) {
			p.push_back( mesh_.point( vh ) );
		}
		auto q = QualityMetrics::conditionSurf( p );

		if ( q > bestQuality ) {
			bestQuality = q;
			bestCase = i;
		}
	}
	LOG_ASSERT( bestCase != -1 );

	auto quadNewOhs = quadOh[bestCase];
	quadNewOhs.push_back( hNew.idx() );

	// add new halfedge to tri halfedges (the correct position in the vector is important!
	bool hNewOppInserted = false;
	decltype( quadNewOhs ) triNewOhs;
	triNewOhs.reserve( 3 );
	for ( const auto& oh : ohs ) {
		if ( std::find( quadNewOhs.begin(), quadNewOhs.end(), oh) != quadNewOhs.end() ) {
			if ( !hNewOppInserted ) {
				triNewOhs.push_back( hNewOpp.idx() );
				hNewOppInserted = true;
			}
		} else {
			triNewOhs.push_back( oh );
		}
	}

	// perform split
	for ( auto& h : quadNewOhs ) {
		halfedges[h].from() = nodeQuad.idx();
		halfedges[halfedges[h].oppositeIdx()].to() = nodeQuad.idx();
	}

	for ( auto& h : triNewOhs ) {
		halfedges[h].from() = nodeTri.idx();
		halfedges[halfedges[h].oppositeIdx()].to() = nodeTri.idx();
	}

	for ( auto i = 0; i < quadNewOhs.size(); ++i ) {
		auto nextOh = ( i + 1 ) % quadNewOhs.size();
		auto prevOh = ( i - 1 + quadNewOhs.size() ) % quadNewOhs.size();

		halfedges[quadNewOhs[i]].prevIdx() = halfedges[quadNewOhs[nextOh]].oppositeIdx();
		halfedges[halfedges[quadNewOhs[i]].oppositeIdx()].nextIdx() = halfedges[quadNewOhs[prevOh]].idx();
	}
	for ( auto i = 0; i < triNewOhs.size(); ++i ) {
		auto nextOh = ( i + 1 ) % triNewOhs.size();
		auto prevOh = ( i - 1 + triNewOhs.size() ) % triNewOhs.size();

		halfedges[triNewOhs[i]].prevIdx() = halfedges[triNewOhs[nextOh]].oppositeIdx();
		halfedges[halfedges[triNewOhs[i]].oppositeIdx()].nextIdx() = halfedges[triNewOhs[prevOh]].idx();
	}
	// add vertices to new halfedges
	hNew.vertex() = halfedges[hNew.prevIdx()].vertex();
	hNewOpp.vertex() = halfedges[hNewOpp.prevIdx()].vertex();
		
	// insert multiple edges
	if ( nodeTri.vertices().size() == 1 ) {
		auto hMultIdx = halfedges.add().idx();
		auto hMultOppIdx = halfedges.add().idx();
		auto& hMult = halfedges[hMultIdx];
		auto& hMultOpp = halfedges[hMultOppIdx];
		hMult.nextIdx() = hMultOpp.idx();
		hMult.prevIdx() = hMultOpp.idx();
		hMultOpp.nextIdx() = hMult.idx();
		hMultOpp.prevIdx() = hMult.idx();
		hMult.from() = halfedges[hNewIdx].from();
		hMult.to() = halfedges[hNewIdx].to();
		hMultOpp.from() = halfedges[hNewOppIdx].from();
		hMultOpp.to() = halfedges[hNewOppIdx].to();
		hMult.oppositeIdx() = hNewOppIdx;
		hMultOpp.oppositeIdx() = hNewIdx;
		hMult.vertex() = nodes[h.from()].vertices()[0];
		hMultOpp.vertex() = nodes[h.from()].vertices()[0];
		nodes[h.from()].vertices().clear();

		halfedges[hNewIdx].oppositeIdx() = hMultOppIdx;
		halfedges[hNewOppIdx].oppositeIdx() = hMultIdx;
	}

}

void DualHybridGraph::garbageCollection() {

	// eliminate halfedges
	std::vector<HalfedgeIdxT> newHalfedgeIdxs( halfedges.size(), -1 );
	decltype( halfedges ) newHalfedges;
	for( auto& h : halfedges ) {
		if( !h.isValid() ) continue;
		Halfedge hNew( newHalfedges.size() );
		hNew.from() = h.from();
		hNew.to() = h.to();
		hNew.nextIdx() = h.nextIdx();
		hNew.oppositeIdx() = h.oppositeIdx();
		hNew.prevIdx() = h.prevIdx();
		hNew.vertex() = h.vertex();
		newHalfedges.push_back( hNew );
		newHalfedgeIdxs[h.idx()] = hNew.idx();
	}
	halfedges = newHalfedges;

	// update connections in halfedges
	for( auto& h : halfedges ) {
		if(h.nextIdx() != -1)
			h.nextIdx() = newHalfedgeIdxs[h.nextIdx()];
		if( h.prevIdx() != -1 )
			h.prevIdx() = newHalfedgeIdxs[h.prevIdx()];
		if( h.oppositeIdx() != -1 )
			h.oppositeIdx() = newHalfedgeIdxs[h.oppositeIdx()];
	}

	// update halfedge-idxs in nodes
	for( auto& n : nodes ) {
		n.halfedgeIdx() = newHalfedgeIdxs[n.halfedgeIdx()];
	}

	// eliminate nodes
	std::vector<HalfedgeIdxT> newNodeIdxs( nodes.size(), -1 );
	decltype( nodes ) newNodes;
	for( auto& n : nodes ) {
		if( !n.isValid() ) continue;
		Node nNew( newNodes.size() );
		nNew.halfedgeIdx() = n.halfedgeIdx();
		nNew.vertices() = n.vertices();
		newNodes.push_back( nNew );
		newNodeIdxs[n.idx()] = nNew.idx();
	}
	nodes = newNodes;

	// update node-idxs in halfedges
	for( auto& h : halfedges ) {
		if(h.to() != -1)
			h.to() = newNodeIdxs[h.to()];
		if( h.from() != -1 )
			h.from() = newNodeIdxs[h.from()];
	}
}

// TODO refactor this function
bool DualHybridGraph::isTriangle( const NodeIdxT & idx ) {
	int i = 0;
	for( const auto& oh : outgoingHalfedges(idx) ) {
		++i;
	}
	return i == 3;
}

bool DualHybridGraph::hasTriangles() {
	for ( auto& n : nodes ) {
		if ( isTriangle( n.idx() ) ) return true;
	}
	return false;
}

int DualHybridGraph::nTriangles() {
	int count = 0;
	for( const auto& b : isTriangleVec_ ) {
		if( b )
			++count;
	}
	return count;
}

void DualHybridGraph::connectTriangles() {

	//auto nTriangles = this->nTriangles();
	//LOG_ASSERT( nTriangles % 2 == 0 );

	const auto nNodes = nodes.size();
	std::vector<DualHybridGraph::NodeIdxT> root( nNodes, -1 );

	for( int nBeginIdx = nNodes - 1; nBeginIdx >= 0; --nBeginIdx ) {
		if( !isTriangleVec_[nBeginIdx] )
			continue;

		std::fill( root.begin(), root.end(), -1 );
		//LOG( INFO ) << nTriangles;
		
		LOG_ASSERT( nBeginIdx != -1 );

		NodeIdxT nEndIdx = -1;

		std::queue<long> q;
		q.push( nBeginIdx );

		while( !q.empty() ) {
			auto idx = q.front();
			q.pop();
			if( isTriangleVec_[idx] && idx != nBeginIdx ) {
				nEndIdx = idx;
				break;
			}

			//for ( auto neighIdx : innerNeighborIdxs(idx) ) {
			//	if ( root[neighIdx] != -1 ) continue;
			//	q.push( neighIdx );
			//	root[neighIdx] = idx;
			//}

			// add neighbors to queue
			auto gheh = halfedges[nodes[idx].halfedgeIdx()];
			const auto ghehInit = gheh;
			do {
				//if( gheh.to() != -1 )
				//	neighs.push_back( gheh.to() );
				const auto& neighIdx = gheh.to();
				if( neighIdx != -1 && root[neighIdx] == -1 ) {
					q.push( neighIdx );
					root[neighIdx] = idx;
				}
				gheh = halfedges[gheh.prevIdx()];
				gheh = halfedges[gheh.oppositeIdx()];
			} while( gheh.idx() != ghehInit.idx() );
		}

		LOG_ASSERT( nEndIdx != -1 );

		// use backtracking to construct path
		std::vector<NodeIdxT> path;
		auto pathIter = nEndIdx;
		while( pathIter != nBeginIdx ) {
			path.push_back( pathIter );
			pathIter = root[pathIter];
		}
		path.push_back( nBeginIdx );

		// perform rat-race
		for( auto i = 0; i < path.size() - 2; ++i ) {
			auto h = findHalfedge( path[i], path[i + 1] );
			auto idx = halfedgeCollapse( h );

			auto hNext = findHalfedge( path[i + 1], path[i + 2] );
			split3to4( hNext, idx );
			//split3to4( hNext, nodes.add().idx() );
		}
		halfedgeCollapse( findHalfedge( path[path.size() - 2], path[path.size() - 1] ) );

		isTriangleVec_[nBeginIdx] = false;
		isTriangleVec_[nEndIdx] = false;

		//nTriangles = dhg.nTriangles();
		//if( nTriangles == nTrianglesOld ) {
		//	LOG( ERROR ) << "Cannot reduce number of triangles. nTriangles = " << nTriangles;
		//} else {
		//	nTrianglesOld = nTriangles;
		//}
	}
	this->garbageCollection();

	
}

void DualHybridGraph::printGridInfo() {
	std::cout << "Halfedges:" << std::endl;
	for ( const auto& h : halfedges ) {
		if ( !h.isValid() ) continue;
		std::cout << h.idx() << "    " << h.from() << " --> " << h.to() << "   |  vertex " << h.vertex().idx() << std::endl;
	}
}
