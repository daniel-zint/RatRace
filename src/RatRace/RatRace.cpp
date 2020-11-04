#include "RatRace.h"

#include "MeshOptimization.h"
#include "QualityMetrics.h"
#include "blossom_quad.h"

#include <queue>
#include <array>

template<int = 0> void discreteMeshOptimization( PolyMesh& mesh, const float grid_scale = 0.5f, int n_iter = 100 );

/// <summary>
/// Compute the angle of the boundary at a given vertex
/// </summary>
/// <param name="mesh">quad mesh</param>
/// <param name="vh">vertex handle</param>
/// <returns>boundary angle in degrees</returns>
inline real calcBoundaryAngle( PolySurfMesh& mesh, const PolySurfMesh::VertexHandle& vh ) {
	LOG_ASSERT( mesh.is_boundary( vh ) );
	PolySurfMesh::HalfedgeHandle hehBound;
	for( auto voh : mesh.voh_range( vh ) ) {
		if( mesh.is_boundary( voh ) ) {
			hehBound = voh;
			break;
		}
	}
	LOG_ASSERT( hehBound.is_valid() );
	auto p = mesh.point( vh );
	auto b = mesh.point( mesh.to_vertex_handle( hehBound ) ) - p;
	b.normalize();
	auto a = mesh.point( mesh.from_vertex_handle( mesh.prev_halfedge_handle( hehBound ) ) ) - p;
	a.normalize();
	PolyMesh::Point n = { 0,0,1 };

	auto dot = a | b;
	auto det = ( a % b ) | n;
	auto angle = std::atan2( det, dot ) * 180.f / M_PI;
	if( angle < 0 )
		angle = 360 + angle;

	return angle;
}

/// <summary>
/// Triangle merging as described by Rank et al. in "Adaptive mesh generation and transformation of triangular to quadrilateral meshes".
/// </summary>
/// <param name="pmesh">triangle mesh</param>
inline void tri2hybridSurf( PolySurfMesh& pmesh ) {
	
	struct prioEdge
	{
		real w;
		OpenMesh::SmartEdgeHandle eh;
	};

	// map mit weight & edge
	std::vector<prioEdge> edges;
	edges.reserve( pmesh.n_edges() );
	std::vector<std::vector<PolyMesh::VertexHandle>> quadVertices( pmesh.n_edges() );
		
	for( const auto& eh : pmesh.edges() ) {

		if ( eh.is_boundary() )
			continue;

		const auto& heh1 = eh.h0();
		const auto& heh2 = eh.h1();

		if( heh1.face().valence() == 4 || heh2.face().valence() == 4 )
			continue;

		/*
		1-----0
		|   / |
		| /   |
		2-----3
		*/

		std::vector<PolyMesh::VertexHandle>vVec( 4 );
		vVec[0] = heh1.to();
		vVec[1] = heh1.next().to();
		vVec[2] = heh2.to();
		vVec[3] = heh2.next().to();

		std::vector<PolyMesh::Point>pVec( 4 );
		pVec[0] = pmesh.point( vVec[0] );
		pVec[1] = pmesh.point( vVec[1] );
		pVec[2] = pmesh.point( vVec[2] );
		pVec[3] = pmesh.point( vVec[3] );

		auto weight = 1.f - QualityMetrics::conditionSurf( pVec );
		
		//prioQueue.insert( { weight, eh } );
		edges.push_back( { weight, eh } );
		quadVertices[eh.idx()] = vVec;
	}


	// delete edges (tracke welche faces zusammengefasst wurden)

	std::vector<bool>faceDeleted( pmesh.n_faces(), false );

	pmesh.request_face_status();
	pmesh.request_edge_status();
	pmesh.request_vertex_status();

	std::vector<PolyMesh::EdgeHandle> edgesToDelete;
	edgesToDelete.reserve( pmesh.n_edges() / 2 );
		
	std::sort( edges.begin(), edges.end(), []( const prioEdge& a, const prioEdge& b ) {return a.w < b.w; } );

	for( const auto& [w, eh] : edges ) {
		const auto& fh1 = eh.h0().face();
		const auto& fh2 = eh.h1().face();

		// check if one of the faces is already part of a quad
		if( faceDeleted[fh1.idx()] || faceDeleted[fh2.idx()] )
			continue;

		faceDeleted[fh1.idx()] = true;
		faceDeleted[fh2.idx()] = true;

		edgesToDelete.push_back( eh );
	}

	for ( const auto& e : edgesToDelete ) {
		pmesh.delete_edge( e, false );
		pmesh.add_face( quadVertices[e.idx()] );
	}

	pmesh.garbage_collection();

}

/// <summary>
/// More efficient triangle merging which uses a breadth first search for generating a hybrid mesh with less remaining triangles.
/// </summary>
/// <param name="pmesh"></param>
inline void tri2hybridBFS( PolySurfMesh& pmesh ) {

	struct prioEdge
	{
		real w;
		OpenMesh::SmartEdgeHandle eh;
		int rank = -1;
	};

	// map mit weight & edge
	std::vector<prioEdge> edges;
	edges.reserve( pmesh.n_edges() );
	std::vector<std::array<PolyMesh::VertexHandle, 4>> quadVertices( pmesh.n_edges() );

	for( const auto& eh : pmesh.edges() ) {

		if( eh.is_boundary() )
			continue;

		const auto& heh1 = eh.h0();
		const auto& heh2 = eh.h1();

		if( heh1.face().valence() == 4 || heh2.face().valence() == 4 )
			continue;

		/*
		1-----0
		|   / |
		| /   |
		2-----3
		*/

		std::array<PolyMesh::VertexHandle, 4>vVec;
		vVec[0] = heh1.to();
		vVec[1] = heh1.next().to();
		vVec[2] = heh2.to();
		vVec[3] = heh2.next().to();

		std::array<PolyMesh::Point, 4>pVec;
		pVec[0] = pmesh.point( vVec[0] );
		pVec[1] = pmesh.point( vVec[1] );
		pVec[2] = pmesh.point( vVec[2] );
		pVec[3] = pmesh.point( vVec[3] );

		auto weight = QualityMetrics::conditionSurf( std::vector<PolyMesh::Point>(pVec.begin(), pVec.end()) );

		//prioQueue.insert( { weight, eh } );
		edges.push_back( { weight, eh } );
		quadVertices[eh.idx()] = vVec;
	}


	// delete edges (tracke welche faces zusammengefasst wurden)

	std::vector<bool>faceDeleted( pmesh.n_faces(), false );

	pmesh.request_face_status();
	pmesh.request_edge_status();
	pmesh.request_vertex_status();

	std::vector<PolyMesh::EdgeHandle> edgesToDelete;
	edgesToDelete.reserve( pmesh.n_edges() / 2 );

	std::sort( edges.begin(), edges.end(), []( const prioEdge& a, const prioEdge& b ) {return a.w > b.w; } );

	OpenMesh::EPropHandleT<int> eRank;
	pmesh.add_property( eRank );
	for( auto e : pmesh.edges() ) {
		pmesh.property( eRank, e ) = INT_MAX;
	}
	for( int i = 0; i < edges.size(); ++i ) {
		pmesh.property( eRank, edges[i].eh ) = i;
	}

	for( const auto& [w, eh, rank] : edges ) {
		// perform bfs
		std::queue<OpenMesh::SmartEdgeHandle> queue;
		queue.push( eh );
		while( !queue.empty() ) {
			auto e = queue.front();
			queue.pop();

			const auto& rank = pmesh.property( eRank, e );
			if( rank == INT_MAX )
				continue;

			const auto& fh1 = e.h0().face();
			const auto& fh2 = e.h1().face();

			// check if one of the faces is already part of a quad
			if( faceDeleted[fh1.idx()] || faceDeleted[fh2.idx()] )
				continue;
			if( edges[rank].w < 0.1 )
				continue;

			faceDeleted[fh1.idx()] = true;
			faceDeleted[fh2.idx()] = true;

			edgesToDelete.push_back( e );

			std::vector<OpenMesh::SmartEdgeHandle> nextEhs;

			nextEhs.push_back( e.h0().next().opp().prev().edge() );
			nextEhs.push_back( e.h0().next().opp().next().edge() );
			nextEhs.push_back( e.h0().prev().opp().prev().edge() );
			nextEhs.push_back( e.h0().prev().opp().next().edge() );

			nextEhs.push_back( e.h1().next().opp().prev().edge() );
			nextEhs.push_back( e.h1().next().opp().next().edge() );
			nextEhs.push_back( e.h1().prev().opp().prev().edge() );
			nextEhs.push_back( e.h1().prev().opp().next().edge() );

			std::sort( nextEhs.begin(), nextEhs.end(), [&pmesh, &eRank]( const OpenMesh::EdgeHandle& a, const OpenMesh::EdgeHandle& b ) {return pmesh.property( eRank, a ) < pmesh.property( eRank, b ); } );

			for( auto ne : nextEhs ) {
				queue.push( ne );
			}

			//for( auto ne : e.h0().to().edges() ) {
			//	queue.push( ne );
			//}
			//for( auto ne : e.h0().next().to().edges() ) {
			//	queue.push( ne );
			//}
			//for( auto ne : e.h1().to().edges() ) {
			//	queue.push( ne );
			//}
			//for( auto ne : e.h1().next().to().edges() ) {
			//	queue.push( ne );
			//}
		}
	}

	pmesh.remove_property( eRank );

	for( const auto& e : edgesToDelete ) {
		pmesh.delete_edge( e, false );
		const std::vector<PolyMesh::VertexHandle> vhs( quadVertices[e.idx()].begin(), quadVertices[e.idx()].end() );
		pmesh.add_face( vhs );
	}

	pmesh.garbage_collection();

}

/// <summary>
/// Perform quad generation with Blossom-Quad. The result is a hybrid mesh.
/// </summary>
/// <param name="tmesh">triangle mesh</param>
void RatRace::initBlossomQuad( TriMesh& tmesh ) {
	LOG_ASSERT( tmesh.n_vertices() > 0 );
	LOG_ASSERT( tmesh.n_faces() % 2 == 0 ) << " perfect matching can only be done for an even number of faces";
	LOG( INFO ) << "Use blossom quad";
	blossom_quad bloss;
	mesh_ = bloss.do_blossom_algo( (PolyMesh)tmesh );
}

/// <summary>
/// Perform quad generation with triangle merging. The result is a hybrid mesh.
/// </summary>
/// <param name="tmesh">triangle mesh</param>
void RatRace::initTriangleMerging( TriMesh& tmesh ) {
	LOG_ASSERT( tmesh.n_vertices() > 0 );
	LOG_ASSERT( tmesh.n_faces() % 2 == 0 ) << " perfect matching can only be done for an even number of faces";
	LOG( INFO ) << "Use triangle merging";
	mesh_ = tmesh;
	//tri2hybridSurf( mesh_ );
	tri2hybridBFS( mesh_ );
}

/// <summary>
/// Perform quad generation with a hybrid approach which combines Blossom-Quad and triangle merging. The result is a hybrid mesh.
/// </summary>
/// <param name="tmesh">triangle mesh</param>
void RatRace::initHybrid( TriMesh& tmesh ) {
	LOG_ASSERT( tmesh.n_vertices() > 0 );
	LOG_ASSERT( tmesh.n_faces() % 2 == 0 ) << " perfect matching can only be done for an even number of faces";
	LOG( INFO ) << "Use fast hybrid mesh generation";
	LOG( INFO ) << "Perform Blossom at boundaries";
	blossom_quad bl;
	// --- Perform blossom quad on the boundaries ---
	mesh_ = bl.quadrangleBoundaries( tmesh );
	//mesh_ = tmesh;
	LOG( INFO ) << "Greedy triangle merge on interior";
	//tri2hybridSurf( mesh_ );
	tri2hybridBFS( mesh_ );
}

/// <summary>
/// Merge remaining triangles in a greedy fashion.
/// </summary>
void RatRace::run() {
	DualHybridGraph dhg( mesh_ );
	LOG( INFO ) << "Number of triangles = " << dhg.nTriangles();

	dhg.connectTriangles();

	mesh_ = dhg.getMesh();
}

/// <summary>
/// Post processing on a quad mesh
/// </summary>
/// <param name="nQuads">Number of desired quads</param>
/// <param name="sizegrid">size function</param>
/// <param name="outputPath">output path for intermediate results (optional)</param>
void RatRace::postProcessing( const int& nQuads, SizeGrid* sizegrid, const std::experimental::filesystem::path& outputPath ) {
	while( true ) {
		int removed = removeInteriorValence2Vertices();
		if( removed == 0 ) break;
		LOG( INFO ) << "Removed interior val2: " << removed;
	}
	mesh_.garbage_collection();
	LOG( INFO ) << "Splitted boundary val2: " <<
		valence2SplitBoundary();
	LOG( INFO ) << "Removed boundary val2: " <<
		removeBoundaryValence2Vertices();
	LOG( INFO ) << "Removed boundary val2: " <<
		removeBoundaryValence2Vertices();
	while( true ) {
		int removed = removeInteriorValence2Vertices();
		if( removed == 0 ) break;
		LOG( INFO ) << "Removed interior val2: " << removed;
	}
	mesh_.garbage_collection();
	if(!outputPath.empty())
		OpenMesh::IO::write_mesh( mesh_, ( outputPath / "basicPostProcessingNoOpt.off" ).string(), OpenMesh::IO::Options::Default, 10 );

	discreteMeshOptimization<1>( mesh_, 0.5f, 5 );
	if( !outputPath.empty() )
		OpenMesh::IO::write_mesh( mesh_, ( outputPath / "basicPostProcessing.off" ).string(), OpenMesh::IO::Options::Default, 10 );

	while( true ) {
		int removed = lowQualSplitBoundary();
		if( removed == 0 ) break;
		LOG( INFO ) << "Removed low quality on boundary: " << removed;
	}
	//OpenMesh::IO::write_mesh( rr.mesh(), ( outputPath / "basicPostProcessing.off" ).string(), OpenMesh::IO::Options::Default, 10 );

	//SizeGrid sizegrid = loadSizeGrid( initMesh, cachePath, meshFile, 1000, 1000 );
	int nQuadsCurr = mesh_.n_faces();
	LOG( INFO ) << "Number of quads = " << nQuadsCurr;

	initRemeshProps();
	for( int k = 0; k < 50; ++k ) {
		int nSwaps = swapValence();
		int nColl = collapse3x3y( sizegrid );
		nQuadsCurr -= nColl;
		int nSplit = valenceSplit( nQuads - nQuadsCurr, sizegrid );
		nQuadsCurr += nSplit;
		LOG( INFO ) << k << ": swap / collapse / split: " << nSwaps << " / " << nColl << " / " << nSplit;
		if( nColl + nSplit == 0 )
			break;

		swapRemeshProps();
	}
	mesh_.garbage_collection();
	deleteRemeshProps();
	discreteMeshOptimization<1>( mesh_, 0.5f, 5 );
}

/// <summary>
/// Reduce irregular vertices by swapping edges
/// </summary>
/// <returns>number of edge swaps</returns>
int RatRace::swapValence() {
	int nSwaps = 0;

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	bool flippedSomething = false;
	do {
		flippedSomething = false;

		for( auto eh : mesh_.edges() ) {
			if( eh.is_boundary() ) continue;
			auto heh = eh.h0();
			
			//if( mesh_.property( vValence, heh.from() ) == 2 || mesh_.property( vValence, heh.to() ) == 2 ) {
			//	continue; 
			//}
			if( heh.from().valence() == 2 || heh.to().valence() == 2 ) {
				continue;
			}

			//	1 --- 0 --- 5
			//	|     ^     |
			//	|     |     |
			//	2 --- 3 --- 4

			auto hehInit = heh;
			std::array<OpenMesh::SmartVertexHandle, 6> vhs;
			vhs[0] = heh.to();
			heh = heh.next();
			vhs[1] = heh.to();
			heh = heh.next();
			vhs[2] = heh.to();
			heh = heh.next();
			vhs[3] = heh.to();
			heh = hehInit.opp().next();
			vhs[4] = heh.to();
			heh = heh.next();
			vhs[5] = heh.to();
			heh = heh.next();

			// check if edge is of interest
			if( !mesh_.property( vRemesh, vhs[0] ) && 
				!mesh_.property( vRemesh, vhs[1] ) &&
				!mesh_.property( vRemesh, vhs[2] ) &&
				!mesh_.property( vRemesh, vhs[3] ) &&
				!mesh_.property( vRemesh, vhs[4] ) &&
				!mesh_.property( vRemesh, vhs[5] ) ) {
				continue;
			}

			std::vector<int> vals;
			vals.reserve( 6 );
			for( const auto& vh : vhs ) {
				auto val = mesh_.property( vValence, vh );
				//auto val = boundaryValence(vh);
				//if( val != mesh_.property( vValence, vh ) ) {
				//	LOG( INFO ) << "Difference in valence found: " << val << " --- " << mesh_.property( vValence, vh );
				//}
				vals.push_back( val );
			}

			//if( vals[0] == 3 + vhs[0].is_boundary() || vals[3] == 3 + vhs[3].is_boundary() ) {
			//	continue;
			//}
			if( vals[0] <= 3 + vhs[0].is_boundary() || vals[3] <= 3 + vhs[3].is_boundary() ) {
				continue;
			}

			auto valOld = vals[0] + vals[3];
			auto valRight = vals[5] + vals[2];
			auto valLeft = vals[1] + vals[4];
			
			if( valOld - valRight >= valOld - valLeft && valOld - valRight >= 3 ) {
				auto l1 = ( mesh_.point( vhs[0] ) - mesh_.point( vhs[3] ) ).length();
				auto l2 = ( mesh_.point( vhs[2] ) - mesh_.point( vhs[5] ) ).length();
				if( l2 > 1.5 * l1 )
					continue;
				// swap to 2,5
				mesh_.edgeSwapPrev( eh );

				for( auto v : vhs ) {
					mesh_.property( vRemesh, v ) = true;
					mesh_.property( vRemeshNew, v ) = true;
					for( auto voh : v.outgoing_halfedges() ) {
						mesh_.property( vRemeshNew, voh.to() ) = true;
						mesh_.property( vRemeshNew, voh.next().to() ) = true;
						mesh_.property( vRemesh, voh.to() ) = true;
						mesh_.property( vRemesh, voh.next().to() ) = true;
					}
				}
				mesh_.property( vValence, vhs[0] ) -= 1;
				mesh_.property( vValence, vhs[3] ) -= 1;
				mesh_.property( vValence, vhs[2] ) += 1;
				mesh_.property( vValence, vhs[5] ) += 1;

				flippedSomething = true;
				++nSwaps;
			} else if( valOld - valLeft > valOld - valRight && valOld - valLeft >= 3 ) {
				auto l1 = ( mesh_.point( vhs[0] ) - mesh_.point( vhs[3] ) ).length();
				auto l2 = ( mesh_.point( vhs[1] ) - mesh_.point( vhs[4] ) ).length();
				if( l2 > 1.5 * l1 )
					continue;
				// swap to 1,4
				mesh_.edgeSwapNext( eh );

				for( auto v : vhs ) {
					mesh_.property( vRemesh, v ) = true;
					mesh_.property( vRemeshNew, v ) = true;
					for( auto voh : v.outgoing_halfedges() ) {
						mesh_.property( vRemeshNew, voh.to() ) = true;
						mesh_.property( vRemeshNew, voh.next().to() ) = true;
						mesh_.property( vRemesh, voh.to() ) = true;
						mesh_.property( vRemesh, voh.next().to() ) = true;
					}
				}
				mesh_.property( vValence, vhs[0] ) -= 1;
				mesh_.property( vValence, vhs[3] ) -= 1;
				mesh_.property( vValence, vhs[1] ) += 1;
				mesh_.property( vValence, vhs[4] ) += 1;
				
				flippedSomething = true;
				++nSwaps;
			}

		}
	} while( flippedSomething );

	//mesh_.garbage_collection();

	return nSwaps;
}

/// <summary>
/// Remove interior vertices with only two incident edges
/// </summary>
/// <returns>number of removed vertices</returns>
int RatRace::removeInteriorValence2Vertices() {
	int nCollapses = 0;

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();
	for( auto vh : mesh_.vertices() ) {
		if( vh.is_boundary() )
			continue;
		if( vh.valence() == 2 ) {
			auto heh = vh.halfedge();
			if( heh.next().to().valence() == 2 )
				mesh_.diagonalCollapse( heh );
			else
				mesh_.diagonalCollapseDirect( heh );
			++nCollapses;
		}
	}
	//mesh_.garbage_collection();

	return nCollapses;
}

/// <summary>
/// Vertex split on boundary vertices which have only two incident edges and low quality.
/// </summary>
/// <returns>number of splits</returns>
int RatRace::valence2SplitBoundary() {
	int nSplits = 0;

	auto isLegalSplit = [this]( OpenMesh::SmartHalfedgeHandle& heh ) {
		if( heh.edge().is_boundary() )
			return false;

		auto vh = heh.from();

		if( vh.valence() >= 4 )
			return true;

		if( vh.is_boundary() ) {
			real angle = calcBoundaryAngle( mesh_, vh );
			if( angle < 130 ) {
				return true;
			}
		}

		return false;
	};

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	for( auto vhBegin : mesh_.vertices() ) {
		if( !vhBegin.is_boundary() || vhBegin.valence() != 2 )
			continue;
	
		// check for vertex quality
		auto heh = vhBegin.halfedge();
		if( heh.is_boundary() )
			heh = heh.opp().next();
		auto p = mesh_.point( vhBegin );
		auto pNext = mesh_.point( heh.to() );
		auto pPrev = mesh_.point( heh.prev().from() );
		auto qv = QualityMetrics::conditionSurfVertex( { pPrev, p, pNext } );
	
		if( qv >= 0.1 )
			continue;
	
		// find all possible splits
		std::vector<OpenMesh::SmartHalfedgeHandle> hehCand;
		for( auto voh : vhBegin.outgoing_halfedges() ) {
			if( !voh.is_boundary() ) {
				auto cand = voh.next();
				if( isLegalSplit( cand ) )
					hehCand.push_back( cand );
			}
			auto vih = voh.opp();
			if( !vih.is_boundary() ) {
				auto cand = vih.prev().opp();
				if( isLegalSplit( cand ) )
					hehCand.push_back( cand );
			}
		}
	
		if( hehCand.empty() )
			continue;
	
		OpenMesh::SmartHalfedgeHandle hehSplit;
		auto minValence = INT_MAX;
		for( const auto& heh : hehCand ) {
			auto vh1 = heh.opp().next().to();
			auto vh2 = heh.prev().from();
			auto val = vh1.valence() + vh2.valence();
			if( val < minValence ) {
				minValence = val;
				hehSplit = heh;
			}
		}
		mesh_.vertexSplit( hehSplit );
		++nSplits;
	}
	mesh_.garbage_collection();
	return nSplits;
}

/// <summary>
/// Vertex split on boundary vertices with low quality, independent of their valence.
/// </summary>
/// <returns>number of splits</returns>
int RatRace::lowQualSplitBoundary() {
	int nSplits = 0;

	auto isLegalSplit = [this]( OpenMesh::SmartHalfedgeHandle& heh ) {
		if( heh.edge().is_boundary() )
			return false;

		auto vh = heh.from();

		if( vh.valence() >= 4 )
			return true;

		if( vh.is_boundary() ) {
			real angle = calcBoundaryAngle( mesh_, vh );
			if( angle < 130 ) {
				return true;
			}
		}

		return false;
	};

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	for( auto fh : mesh_.faces() ) {
		if( !fh.is_boundary() )
			continue;

		std::vector<PolyMesh::Point> points;
		points.reserve( 4 );
		for( auto fv : mesh_.fv_range( fh ) ) {
			points.push_back( mesh_.point( fv ) );
		}
		auto q = QualityMetrics::conditionSurf( points );

		if( q >= 0.1 )
			continue;

		// find low quality vertex
		auto heh = fh.halfedge();
		auto hehIt = heh;

		while( true ) {
			auto p = mesh_.point( hehIt.from() );
			auto pNext = mesh_.point( hehIt.to() );
			auto pPrev = mesh_.point( hehIt.prev().from() );
			auto qv = QualityMetrics::conditionSurfVertex( { pPrev, p, pNext } );

			if( qv < 0.1 )
				break;
			hehIt = hehIt.next();
		}
		heh = hehIt;

		// find all possible splits
		std::vector<OpenMesh::SmartHalfedgeHandle> hehCand;
		if( isLegalSplit( heh.next() ) )
			hehCand.push_back( heh.next() );
		if( isLegalSplit( heh.prev().prev().opp() ) )
			hehCand.push_back( heh.prev().prev().opp() );

		if( hehCand.empty() )
			continue;

		OpenMesh::SmartHalfedgeHandle hehSplit;
		auto minValence = INT_MAX;
		for( const auto& heh : hehCand ) {
			auto vh1 = heh.opp().next().to();
			auto vh2 = heh.prev().from();
			auto val = vh1.valence() + vh2.valence();
			if( val < minValence ) {
				minValence = val;
				hehSplit = heh;
			}
		}
		mesh_.vertexSplit( hehSplit );
		++nSplits;
	}

	mesh_.garbage_collection();
	return nSplits;
}

/// <summary>
/// Diagonal collapse on low quality valence two vertices on the boundary
/// </summary>
/// <returns>number of collapses</returns>
int RatRace::removeBoundaryValence2Vertices() {
	int nCollapses = 0;

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_halfedge_status();
	mesh_.request_vertex_status();

	for( auto vh : mesh_.vertices() ) {
		if( !vh.is_boundary() )
			continue;
		if( vh.valence() != 2 )
			continue;

		auto heh = vh.halfedge();
		if( heh.is_boundary() )
			heh = heh.opp().next();

		auto p = mesh_.point( vh );
		auto pNext = mesh_.point( heh.to() );
		auto pPrev = mesh_.point( heh.prev().from() );
		auto qv = QualityMetrics::conditionSurfVertex( { pPrev, p, pNext } );

		if( qv < 0.1 ) {
			// remove this vertex
			auto hehColl = heh.prev();
			if( hehColl.from().valence() == 2 || hehColl.prev().from().valence() == 2 || hehColl.next().to().valence() == 2 ) {
				if( hehColl.from().valence() == 2 ) {
					auto heh2 = hehColl.opp();
					auto heh1 = heh2.prev();
					mesh_.collapse( heh1 );
					mesh_.collapse( heh2 );
				} else if( hehColl.next().to().valence() == 2 ) {
					auto heh1 = hehColl;
					auto heh2 = hehColl.next();
					mesh_.collapse( heh1 );
					mesh_.collapse( heh2 );
				} else {
					LOG( WARNING ) << "Quad with two valence 2 vertices on boundary";
				}
				continue;
			}
			auto vhRemain = hehColl.next().to();
			// decide on position
			auto a = abs( calcBoundaryAngle( mesh_, vh ) - 180 );
			auto aNext = abs( calcBoundaryAngle( mesh_, heh.to() ) - 180 );
			auto aPrev = abs( calcBoundaryAngle( mesh_, heh.prev().from() ) - 180 );
			if( aNext < 10 && aPrev < 10 )
				mesh_.set_point( vhRemain, p );
			else if( aNext > aPrev )
				mesh_.set_point( vhRemain, pNext );
			else if( aPrev > aNext )
				mesh_.set_point( vhRemain, pPrev );
			mesh_.diagonalCollapseDirect( hehColl );
			++nCollapses;
		}

	}
	mesh_.garbage_collection();
	return nCollapses;
}

/// <summary>
/// Diagonal collapse on the interior constraint by a size function
/// </summary>
/// <param name="sizegrid">size function</param>
/// <returns>number of collapses</returns>
int RatRace::collapse3x3y( SizeGrid* sizegrid ) {
	int nCollapses = 0;

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	for( auto fh : mesh_.faces() ) {
		if( fh.is_boundary() )
			continue;

		std::array<OpenMesh::SmartVertexHandle, 4> vhs;
		std::array<int, 4> val;
		
		auto h = fh.halfedge();
		vhs[0] = h.from();
		vhs[1] = h.to();
		vhs[2] = h.next().to();
		vhs[3] = h.prev().from();
		//val[0] = boundaryValence( vhs[0] );
		//val[1] = boundaryValence( vhs[1] );
		//val[2] = boundaryValence( vhs[2] );
		//val[3] = boundaryValence( vhs[3] );
		val[0] = mesh_.property(vValence, vhs[0]);
		val[1] = mesh_.property(vValence, vhs[1]);
		val[2] = mesh_.property(vValence, vhs[2]);
		val[3] = mesh_.property(vValence, vhs[3]);

		// check if edge is of interest
		//bool isMarked = false;
		if( !mesh_.property( vRemesh, vhs[0] ) &&
			!mesh_.property( vRemesh, vhs[1] ) &&
			!mesh_.property( vRemesh, vhs[2] ) &&
			!mesh_.property( vRemesh, vhs[3] ) ) {
			continue;
			//isMarked = true;
		}

		int valCurr = ( val[0] == 4 ) + ( val[1] == 4 ) + ( val[2] == 4 ) + ( val[3] == 4 );

		if( ( val[0] + val[2] - 2 == 4 ) + ( val[1] - 1 == 4 ) + ( val[3] - 1 == 4 ) > valCurr ) {
			if( vhs[0].is_boundary() || vhs[2].is_boundary() )
				continue;
			if( val[1] == 3 || val[3] == 3 )
				continue;
			auto v0 = vhs[0];
			auto v2 = vhs[2];
			auto p0 = mesh_.point( v0 );
			auto p2 = mesh_.point( v2 );
			if( sizegrid != nullptr ) {
				auto l = ( ( p0 - p2 ).length() / sizegrid->getSize( p0, p2 ) );
				if( l > 1.5f )
					continue;
			}
			//LOG_ASSERT( !isMarked );
			for( auto v : vhs ) {
				mesh_.property( vRemesh, v ) = true;
				mesh_.property( vRemeshNew, v ) = true;
				for( auto voh : v.outgoing_halfedges() ) {
					mesh_.property( vRemeshNew, voh.to() ) = true;
					mesh_.property( vRemeshNew, voh.next().to() ) = true;
					mesh_.property( vRemesh, voh.to() ) = true;
					mesh_.property( vRemesh, voh.next().to() ) = true;
				}
			}
			mesh_.property( vValence, vhs[0] ) = val[0] + val[2] - 2;
			mesh_.property( vValence, vhs[2] ) = val[0] + val[2] - 2;
			mesh_.property( vValence, vhs[1] ) -= 1;
			mesh_.property( vValence, vhs[3] ) -= 1;

			mesh_.set_point( v2, 0.5f * ( p0 + p2 ) );
			mesh_.diagonalCollapseDirect( h );
			++nCollapses;
		} else if( ( val[1] + val[3] - 2 == 4 ) + ( val[0] - 1 == 4 ) + ( val[2] - 1 == 4 ) > valCurr ) {
			if( vhs[1].is_boundary() || vhs[3].is_boundary() )
				continue;
			if( val[0] == 3 || val[2] == 3 )
				continue;
			auto v1 = vhs[1];
			auto v3 = vhs[3];
			auto p1 = mesh_.point( v1 );
			auto p3 = mesh_.point( v3 );
			if( sizegrid != nullptr ) {
				auto l = ( ( p1 - p3 ).length() / sizegrid->getSize( p1, p3 ) );
				if( l > 1.5f )
					continue;
			}
			//LOG_ASSERT( !isMarked );
			for( auto v : vhs ) {
				mesh_.property( vRemesh, v ) = true;
				mesh_.property( vRemeshNew, v ) = true;
				for( auto voh : v.outgoing_halfedges() ) {
					mesh_.property( vRemeshNew, voh.to() ) = true;
					mesh_.property( vRemeshNew, voh.next().to() ) = true;
					mesh_.property( vRemesh, voh.to() ) = true;
					mesh_.property( vRemesh, voh.next().to() ) = true;
				}
			}
			mesh_.property( vValence, vhs[1] ) = val[1] + val[3] - 2;
			mesh_.property( vValence, vhs[3] ) = val[1] + val[3] - 2;
			mesh_.property( vValence, vhs[0] ) -= 1;
			mesh_.property( vValence, vhs[2] ) -= 1;

			mesh_.set_point( v3, 0.5f * ( p1 + p3 ) );
			mesh_.diagonalCollapseDirect( h.next() );
			++nCollapses;
		}

	}

	//mesh_.garbage_collection();

	return nCollapses;
}

/// <summary>
/// Vertex split on the interior constraint by a size function
/// </summary>
/// <param name="nMaxSplits">maximal number of splits</param>
/// <param name="sizegrid">size function</param>
/// <returns>number of splits</returns>
int RatRace::valenceSplit( int nMaxSplits, SizeGrid* sizegrid ) {
	int nSplits = 0;

	if( nMaxSplits == 0 ) {
		return 0;
	}

	auto isLegalSplit = [this]( OpenMesh::SmartHalfedgeHandle& heh ) {
		if( heh.edge().is_boundary() )
			return false;

		auto vh = heh.from();

		//if( boundaryValence( vh ) > 4 )
		//	return true;
		if( mesh_.property( vValence, vh ) > 4 )
			return true;

		return false;
	};

	mesh_.request_face_status();
	mesh_.request_edge_status();
	mesh_.request_vertex_status();

	for( auto vh : mesh_.vertices() ) {

		if( vh.is_boundary() )
			continue;

		//auto val = boundaryValence( vh );
		auto val = mesh_.property( vValence, vh );


		if( val > 3 )
			continue;

		// find all possible splits
		std::vector<OpenMesh::SmartHalfedgeHandle> hehCand;
		for( auto voh : vh.outgoing_halfedges() ) {
			if( !voh.is_boundary() ) {
				auto cand = voh.next();
				if( isLegalSplit( cand ) )
					hehCand.push_back( cand );
			}
			auto vih = voh.opp();
			if( !vih.is_boundary() ) {
				auto cand = vih.prev().opp();
				if( isLegalSplit( cand ) )
					hehCand.push_back( cand );
			}
		}

		if( hehCand.empty() )
			continue;

		OpenMesh::SmartHalfedgeHandle hehSplit;
		auto minValence = INT_MAX;
		for( const auto& heh : hehCand ) {
			auto vh1 = heh.opp().next().to();
			auto vh2 = heh.prev().from();
			//auto val = boundaryValence( vh1 ) + boundaryValence( vh2 );
			auto val = mesh_.property( vValence, vh1 ) + mesh_.property( vValence, vh2 );
			if( val < minValence ) {
				minValence = val;
				hehSplit = heh;
			}
		}

		std::array<OpenMesh::SmartVertexHandle, 4> vhs;
		vhs[0] = hehSplit.from();
		vhs[1] = hehSplit.opp().next().to();
		vhs[2] = hehSplit.prev().from();
		// check if edge is of interest
		if( !mesh_.property( vRemesh, vhs[0] ) &&
			!mesh_.property( vRemesh, vhs[1] ) &&
			!mesh_.property( vRemesh, vhs[2] ) ) {
			continue;
		}

		if( sizegrid != nullptr ) {
			//auto l = mesh_.calc_edge_length( hehSplit );
			auto l = ( mesh_.point( hehSplit.from() ) - mesh_.point( hehSplit.to() ) ).length();
			//auto h = sizegrid->getSize( mesh_.point( vh ) );
			auto h = sizegrid->getSize( mesh_.point( hehSplit.from() ), mesh_.point( hehSplit.to() ) );
			auto r = ( l / h );
			if( r < 0.9f )
				continue;
		}

		auto vhNew = mesh_.vertexSplit( hehSplit );
		vhs[3] = vhNew;
		for( auto v : vhs ) {
			mesh_.property( vRemesh, v ) = true;
			mesh_.property( vRemeshNew, v ) = true;
			for( auto voh : v.outgoing_halfedges() ) {
				mesh_.property( vRemeshNew, voh.to() ) = true;
				mesh_.property( vRemeshNew, voh.next().to() ) = true;
				mesh_.property( vRemesh, voh.to() ) = true;
				mesh_.property( vRemesh, voh.next().to() ) = true;
			}
		}
		mesh_.property( vValence, vhs[0] ) -= 1;
		mesh_.property( vValence, vhs[1] ) += 1;
		mesh_.property( vValence, vhs[2] ) += 1;
		mesh_.property( vValence, vhs[3] ) = 3;

		++nSplits;
		if( nSplits == nMaxSplits ) {
			break;
		}
	}
	//mesh_.garbage_collection();

	return nSplits;
}

/// <summary>
/// Initialize properties for remeshing. This method must be called before remeshing can be performed!
/// </summary>
void RatRace::initRemeshProps() {
	mesh_.add_property( vRemesh );
	mesh_.add_property( vRemeshNew );
	mesh_.add_property( vValence );

	for( auto vh : mesh_.vertices() ) {
		if(vh.valence() != 4)
			mesh_.property( vRemesh, vh ) = true;
		else
			mesh_.property( vRemesh, vh ) = false;
		mesh_.property( vRemeshNew, vh ) = false;

		mesh_.property( vValence, vh ) = boundaryValence( vh );
	}
}

/// <summary>
/// Swap remeshing properties. This method must be called after each remeshing iteration
/// </summary>
void RatRace::swapRemeshProps() {
	for( auto vh : mesh_.vertices() ) {
		mesh_.property( vRemesh, vh ) = mesh_.property( vRemeshNew, vh );
		//mesh_.property( vRemesh, vh ) = true;	// for debug
		mesh_.property( vRemeshNew, vh ) = false;

		//mesh_.property( vValence, vh ) = boundaryValence( vh );
	}
}

/// <summary>
/// Delete remesh properties. This method should be called after remeshing.
/// </summary>
void RatRace::deleteRemeshProps() {
	mesh_.remove_property( vRemesh );
	mesh_.remove_property( vRemeshNew );
	mesh_.remove_property( vValence );
}

/// <summary>
/// Compute valence considering the boundary shape if the vertex is on the boundary.
/// </summary>
/// <param name="vh">vertex handle</param>
/// <returns>valence</returns>
int RatRace::boundaryValence( const OpenMesh::SmartVertexHandle& vh ) {
	auto val = vh.valence();
	if( vh.is_boundary() ) {
		auto angle = 360.f - calcBoundaryAngle( mesh_, vh );
		if( angle < 30 )
			val -= 1;
		else if( angle < 120 )
			;	// do not change valence
		else if( angle < 210 )
			val += 1;
		else
			val += 2;
	}
	return val;
}
