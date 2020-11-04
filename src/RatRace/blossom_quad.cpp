#include "blossom_quad.h"

#include "QualityMetrics.h"
#include <queue>
#include "workingDirectory.h"

using namespace std;

PolyMesh blossom_quad::do_blossom_algo(PolyMesh loadedMesh ) {

	// raw triangle mesh
	workMesh = loadedMesh;
	// property for edge in connection graph
	OpenMesh::HPropHandleT<int>	hprop_graph_edge_idx;
	workMesh.add_property( hprop_graph_edge_idx, "hprop_graph_edge_idx" );
	for( auto h : workMesh.halfedges() ) {
		workMesh.property( hprop_graph_edge_idx, h ) = -1;
	}
	// property for external edge in connection graph
	OpenMesh::VPropHandleT<int> vprop_external_edge_idx;
	workMesh.add_property( vprop_external_edge_idx, "vprop_external_edge_idx" );
	for( auto v : workMesh.vertices() ) {
		workMesh.property( vprop_external_edge_idx, v ) = -1;
	}

	// connectivity graph 
	con_graph graph;

	int node_number = workMesh.n_faces();
	LOG_ASSERT( node_number % 2 == 0 );
	
	/*
	 * go through all halfedges of the mesh and create a connectivity graph containing all faces
	 * this graph is needed for the blossom algorithm
	 */
	for( auto heh : workMesh.halfedges() ) {
		auto fh = heh.face();
		auto fhOpp = heh.opp().face();

		/*
		 * check, if edge is border edge
		 * check, if edge is already in connectivity graph
		 */
		if( !heh.edge().is_boundary() && workMesh.property( hprop_graph_edge_idx, heh ) == -1 ) {
			// create new graph edge in connectivity graph
			graph_edge gEdge;
			gEdge.face0 = fh;
			gEdge.face1 = fhOpp;
			// calculate vertices for possible quad
			auto quadVhs = connect_triangles( heh );
			gEdge.quad_vertices = quadVhs;

			// calculate cost for graph edge
			std::vector<PolyMesh::Point> points{ workMesh.point( quadVhs[0] ), workMesh.point( quadVhs[1] ), workMesh.point( quadVhs[2] ), workMesh.point( quadVhs[3] ) };
			real qual = 1.f - QualityMetrics::condition( points );
			gEdge.cost = (int)( qual * 100 );

			// add edge to connection graph
			graph.edges.push_back( gEdge );
			// mark mesh edge as connected
			workMesh.property( hprop_graph_edge_idx, heh ) = graph.edges.size() - 1;
			workMesh.property( hprop_graph_edge_idx, heh.opp() ) = graph.edges.size() - 1;
		}
		/*
		 * check, if halfedge is boundary
		 * if so, add an external edge to the connectivity graph
		 */
		else if( workMesh.is_boundary( heh ) ) {
			auto vhFrom = heh.from();
			auto vhTo = heh.to();

			// add external edge to the from vertex
			if( workMesh.property( vprop_external_edge_idx, vhFrom ) == -1 ) {
				//PolyMesh::FaceHandle externalConnectionface;
				int valence = vhFrom.valence();

				//for valence < 4, do not add external edge
				if( valence < 4 )
					continue;

				//border face, connected with this vertex 
				//auto fhPrevOpp = workMesh.opposite_face_handle( workMesh.prev_halfedge_handle( heh ) );
				auto fhPrevOpp = heh.prev().opp().face();
				// create new external edge for connectivity graph
				graph_edge external_edge;
				external_edge.face0 = fhOpp;
				external_edge.face1 = fhPrevOpp;
				external_edge.connecting_vertex = vhFrom;
				external_edge.cost = EXTERNAL_EDGE_COST;

				workMesh.property( vprop_external_edge_idx, vhFrom ) = graph.external_edges.size();
				graph.external_edges.push_back( external_edge );
			}
		}
	}

	/* 
	 * create data structures needed for the blossom algorithm
	 * in the edges array there are the indices of the start and end index of every edge in the connectivity graph, so every index is a face index
	 * in the weights array there are the calculated weights of every edge in the connectivity graph
	 */
	std::vector<int> edges, weights;
	int edge_number = graph.edges.size() + graph.external_edges.size();
	edges.reserve( 2 * edge_number );
	weights.reserve( edge_number );
	for( const auto& e : graph.edges ) {
		weights.push_back( e.cost );
		edges.push_back( e.face0.idx() );
		edges.push_back( e.face1.idx() );
	}
	for( const auto& e : graph.external_edges ) {
		weights.push_back( e.cost );
		edges.push_back( e.face0.idx() );
		edges.push_back( e.face1.idx() );
	}

	/*
	 * beginning of creation of the result mesh
	 * first all points are copied to the new mesh
	 */
	int point_count = 0;
	for( auto vh : workMesh.vertices() ) {
		blossomMesh_step1_vertex_vec.push_back( blossomMesh_step1.add_vertex( workMesh.point( vh ) ) );
		point_count++;
	}

	// use the blossom V algorithm to find a cost perfect minimal matching
	PerfectMatching pm( node_number, edge_number );
	for( int e = 0; e < edge_number; e++ )
		pm.AddEdge( edges[2 * e], edges[2 * e + 1], weights[e] );
	pm.options.update_duals_before = true;
	pm.options.fractional_jumpstart = false;
	pm.Solve();

	constexpr bool check_perfect_matching = true;
	if constexpr (check_perfect_matching)
	{
		int res = CheckPerfectMatchingOptimality(node_number, edge_number, edges.data(), weights.data(), &pm);
		std::printf("check optimality: res=%d (%s)\n", res, (res == 0) ? "ok" : ((res == 1) ? "error" : "fatal error"));
	}
	double cost = ComputePerfectMatchingCost(node_number, edge_number, edges.data(), weights.data(), &pm);
	std::printf("cost = %.1f\n", cost);
	
	/*
	 * after a matching is found, use it to create a quad mesh from the connectivity graph
	 * if an edge is used, add its quad face to the match
	 * if an external edge is used, mark it for later usage
	 */
	for ( int firstFaceIdx = 0; firstFaceIdx < node_number; firstFaceIdx++) {
		int secondFaceIdx = pm.GetMatch( firstFaceIdx );
		
		if (firstFaceIdx < secondFaceIdx) {
			PolyMesh::FaceHandle firstFace = workMesh.face_handle(firstFaceIdx);
			for( auto firstFaceHalfedge : workMesh.fh_range( firstFace ) ) {
				int edgeId = workMesh.property( hprop_graph_edge_idx, firstFaceHalfedge );
				if( edgeId < 0 )
					continue;
				graph_edge edgeToUse = graph.edges[edgeId];
				// add face to mesh
				if( ( edgeToUse.face0.idx() == firstFaceIdx && edgeToUse.face1.idx() == secondFaceIdx ) || ( edgeToUse.face0.idx() == secondFaceIdx && edgeToUse.face1.idx() == firstFaceIdx ) ) {
					blossomMesh_step1.add_face( edgeToUse.quad_vertices );
					break;
				}
				// mark external edge for usage
				PolyMesh::VertexHandle fromVertex = workMesh.from_vertex_handle( firstFaceHalfedge );
				if( workMesh.property( vprop_external_edge_idx, fromVertex ) != -1 ) {
					graph_edge externalEdgeToMark = graph.external_edges[workMesh.property( vprop_external_edge_idx, fromVertex )];
					if( ( externalEdgeToMark.face0.idx() == firstFaceIdx && externalEdgeToMark.face1.idx() == secondFaceIdx ) || ( externalEdgeToMark.face0.idx() == secondFaceIdx && externalEdgeToMark.face1.idx() == firstFaceIdx ) ) {
						graph.external_edges[workMesh.property( vprop_external_edge_idx, fromVertex )].to_be_used = true;
						break;
					}
				}
				PolyMesh::VertexHandle toVertex = workMesh.to_vertex_handle( firstFaceHalfedge );
				if( workMesh.property( vprop_external_edge_idx, toVertex ) != -1 ) {
					graph_edge externalEdgeToMark = graph.external_edges[workMesh.property( vprop_external_edge_idx, toVertex )];
					if( ( externalEdgeToMark.face0.idx() == firstFaceIdx && externalEdgeToMark.face1.idx() == secondFaceIdx ) || ( externalEdgeToMark.face0.idx() == secondFaceIdx && externalEdgeToMark.face1.idx() == firstFaceIdx ) ) {
						graph.external_edges[workMesh.property( vprop_external_edge_idx, toVertex )].to_be_used = true;
						break;
					}
				}
			}
		}		
	}
	
	// add faces which are connected by external edges
	for( int i = 0; i < graph.external_edges.size(); i++ ) {
		// find edges, that are marked for use
		if( graph.external_edges[i].to_be_used ) {
			graph.external_edges[i].to_be_used = false;

			// get the two faces to connect and the connecting vertex from the external edge
			PolyMesh::FaceHandle face0 = graph.external_edges[i].face0;
			PolyMesh::FaceHandle face1 = graph.external_edges[i].face1;

			/*
			 * create and add new quad face
			 * use vertices from the first triangle face and the duplicated vertex
			 */
			vector<PolyMesh::VertexHandle> quad_face0;
			for( auto vh_face0 : workMesh.fv_range( face0 ) ) {
				quad_face0.push_back( blossomMesh_step1_vertex_vec[vh_face0.idx()] );
			}
			blossomMesh_step1.add_face( quad_face0 );
			/*
			* create and add new quad face
			* use vertices from the second triangle face and the duplicated vertex
			*/
			vector<PolyMesh::VertexHandle> quad_face1;
			for( auto vh_face1 : workMesh.fv_range( face1 ) ) {
				quad_face1.push_back( blossomMesh_step1_vertex_vec[vh_face1.idx()] );
			}
			blossomMesh_step1.add_face( quad_face1 );
		}
	}

	// return the created hybrid mesh
	return blossomMesh_step1;
}

void testPrint( TriMesh& mesh, const OpenMesh::FPropHandleT<int>& fBoundaryId, const OpenMesh::VPropHandleT<int>& vBloss ) {
	mesh.request_face_colors();
	mesh.request_vertex_colors();
	for( auto fh : mesh.faces() ) {
		switch( mesh.property( fBoundaryId, fh ) ) {
		case -1:
			mesh.set_color( fh, { 255,255,255 } );
			break;
		case 0:
			mesh.set_color( fh, { 255,0,0 } );
			break;
		case 1:
			mesh.set_color( fh, { 0,255,0 } );
			break;
		case 2:
			mesh.set_color( fh, { 0,0,255 } );
			break;
		case 3:
			mesh.set_color( fh, { 255,255,0 } );
			break;
		case 4:
			mesh.set_color( fh, { 255,0,255 } );
			break;
		case 5:
			mesh.set_color( fh, { 0,255,255 } );
			break;
		case 6:
			mesh.set_color( fh, { 128,255,0 } );
			break;
		case 7:
			mesh.set_color( fh, { 0,128,255 } );
			break;
		case 8:
			mesh.set_color( fh, { 255,0,128 } );
			break;
		default:
			mesh.set_color( fh, { mesh.property( fBoundaryId, fh ),mesh.property( fBoundaryId, fh ),mesh.property( fBoundaryId, fh ) } );
			break;
		}
	}

	for( auto vh : mesh.vertices() ) {
		switch( mesh.property( vBloss, vh ) ) {
		case -1:
			mesh.set_color( vh, { 255,255,255 } );
			break;
		case 0:
			mesh.set_color( vh, { 255,0,0 } );
			break;
		case 1:
			mesh.set_color( vh, { 0,255,0 } );
			break;
		case 2:
			mesh.set_color( vh, { 0,0,255 } );
			break;
		case 3:
			mesh.set_color( vh, { 255,255,0 } );
			break;
		case 4:
			mesh.set_color( vh, { 255,0,255 } );
			break;
		case 5:
			mesh.set_color( vh, { 0,255,255 } );
			break;
		case 6:
			mesh.set_color( vh, { 128,255,0 } );
			break;
		case 7:
			mesh.set_color( vh, { 0,128,255 } );
			break;
		case 8:
			mesh.set_color( vh, { 255,0,128 } );
			break;
		default:
			mesh.set_color( vh, { mesh.property( vBloss, vh ),mesh.property( vBloss, vh ),mesh.property( vBloss, vh ) } );
			break;
		}
	}
	OpenMesh::IO::write_mesh( mesh, std::string(WORKING_DIRECTORY) + "output/ring.off", OpenMesh::IO::Options::FaceColor | OpenMesh::IO::Options::VertexColor );
}

inline void fillDomainCorners( TriMesh& mesh, const OpenMesh::VPropHandleT<int>& vBoundaryId, const int& nBoundaries ) {
	// mark vertices in domain corners
	bool addSomething = false;
	do {
		addSomething = false;
		for( auto vh : mesh.vertices() ) {
			if( mesh.property( vBoundaryId, vh ) != -1 )
				continue;
			std::map<int, int> idOcc;
			for( auto vv : vh.vertices() ) {
				const auto& id = mesh.property( vBoundaryId, vv );
				if( id != -1 ) {
					auto it = idOcc.find( id );
					if( it != idOcc.end() ) {
						++(it->second);
					} else {
						idOcc[id] = 1;
					}
				}
			}
			
			//int maxRingNeigh = 3;	// at least 3 vertices must have same boundary id to be considered
			int maxRingNeigh = vh.valence() / 2;
			for( const auto& [id, n] : idOcc ) {
				if( n > maxRingNeigh ) {
					maxRingNeigh = id;
					mesh.property( vBoundaryId, vh ) = id;
					addSomething = true;
				}
			}

			// TODO this version might be faster but must be checked thoroughly for correctness
			//int nNeighs = 0;
			//int bId = -1;
			//for( auto vv : vh.vertices() ) {
			//	if( mesh.property( vBoundaryId, vv ) != -1 ) {
			//		++nNeighs;
			//		bId = mesh.property( vBoundaryId, vv );
			//	}
			//}
			//if( nNeighs > 3 ) {
			//	mesh.property( vBoundaryId, vh ) = bId;
			//	addSomething = true;
			//}
		}

	} while( addSomething );

	//bool addSomething = false;
	//do {
	//	addSomething = false;
	//	for( auto vh : mesh.vertices() ) {
	//		if( mesh.property( vBoundaryId, vh ) != -1 )
	//			continue;
	//		int neighCol = -1;
	//		int nNeighs = 0;
	//		for( auto vv : vh.vertices() ) {
	//			if( mesh.property( vBoundaryId, vv ) != -1 ) {
	//				++nNeighs;
	//				neighCol = mesh.property( vBoundaryId, vv );
	//			}
	//		}
	//
	//		if( nNeighs > 3 ) {
	//			mesh.property( vBoundaryId, vh ) = neighCol;
	//			addSomething = true;
	//		}
	//	}
	//
	//} while( addSomething );
}

PolyMesh blossom_quad::quadrangleBoundaries( TriMesh& mesh ) {

	struct BlossomEdge
	{
		OpenMesh::SmartFaceHandle f0, f1;
		OpenMesh::SmartEdgeHandle e;
		int w;
		int id = -1;
	};

	OpenMesh::VPropHandleT<int> vBoundaryId;
	OpenMesh::FPropHandleT<int> fBoundaryId;
	mesh.add_property( vBoundaryId );
	mesh.add_property( fBoundaryId );

	int nBoundaries = 0;
	for( auto v : mesh.vertices() ) {
		mesh.property( vBoundaryId, v ) = -1;
	}
	for( auto e : mesh.edges() ) {
		if( !e.is_boundary() )
			continue;

		auto h = e.h0();
		if( !h.is_boundary() )
			h = e.h1();

		if( mesh.property( vBoundaryId, h.from() ) != -1 )
			continue;

		auto hIt = h;
		do {
			mesh.property( vBoundaryId, hIt.from() ) = nBoundaries;
			hIt = hIt.next();
		} while( hIt != h );
		++nBoundaries;
	}
	for( auto f : mesh.faces() ) {
		mesh.property( fBoundaryId, f ) = -1;
	}

	OpenMesh::FPropHandleT<int> fBloss;
	mesh.add_property( fBloss );

	// mark vertices near boundaries
	for( auto vh : mesh.vertices() ) {
		if( !vh.is_boundary() )
			continue;
		auto id = mesh.property( vBoundaryId, vh );
		for( auto vv : vh.vertices() ) {
			mesh.property( vBoundaryId, vv ) = id;
			for( auto vvv : vv.vertices() ) {
				mesh.property( vBoundaryId, vvv ) = id;
				//for( auto vvvv : vvv.vertices() ) {
				//	mesh.property( vBoundaryId, vvvv ) = id;
				//}
			}
		}
	}

	//testPrint( mesh, fBoundaryId, vBoundaryId );

	fillDomainCorners( mesh, vBoundaryId, nBoundaries );
	
	//testPrint( mesh, fBoundaryId, vBoundaryId );

	// if boundaries touch, fuse them
	for( auto fh : mesh.faces() ) {
		std::vector<int> id;
		id.reserve( 3 );
		for( auto fv : fh.vertices() ) {
			id.push_back( mesh.property( vBoundaryId, fv ) );
		}
		std::sort( id.begin(), id.end() );

		if( id[0] == -1 )
			continue;
		if( id[0] == id[2] )
			continue;

		// use bfs to fuse boundaries
		std::queue<OpenMesh::SmartVertexHandle> q;
		for( auto fv : fh.vertices() ) {
			q.push( fv );
		}

		while( !q.empty() ) {
			auto v = q.front();
			q.pop();

			int& vId = mesh.property( vBoundaryId, v );
			if( vId == id[0] || vId == -1 )
				continue;

			vId = id[0];
			//for( auto vv : v.vertices() )
			//	q.push( vv );
			for( auto voh : v.outgoing_halfedges() ) {
				if( mesh.property( vBoundaryId, voh.next().to() ) != -1 || mesh.property( vBoundaryId, voh.opp().next().to() ) != -1 ) {
					q.push( voh.to() );
				}
			}
		}
	}

	fillDomainCorners( mesh, vBoundaryId, nBoundaries );

	//testPrint( mesh, fBoundaryId, vBoundaryId );

	std::vector<int> nBlossFaces(nBoundaries, 0);
	for( auto fh : mesh.faces() ) {
		bool blossFace = true;
		int id = -1;
		for( auto fv : fh.vertices() ) {
			if( mesh.property( vBoundaryId, fv ) == -1 ) {
				blossFace = false;
				break;
			} else if(id == -1) {
				id = mesh.property( vBoundaryId, fv );
			} else if( id != mesh.property( vBoundaryId, fv ) ) {
				blossFace = false;
				break;
			}
		}
		
		if( blossFace ) {
			mesh.property( fBoundaryId, fh ) = id;
			mesh.property( fBloss, fh ) = nBlossFaces[id]++;
		} else {
			mesh.property( fBloss, fh ) = -1;
			mesh.property( fBoundaryId, fh ) = -1;
		}
	}

	//testPrint( mesh, fBoundaryId, vBoundaryId );

	// fill small "interior islands"
	std::vector<std::vector<OpenMesh::SmartFaceHandle>> segments;
	OpenMesh::FPropHandleT<bool> isSeg;
	mesh.add_property( isSeg );
	for( auto fh : mesh.faces() ) {
		mesh.property( isSeg, fh ) = false;
	}
	for( auto fh : mesh.faces() ) {
		if( mesh.property( isSeg, fh ) || mesh.property(fBoundaryId, fh) != -1 )
			continue;

		// use bfs to find all elements of seg
		std::vector<OpenMesh::SmartFaceHandle> seg;
		std::queue<OpenMesh::SmartHalfedgeHandle> q;

		seg.push_back( fh );
		mesh.property( isSeg, fh ) = true;
		for( auto h : fh.halfedges() ) {
			q.push( h.opp() );
		}

		//LOG( INFO ) << "fh : " << fh.idx();
		while( !q.empty() ) {
			auto h = q.front();
			q.pop();
			auto f = h.face();
			//LOG( INFO ) << "f = " << f.idx() << "  |  p = " << mesh.point( h.from() );
			if( mesh.property( isSeg, f ) || mesh.property(fBoundaryId, f) != -1 || h.is_boundary() )
				continue;

			seg.push_back( f );
			mesh.property( isSeg, f ) = true;
			q.push( h.next().opp() );
			q.push( h.prev().opp() );
		}

		segments.push_back( seg );
	}
	// add small segments to boundary
	decltype(segments.size()) maxSegSize = 0;
	for( const auto& seg : segments ) {
		maxSegSize = std::max( seg.size(), maxSegSize );
	}

	for( const auto& seg : segments ) {
		if( seg.size() == maxSegSize )
			continue;

		// find dominant boundary type
		std::vector<int> nSegBoundaries( nBoundaries, 0 );
		for( const auto f : seg ) {
			for( auto h : f.halfedges() ) {
				auto neigh = h.opp().face();
				auto id = mesh.property( fBoundaryId, neigh );
				if( id != -1 ) {
					++nSegBoundaries[id];
				}
			}
		}
		int boundaryId = -1;
		int nSegs = 0;
		for( int i = 0; i < nSegBoundaries.size(); ++i ) {
			if( nSegBoundaries[i] > nSegs ) {
				nSegs = nSegBoundaries[i];
				boundaryId = i;
			}
		}

		// make all elements this id
		for( const auto f : seg ) {
			mesh.property( fBoundaryId, f ) = boundaryId;
			mesh.property( fBloss, f ) = nBlossFaces[boundaryId]++;
		}
	}
	mesh.remove_property( isSeg );

	// make number of faces even
	for( int i = 0; i < nBlossFaces.size(); ++i ) {
		if( ( nBlossFaces[i] & 1 ) == 0 )
			continue;

		for( auto eh : mesh.edges() ) {
			if( eh.is_boundary() )
				continue;

			auto fh1 = eh.h0().face();
			auto fh2 = eh.h0().opp().face();
			auto bloss1 = mesh.property( fBoundaryId, fh1 );
			auto bloss2 = mesh.property( fBoundaryId, fh2 );

			if( bloss1 == i && bloss2 == -1 ) {
				mesh.property( fBoundaryId, fh2 ) = i;
				mesh.property( fBloss, fh2 ) = nBlossFaces[i]++;
				break;
			}
			if( bloss1 == -1 && bloss2 == i ) {
				mesh.property( fBoundaryId, fh1 ) = i;
				mesh.property( fBloss, fh1 ) = nBlossFaces[i]++;
				break;
			}
		}
	}

	//testPrint( mesh, fBoundaryId, vBoundaryId );
	
	// add edges to perfect matching
	std::vector<std::vector<BlossomEdge>> blossEdges( nBoundaries );
	for( auto eh : mesh.edges() ) {
		auto heh = eh.h0();
		const auto fh1 = heh.face();
		const auto fh2 = heh.opp().face();
		if( eh.is_boundary() ) {
			if( heh.is_boundary() )
				heh = eh.h1();
			if( heh.opp().from().valence() < 4 )
				continue;
			blossEdges[mesh.property( fBoundaryId, heh.face() )].push_back( { heh.face(),heh.opp().prev().opp().face(), heh.edge(), 1000 } );
			continue;
		}

		const auto bloss1 = mesh.property( fBoundaryId, fh1 );
		const auto bloss2 = mesh.property( fBoundaryId, fh2 );

		if( bloss1 == -1 && bloss2 == -1 )
			continue;

		auto quadVhs = connect_triangles( heh );
		std::vector<PolyMesh::Point> points{ mesh.point( quadVhs[0] ), mesh.point( quadVhs[1] ), mesh.point( quadVhs[2] ), mesh.point( quadVhs[3] ) };

		real q = 1.f - QualityMetrics::condition( points );

		if( bloss1 == bloss2 ) {
			// internal
			blossEdges[bloss1].push_back( { fh1,fh2,eh,(int)( 100 * q ) } );
		} else {
			// add an external edge like on boundary
			if( bloss1 != -1 ) {
				auto hIt = heh.opp();
				do {
					hIt = hIt.prev().opp();
				} while( mesh.property( fBoundaryId, hIt.face() ) != bloss1 );
				if( fh1 != hIt.face() ) {
					blossEdges[bloss1].push_back( { fh1,hIt.face(), OpenMesh::SmartEdgeHandle(),1000 } );
				}
			}
			if( bloss2 != -1 ) {
				auto hIt = heh;
				do {
					hIt = hIt.prev().opp();
				} while( mesh.property( fBoundaryId, hIt.face() ) != bloss2 );
				if( hIt.face() != fh2 ) {
					blossEdges[bloss2].push_back( { hIt.face(),fh2, OpenMesh::SmartEdgeHandle(),1000 } );
				}
			}
		}
	}

	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();

	//////////////////////////////
	// perform perfect matching //
	PolyMesh hmesh;
	hmesh.reserve( mesh.n_vertices(), mesh.n_edges(), mesh.n_faces() );
	for( auto vh : mesh.vertices() ) {
		hmesh.add_vertex( mesh.point( vh ) );
	}
	// add interior triangles
	for( auto fh : mesh.faces() ) {
		if( mesh.property( fBoundaryId, fh ) == -1 ) {
			auto heh = fh.halfedge();
			hmesh.add_face( { mesh.vertex_handle( heh.from().idx() ), mesh.vertex_handle( heh.to().idx() ), mesh.vertex_handle( heh.next().to().idx() ) } );
		}
	}

	for( int i = 0; i < nBlossFaces.size(); ++i ) {
		LOG_ASSERT( (nBlossFaces[i] & 1) == 0 );
		
		PerfectMatching pm( nBlossFaces[i], blossEdges[i].size() );
		for( auto& bEdge : blossEdges[i] ) {
			const auto& f0 = bEdge.f0;
			const auto& f1 = bEdge.f1;
			const auto& id0 = mesh.property( fBloss, f0 );
			const auto& id1 = mesh.property( fBloss, f1 );
			bEdge.id = pm.AddEdge( id0, id1, bEdge.w );
		}
		pm.options.update_duals_before = true;
		pm.options.fractional_jumpstart = false;
		//pm.options.verbose = false;
		pm.Solve();

		// merge triangles
		for( const auto& bEdge : blossEdges[i] ) {
			int sol = pm.GetSolution( bEdge.id );

			switch( sol ) {
			case 0:
			{
				// do nothing
				break;
			}
			case 1:
			{
				// TODO handle bad boundary vertices directly here!!
				if( !bEdge.e.is_valid() ||		// edge between interior and boundary region
					bEdge.e.is_boundary() ) {	// external edge	

					auto heh = bEdge.f0.halfedge();
					hmesh.add_face( { mesh.vertex_handle( heh.from().idx() ), mesh.vertex_handle( heh.to().idx() ), mesh.vertex_handle( heh.next().to().idx() ) } );
					heh = bEdge.f1.halfedge();
					hmesh.add_face( { mesh.vertex_handle( heh.from().idx() ), mesh.vertex_handle( heh.to().idx() ), mesh.vertex_handle( heh.next().to().idx() ) } );
					continue;
				}

				// merge
				std::vector<PolyMesh::VertexHandle> vhs( 4 );
				vhs[0] = bEdge.e.h0().from();
				vhs[1] = bEdge.e.h0().opp().next().to();
				vhs[2] = bEdge.e.h0().to();
				vhs[3] = bEdge.e.h0().next().to();
				hmesh.add_face( vhs );
				break;
			}
			default:
				LOG( ERROR ) << " Edge solution = " << sol;
				break;
			}
		}
	}

	//OpenMesh::IO::write_mesh( hmesh, std::string(WORKING_DIRECTORY) + "output/ringQ.off" );

	mesh.remove_property( fBloss );
	mesh.remove_property( vBoundaryId );
	mesh.remove_property( fBoundaryId );

	return hmesh;
}

inline std::vector<PolyMesh::VertexHandle> blossom_quad::connect_triangles( const OpenMesh::SmartHalfedgeHandle& heh ) {
	LOG_ASSERT( !heh.edge().is_boundary() );

	std::vector<PolyMesh::VertexHandle> vhs(4);
	vhs[0] = heh.from();
	vhs[1] = heh.opp().next().to();
	vhs[2] = heh.to();
	vhs[3] = heh.next().to();

	return vhs;
}