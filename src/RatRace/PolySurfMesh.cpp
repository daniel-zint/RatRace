#include "PolySurfMesh.h"
#include <algorithm>

void PolySurfMesh::edgeSwapNext( const EdgeHandle& eh ) {
	LOG_ASSERT( !is_boundary( eh ) );
	auto heh = halfedge_handle( eh, 0 );
	LOG_ASSERT( valence( from_vertex_handle( heh ) ) > 2 );
	LOG_ASSERT( valence( to_vertex_handle( heh ) ) > 2 );

	//	1 --- 0 --- 5
	//	|     ^     |
	//	|     |     |
	//	2 --- 3 --- 4

	auto hehInit = heh;
	std::vector<PolyMesh::VertexHandle> vhs;
	vhs.push_back( to_vertex_handle( heh ) );	// 0
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 1
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 2
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 3
	heh = opposite_halfedge_handle( hehInit );
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 4
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 5
	heh = next_halfedge_handle( heh );

	delete_edge( eh, false );
	auto newFh1 = add_face( { vhs[4], vhs[5], vhs[0], vhs[1] } );
	auto newFh2 = add_face( { vhs[1], vhs[2], vhs[3], vhs[4] } );
}

void PolySurfMesh::edgeSwapPrev( const EdgeHandle& eh ) {
	LOG_ASSERT( !is_boundary( eh ) );
	auto heh = halfedge_handle( eh, 0 );
	LOG_ASSERT( valence( from_vertex_handle( heh ) ) > 2 );
	LOG_ASSERT( valence( to_vertex_handle( heh ) ) > 2 );

	//	1 --- 0 --- 5
	//	|     ^     |
	//	|     |     |
	//	2 --- 3 --- 4

	auto hehInit = heh;
	std::vector<PolyMesh::VertexHandle> vhs;
	vhs.push_back( to_vertex_handle( heh ) );	// 0
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 1
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 2
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 3
	heh = opposite_halfedge_handle( hehInit );
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 4
	heh = next_halfedge_handle( heh );
	vhs.push_back( to_vertex_handle( heh ) );	// 5
	heh = next_halfedge_handle( heh );

	delete_edge( eh, false );
	auto newFh1 = add_face( { vhs[5], vhs[0], vhs[1], vhs[2] } );
	auto newFh2 = add_face( { vhs[2], vhs[3], vhs[4], vhs[5] } );
}

OpenMesh::SmartVertexHandle PolySurfMesh::vertexSplit( const HalfedgeHandle & heh ) {
	OpenMesh::SmartHalfedgeHandle shh( heh.idx(), this );
	LOG_ASSERT( !shh.edge().is_boundary() );

	auto vh = shh.from();
	auto fh1 = shh.face();
	auto fh2 = shh.opp().face();
	auto pNew = 0.5f * ( point( vh ) + point( shh.to() ) );

	decltype(vh) vhNew = add_vertex( pNew );

	std::vector<OpenMesh::SmartVertexHandle> fv1, fv2, fv3;
	for( auto fv : fv_range( fh1 ) ) {
		fv1.push_back( fv );
	}
	for( auto fv : fv_range( fh2 ) ) {
		fv2.push_back( fv );
	}
	std::replace( fv1.begin(), fv1.end(), vh, vhNew );
	std::replace( fv2.begin(), fv2.end(), vh, vhNew );

	// build new face
	fv3.push_back( vhNew );
	fv3.push_back( shh.prev().from() );
	fv3.push_back( vh );
	fv3.push_back( shh.opp().next().to() );

	// TODO replace this by delete_edge
	delete_face( fh1, false );
	delete_face( fh2, false );

	add_face( fv1 );
	add_face( fv2 );
	add_face( fv3 );

	return vhNew;
}

void PolySurfMesh::edge3Split( const HalfedgeHandle& heh ) {
	const auto v0 = from_vertex_handle( heh );
	const auto v3 = to_vertex_handle( heh );
	const auto p0 = point( v0 );
	const auto p3 = point( v3 );
	const auto p1 = ( 2.f / 3.f ) * p0 + ( 1.f / 3.f ) * p3;
	const auto p2 = ( 1.f / 3.f ) * p0 + ( 2.f / 3.f ) * p3;
	const decltype(v0) v1 = add_vertex( p1 );
	const decltype(v0) v2 = add_vertex( p2 );

	const auto hehOpp = opposite_halfedge_handle( heh );

	if( !is_boundary( heh ) ) {
		const auto fh = face_handle( heh );
		const auto v4 = to_vertex_handle( next_halfedge_handle( heh ) );
		const auto v5 = from_vertex_handle( prev_halfedge_handle( heh ) );
		const auto p4 = point( v4 );
		const auto p5 = point( v5 );
		const auto p6 = 0.5 * ( p1 + p4 );
		const auto p7 = 0.5 * ( p2 + p5 );
		const auto v6 = add_vertex( p6 );
		const auto v7 = add_vertex( p7 );

		delete_face( fh, false );
		add_face( { v0,v1,v7,v5 } );
		add_face( { v1,v2,v6,v7 } );
		add_face( { v2,v3,v4,v6 } );
		add_face( { v5,v7,v6,v4 } );
	}
	if( !is_boundary( hehOpp ) ) {
		const auto fhOpp = face_handle( hehOpp );
		const auto v8 = to_vertex_handle( next_halfedge_handle( hehOpp ) );
		const auto v9 = from_vertex_handle( prev_halfedge_handle( hehOpp ) );
		const auto p8 = point( v8 );
		const auto p9 = point( v9 );
		const auto p10 = 0.5 * ( p2 + p8 );
		const auto p11 = 0.5 * ( p1 + p9 );
		const auto v10 = add_vertex( p10 );
		const auto v11 = add_vertex( p11 );

		delete_face( fhOpp, false );
		add_face( { v0,v8,v10,v1 } );
		add_face( { v1,v10,v11,v2 } );
		add_face( { v2,v11,v9,v3 } );
		add_face( { v8,v9,v11,v10 } );
	}
}

void PolySurfMesh::diagonalCollapse( const HalfedgeHandle & heh ) {
	LOG_ASSERT( !is_boundary( edge_handle( heh ) ) );
	auto fh = face_handle( heh );
	auto vh1 = from_vertex_handle( heh );
	auto vh2 = to_vertex_handle( next_halfedge_handle( heh ) );

	std::set<FaceHandle> neighs;
	for( auto fv : fv_range( fh ) ) {
		for( auto vf : vf_range( fv ) ) {
			neighs.insert( vf );
		}
	}
	neighs.erase( fh );

	std::vector<std::vector<VertexHandle>> neighVhs;
	for( auto n : neighs ) {
		std::vector<VertexHandle> vhs;
		for( auto fv : fv_range( n ) ) {
			vhs.push_back( fv );
		}
		std::replace( vhs.begin(), vhs.end(), vh1, vh2 );
		neighVhs.push_back( vhs );
	}

	delete_face( fh, false );
	for( auto n : neighs ) {
		delete_face( n, false );
	}
	delete_vertex( vh1, false );
	for( auto vhs : neighVhs ) {
		add_face( vhs );
	}

}

void PolySurfMesh::diagonalCollapseDirect( const HalfedgeHandle& heh ) {
	auto fh = face_handle( heh );
	
	std::vector<PolySurfMesh::HalfedgeHandle> h( 4 ), o( 4 );
	h[0] = heh;
	h[1] = next_halfedge_handle( h[0] );
	h[2] = next_halfedge_handle( h[1] );
	h[3] = next_halfedge_handle( h[2] );
	LOG_ASSERT( next_halfedge_handle( h[3] ) == h[0] ) << " face is not a quad";
	
	o[0] = opposite_halfedge_handle( h[0] );
	o[1] = opposite_halfedge_handle( h[1] );
	o[2] = opposite_halfedge_handle( h[2] );
	o[3] = opposite_halfedge_handle( h[3] );

	std::vector<PolySurfMesh::EdgeHandle> e( 4 );
	e[0] = edge_handle( h[0] );
	e[1] = edge_handle( h[1] );
	e[2] = edge_handle( h[2] );
	e[3] = edge_handle( h[3] );

	auto v0 = from_vertex_handle( h[0] );
	auto v1 = from_vertex_handle( h[1] );
	auto v2 = from_vertex_handle( h[2] );
	auto v3 = from_vertex_handle( h[3] );

	int nVal2 = (valence(v0) == 2) + ( valence( v1 ) == 2 ) + ( valence( v2 ) == 2 ) + ( valence( v3 ) == 2 );
	LOG_ASSERT( nVal2 < 2 );

	auto o0prev = prev_halfedge_handle( o[0] );
	auto o0next = next_halfedge_handle( o[0] );
	auto o3prev = prev_halfedge_handle( o[3] );
	auto o3next = next_halfedge_handle( o[3] );

	auto fo0 = face_handle( o[0] );
	auto fo3 = face_handle( o[3] );

	if( fo0.is_valid() && halfedge_handle( fo0 ) == o[0] ) {
		set_halfedge_handle( fo0, h[1] );
	}
	if( fo3.is_valid() && halfedge_handle( fo3 ) == o[3] ) {
		set_halfedge_handle( fo3, h[2] );
	}

	for( auto vih : vih_range( v0 ) ) {
		set_vertex_handle( vih, v2 );
	}

	set_next_halfedge_handle( o0prev, h[1] );
	if( fo0 != fo3 ) {
		set_next_halfedge_handle( h[1], o0next );
	}
	
	set_next_halfedge_handle( h[2], o3next );
	if( fo0 != fo3 ) {
		set_next_halfedge_handle( o3prev, h[2] );
	}
	
	if( halfedge_handle( v3 ) == h[3] ) {
		set_halfedge_handle( v3, o[2] );
	}
	if( halfedge_handle( v1 ) == o[0] ) {
		set_halfedge_handle( v1, h[1] );
	}

	if( is_boundary( v0 ) ) {
		// find boundary heh
		for( auto voh : voh_range( v2 ) ) {
			if( voh.is_boundary() ) {
				set_halfedge_handle( v2, voh );
				break;
			}
		}
	}

	set_face_handle( h[1], fo0 );
	set_face_handle( h[2], fo3 );
		

	set_isolated( v0 );
		
	// delete stuff
	status( v0 ).set_deleted( true );
	status( fh ).set_deleted( true );
	status( e[0] ).set_deleted( true );
	status( e[3] ).set_deleted( true );
	if( has_halfedge_status() ) {
		status( h[0] ).set_deleted( true );
		status( h[3] ).set_deleted( true );
		status( o[0] ).set_deleted( true );
		status( o[3] ).set_deleted( true );
	}
	if( is_boundary( h[1] ) && is_boundary( o[1] ) ) {
		auto n = next_halfedge_handle( o[0] );
		auto p = prev_halfedge_handle( o[1] );
		set_next_halfedge_handle( p, n );
		if( halfedge_handle( v2 ) == o[1] ) {
			set_halfedge_handle( v2, n );
		}
		delete_edge( e[1], false );
		set_isolated( v1 );
		status( v1 ).set_deleted( true );
	}
	if( is_boundary( h[2] ) && is_boundary( o[2] ) ) {
		auto n = next_halfedge_handle( o[2] );
		auto p = prev_halfedge_handle( h[2] );
		set_next_halfedge_handle( p, n );
		delete_edge( e[2], false );
		set_isolated( v3 );
		status( v3 ).set_deleted( true );
	}
}

void PolySurfMesh::halfedgeCollapse( const HalfedgeHandle& heh ) {
	// collect vertex handles around vh_from starting at vh_to
	OpenMesh::SmartHalfedgeHandle shh( heh.idx(), this );
	LOG_ASSERT( !shh.edge().is_boundary() );

	auto vh1 = shh.from();
	LOG_ASSERT( !vh1.is_boundary() );
	auto vh2 = shh.to();

	std::vector<FaceHandle> fhs;
	for( auto vf : vf_range( vh1 ) ) {
		fhs.push_back( vf );
	}
	std::vector<VertexHandle> vhs;
	vhs.reserve( fhs.size() );
	auto shhIt = shh;
	do {
		vhs.push_back( shhIt.to() );
		shhIt = shhIt.next();
		vhs.push_back( shhIt.to() );
		shhIt = shhIt.next().next().opp();
	} while( shhIt != heh );

	for( auto fh : fhs ) {
		delete_face( fh, false );
	}
	delete_vertex( vh1 );

	add_face( { vhs[0], vhs[1], vhs[2], vhs[3] } );
	for( int i = 3; i < vhs.size() - 1; i += 2 ) {
		add_face( { vhs[i], vhs[i + 1], vhs[i + 2], vhs[0] } );
	}
}

void PolySurfMesh::vertexRotate( const VertexHandle& vh ) {
	OpenMesh::SmartVertexHandle svh( vh.idx(), this );
	LOG_ASSERT( !svh.is_boundary() );

	std::vector<OpenMesh::SmartVertexHandle> ringA, ringB;
	for( auto voh : svh.outgoing_halfedges() ) {
		ringA.push_back( voh.to() );
		ringB.push_back( voh.next().to() );
	}

	// perform vertex rotate
	std::vector<FaceHandle> fhs;
	for( auto vf : vf_range( vh ) ) {
		fhs.push_back( vf );
	}
	for( auto fh : fhs ) {
		delete_face( fh, false );
	}
	for( int i = 0; i < ringA.size(); ++i ) {
		auto vh1 = ringB[( i + 1 ) % ringA.size()];
		auto vh2 = ringA[i];
		auto vh3 = ringB[i];
		add_face( { vh, vh1, vh2, vh3 } );
	}
}

inline PolyMesh mergeDoubleVertices( PolyMesh& pm ) {
	PolyMesh m;

	// search for vertices with the same point
	std::map<PolyMesh::VertexHandle, PolyMesh::VertexHandle> vhMap;

	// double vertices can only appear on boundaries
	for( auto vh : pm.vertices() ) {
		if( pm.is_boundary( vh ) )
			continue;
		vhMap[vh] = vh;
	}

	for( auto vh : pm.vertices() ) {
		if( !pm.is_boundary( vh ) )
			continue;

		auto p = pm.point( vh );
		std::vector<PolyMesh::VertexHandle> dvh;

		for( auto vv : pm.vertices() ) {
			if( !pm.is_boundary( vv ) )
				continue;

			auto p2 = pm.point( vv );
			if( p == p2 ) {
				dvh.push_back( vv );
			}
		}

		std::sort( dvh.begin(), dvh.end() );
		for( auto d : dvh ) {
			vhMap[d] = dvh[0];
		}
	}

	// find maximal index in vhMap[i].second
	auto idxMax = -1;
	for( auto vhm : vhMap ) {
		idxMax = std::max( idxMax, vhm.second.idx() );
	}

	for( auto i = 0; i <= idxMax; ++i ) {
		m.add_vertex( pm.point( pm.vertex_handle( i ) ) );
	}

	for( auto fh : pm.faces() ) {
		std::vector<PolyMesh::VertexHandle> vhs;
		vhs.reserve( 4 );
		for( auto fv : pm.fv_range( fh ) ) {
			vhs.push_back( vhMap[fv] );
		}
		m.add_face( vhs );
	}

	return m;
}