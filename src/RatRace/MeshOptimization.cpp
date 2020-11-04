#include "MeshOptimization.h"

#include "QualityMetrics.h"	

#include <array>
#include <queue>

#include <glog/logging.h>

namespace MeshOptimization
{
	void getOneRing( PolyMesh& mesh, const TriMesh::VertexHandle& vh, std::vector<TriMesh::VertexHandle>& oneRingVh ) {
		auto heh = mesh.halfedge_handle( vh );
		auto heh_init = heh;
		do {
			oneRingVh.push_back( mesh.to_vertex_handle( heh ) );
			oneRingVh.push_back( mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) ) );
			heh = mesh.opposite_halfedge_handle( mesh.prev_halfedge_handle( heh ) );
		} while ( heh != heh_init );
		oneRingVh.push_back( oneRingVh[0] );
	}

	real jacobianQuadMetric( TriMesh::Point points[4] ) {
		
		TriMesh::Point e[4];
		real e_length_squared[4];
		for( auto i = 0; i < 4; ++i ) {
			auto j = ( i + 1 ) % 4;
			e[i] = points[j] - points[i];
			e_length_squared[i] = e[i].sqrnorm();
		}

		auto detJ1 = e[0][1] * e[3][0] - e[0][0] * e[3][1];
		auto l1 = e_length_squared[0] + e_length_squared[3];
		auto q1 = 2 * detJ1 / l1;

		auto detJ2 = e[1][1] * e[0][0] - e[1][0] * e[0][1];
		auto l2 = e_length_squared[1] + e_length_squared[0];
		auto q2 = 2 * detJ2 / l2;

		//auto detJ3 = e[2][1] * e[1][0] - e[2][0] * e[1][1];
		//auto l3 = e_length_squared[2] + e_length_squared[1];
		//auto q3 = 2 * detJ3 / l3;

		auto detJ4 = e[3][1] * e[2][0] - e[3][0] * e[2][1];
		auto l4 = e_length_squared[3] + e_length_squared[2];
		auto q4 = 2 * detJ4 / l4;

		auto qMin = std::min( { q1, q2, q4 } );
		auto detMin = std::min( { detJ1, detJ2, detJ4 } );

		if( detMin > 0 )
			return qMin;
		else
			return detMin;

	}

	real minQuality( PolyMesh& mesh, const TriMesh::VertexHandle& vh, const std::vector<TriMesh::VertexHandle>& oneRingVh, const TriMesh::Point& p ) {

		// compute quality of vh when moving it to position p
		real q = FLT_MAX;

		std::vector<TriMesh::Point> oneRing( oneRingVh.size() );

		for ( int i = 0; i < oneRingVh.size(); ++i ) {
			oneRing[i] = mesh.point( oneRingVh[i] );
		}

		for ( int k = 0; k < oneRing.size() - 2; k += 2 ) {
			PolyMesh::Point quad[4] = { p, oneRing[k], oneRing[k + 1], oneRing[k + 2] };
			//real quality = meanRatioMetricQuadWithUntangling( quad );
			real quality = jacobianQuadMetric( quad );
			q = std::min( q, quality );
		}

		return q;
	}

	void optimizeHierarchicalCPU( PolyMesh& mesh, const TriMesh::VertexHandle& vh ) {
		constexpr size_t N = 8;
		constexpr real affineFactor = 1.f / (real)( N - 1 );
		constexpr real grid_scale = 0.5f;
		constexpr int depth_max = 3;

		TriMesh::Point pOpt = mesh.point( vh );

		std::vector<TriMesh::VertexHandle> oneRingVh;

		// get oneRingVh
		getOneRing( mesh, vh, oneRingVh );

		if( oneRingVh.size() == 5 ) {
			// vertex has valence 2, just place it in the middle of its neighbors
			auto p1 = mesh.point( oneRingVh[0] );
			auto p2 = mesh.point( oneRingVh[2] );
			auto pNew = 0.5 * ( p1 + p2 );
			mesh.set_point( vh, pNew );
			return;
		}

		// compute grid-size
		real aabbSize = -FLT_MAX;
		for ( auto vh : oneRingVh ) {
			auto d = ( pOpt - mesh.point( vh ) ).length();
			aabbSize = std::max( aabbSize, d );
		}

		aabbSize = grid_scale * aabbSize;

		auto depth_scale = grid_scale;

		auto qualityBefore = minQuality( mesh, vh, oneRingVh, pOpt );
		auto pBefore = pOpt;

		// compute quality
		for ( int depth = 0; depth < depth_max; ++depth ) {
			std::array<real, N*N> q;

			real uMax, uMin, vMax, vMin;
			uMax = +depth_scale * aabbSize;
			uMin = -depth_scale * aabbSize;
			vMax = +depth_scale * aabbSize;
			vMin = -depth_scale * aabbSize;

			//#pragma omp parallel for
			for ( int i = 0; i < N; ++i ) {

				for ( int j = 0; j < N; ++j ) {

					real u = pOpt[0] + affineFactor * ( i * uMin + ( N - 1 - i ) * uMax );
					real v = pOpt[1] + affineFactor * ( j * vMin + ( N - 1 - j ) * vMax );

					// find minimal quality at node (i,j)
					q[i * N + j] = minQuality( mesh, vh, oneRingVh, { u,v,0 } );
				}

			}

			// find max of q
			size_t iOpt = -1;
			size_t jOpt = -1;
			real qOpt = -FLT_MAX;

			for ( size_t i = 0; i < N; ++i ) {
				for ( size_t j = 0; j < N; ++j ) {
					if ( q[i * N + j] > qOpt ) {
						iOpt = i;
						jOpt = j;
						qOpt = q[i * N + j];
					}
				}
			}

			real u = pOpt[0] + affineFactor * ( iOpt * uMin + ( N - 1 - iOpt ) * uMax );
			real v = pOpt[1] + affineFactor * ( jOpt * vMin + ( N - 1 - jOpt ) * vMax );

			pOpt = { u,v,0 };

			//depth dependent scaling factor
			depth_scale = depth_scale * ( 2.f / ( N - 1 ) );
		}

		real qualityAfter = minQuality( mesh, vh, oneRingVh, pOpt );
		if ( qualityAfter > qualityBefore ) {
			mesh.set_point( vh, pOpt );
		}

	}

	void discreteMeshOptimizationCPU( PolyMesh& mesh, const int n_iter ) {

		for ( auto i = 0; i < n_iter; ++i ) {

			for ( auto vh_it = mesh.vertices_begin(); vh_it != mesh.vertices_end(); ++vh_it ) {
				if ( mesh.is_boundary( *vh_it ) ) continue;
				optimizeHierarchicalCPU( mesh, *vh_it );
			}

		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	void laplaceSmoothing( PolySurfMesh& psm, const int n_iter ) {
		for( auto iter = 0; iter < n_iter; ++iter ) {
			for( auto vh : psm.vertices() ) {
				if( psm.is_boundary( vh ) )
					continue;
				std::vector<TriMesh::VertexHandle> oneRingVh;
				for( auto vv : psm.vv_range( vh ) )
					oneRingVh.push_back( vv );
				oneRingVh.push_back( vh );

				
				PolyMesh::Point uvwMean( 0, 0, 0 );
				for( auto vv : oneRingVh ) {
					auto xyz = psm.point( vv );
					uvwMean += xyz;
				}
				uvwMean /= oneRingVh.size();
				
				psm.set_point( vh, uvwMean );
			}
		}

	}

	
	PolyMesh removeInnerValence2Vertices( PolyMesh& mesh ) {
		PolyMesh m = mesh;

		m.request_vertex_status();
		m.request_edge_status();
		m.request_face_status();

		bool foundOne = false;
		do {
			foundOne = false;

			for ( auto vh : m.vertices() ) {
				if ( m.is_boundary( vh ) ) continue;
				auto valence = std::distance( m.vv_begin( vh ), m.vv_end( vh ) );
				if ( valence != 2 ) continue;

				std::vector<PolyMesh::VertexHandle> neighs;
				getOneRing( m, vh, neighs );
				neighs.pop_back();

				LOG_ASSERT( neighs.size() == 4 );

				m.delete_vertex( vh, false );
				
				m.add_face( neighs );

				foundOne = true;
				break;
			}

		} while ( foundOne );

		m.garbage_collection();

		return m;
	}

	/*Do not use!!! Not working properly!!!*/
	PolyMesh removeBoundaryValence2Vertices( PolyMesh& mesh, real qMin ) {
		PolyMesh m = mesh;

		m.request_vertex_status();
		m.request_edge_status();
		m.request_face_status();

		bool foundOne = false;
		do {
			foundOne = false;

			for( auto vh : m.vertices() ) {
				if( !m.is_boundary( vh ) ) continue;
				auto valence = std::distance( m.vv_begin( vh ), m.vv_end( vh ) );
				if( valence != 2 ) continue;

				PolyMesh::HalfedgeHandle heh;
				for( auto voh : m.voh_range( vh ) ) {
					if( !m.is_boundary( voh ) ) {
						heh = voh;
						break;
					}
				}

				PolyMesh::FaceHandle fh = m.face_handle( heh );

				auto q = QualityMetrics::condition( m, fh );
				if( q > qMin ) continue;	// quad has sufficient quality and does not need to be removed

				auto vhNext = m.to_vertex_handle( heh );
				auto vhPrev = m.from_vertex_handle( m.prev_halfedge_handle( heh ) );
				
				auto hehIt = heh;
				hehIt = m.next_halfedge_handle( hehIt );
				hehIt = m.opposite_halfedge_handle( hehIt );
				hehIt = m.next_halfedge_handle( hehIt );
				std::vector<PolyMesh::VertexHandle> verticesNewFh;
				verticesNewFh.reserve( 4 );
				verticesNewFh.push_back( m.to_vertex_handle( hehIt ) );
				hehIt = m.next_halfedge_handle( hehIt );
				verticesNewFh.push_back( m.to_vertex_handle( hehIt ) );
				hehIt = m.next_halfedge_handle( hehIt );
				verticesNewFh.push_back( m.to_vertex_handle( hehIt ) );
				verticesNewFh.push_back( vhPrev );
								
				m.set_point( vhPrev, m.point( vh ) );

				m.delete_vertex( vh );
				m.delete_vertex( vhNext );

				m.add_face( verticesNewFh );

				foundOne = true;
				break;
			}

		} while( foundOne );

		m.garbage_collection();

		return m;
	}

}