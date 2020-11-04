#pragma once

#include <glog/logging.h>

#include <vector>
#include <array>

#include "MeshHeader.h"

namespace QualityMetrics
{
	inline real condition( const std::vector<PolyMesh::Point>& points ) {
		const auto ps = points.size();

		std::vector<PolyMesh::Point> edges( ps, { 0,0,0 } );
		std::vector<real> edgeSqrNorm( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto j = ( i + 1 ) % ps;
			edges[i] = points[j] - points[i];
			edges[i][2] = 0;
			edgeSqrNorm[i] = edges[i].sqrnorm();
		}

		auto determinant = [](const PolyMesh::Point& p1, const PolyMesh::Point& p2 ) {
			return ( p1 % p2 )[2];
		};

		std::vector<real> detJ( ps, -1 );
		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			detJ[i] = determinant( edges[i], -edges[il] );
		}


		std::vector<real> c( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			c[i] = detJ[i] / ( edgeSqrNorm[i] + edgeSqrNorm[il] );
		}

		auto cMin = std::min_element( c.begin(), c.end() );

		return 2 * *cMin;
	}

	inline real condition( PolyMesh& mesh, const PolyMesh::FaceHandle& fh ) {
		std::vector<PolyMesh::Point> p;
		p.reserve( 4 );
		for( auto fv : mesh.fv_range( fh ) ) {
			p.push_back( mesh.point( fv ) );
		}

		return condition( p );
	}
	
	// condition number for vertex '1'
	inline real conditionSurfVertex( const std::vector<PolyMesh::Point>& points ) {
		LOG_ASSERT( points.size() == 3 );

		//			 2
		// 			 |
		//			 e1
		//	 		 |
		// 0---e2----1

		PolyMesh::Point e1 = points[2] - points[1];
		PolyMesh::Point e2 = points[0] - points[1];

		auto l1 = e1.sqrnorm();
		auto l2 = e2.sqrnorm();

		PolyMesh::Point n = { 0,0,1 };
		auto detJ = ( e1 % e2 ) | n;

		if( detJ < 0 ) {
			return detJ;
		} else {
			auto c = 2 * detJ / ( l1 + l2 );
			return c;
		}
	}

	inline real conditionSurf( const std::vector<PolyMesh::Point>& points ) {
		const auto ps = points.size();

		std::vector<PolyMesh::Point> edges( ps, { 0,0,0 } );
		std::vector<real> edgeSqrNorm( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto j = ( i + 1 ) % ps;
			edges[i] = points[j] - points[i];
			edgeSqrNorm[i] = edges[i].sqrnorm();
		}

		auto determinant = []( const PolyMesh::Point& p1, const PolyMesh::Point& p2 ) {
			return p1[0] * p2[1] - p1[1] * p2[0];
		};

		std::vector<real> detJ( ps, -1 );
		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			detJ[i] = determinant( edges[i], -edges[il] );
		}


		std::vector<real> c( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			c[i] = detJ[i] / ( edgeSqrNorm[i] + edgeSqrNorm[il] );
		}

		auto detMin = std::min_element( detJ.begin(), detJ.end() );

		if( *detMin > 0 ) {
			auto cMin = std::min_element( c.begin(), c.end() );
			return 2 * *cMin;
		} else {
			return *detMin;
		}

	}

	inline real meanConditionSurf( const std::vector<PolyMesh::Point>& points ) {
		const auto ps = points.size();

		std::vector<PolyMesh::Point> edges( ps, { 0,0,0 } );
		std::vector<real> edgeSqrNorm( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto j = ( i + 1 ) % ps;
			edges[i] = points[j] - points[i];
			edgeSqrNorm[i] = edges[i].sqrnorm();
		}

		auto determinant = []( const PolyMesh::Point& p1, const PolyMesh::Point& p2 ) {
			return p1[0] * p2[1] - p1[1] * p2[0];
		};

		std::vector<real> detJ( ps, -1 );
		for( auto i = 0; i < ps; ++i ) {
			auto il = ( i - 1 + ps ) % ps;
			detJ[i] = determinant( edges[i], -edges[il] );
		}


		std::vector<real> c( ps, -1 );

		for( auto i = 0; i < ps; ++i ) {
			auto ir = ( i + 1 ) % ps;
			auto il = ( i - 1 + ps ) % ps;

			c[i] = detJ[i] / ( edgeSqrNorm[ir] + edgeSqrNorm[il] );
		}

		real cMean = 0;
		for( auto e : c ) {
			cMean += e;
		}
		return cMean / c.size();

	}

	inline real conditionSurf( PolyMesh& mesh, TriMesh& initMesh, const PolyMesh::FaceHandle& fh ) {
		std::vector<PolyMesh::Point> p;
		p.reserve( 4 );
		for( auto fv : mesh.fv_range( fh ) ) {
			p.push_back( mesh.point( fv ) );
		}

		return conditionSurf( p );
	}

	inline real conditionSurf( PolyMesh& mesh, const PolyMesh::FaceHandle& fh ) {
		std::vector<PolyMesh::Point> p;
		p.reserve( 4 );
		for( auto fv : mesh.fv_range( fh ) ) {
			p.push_back( mesh.point( fv ) );
		}

		return conditionSurf( p );
	}

	inline auto minConditionSurf( PolyMesh& mesh, TriMesh& initMesh ) {
		real qMin = FLT_MAX;
		real qMean = 0;

		for( auto fh : mesh.faces() ) {
			auto q = conditionSurf( mesh, initMesh, fh );
			qMin = std::min( q, qMin );
			qMean += q;

		}
		qMean /= mesh.n_faces();
		return std::make_tuple( qMin, qMean );
	}

	inline auto minConditionSurf( PolyMesh& mesh ) {
		real qMin = FLT_MAX;
		real qMean = 0;

		for( auto fh : mesh.faces() ) {
			auto q = conditionSurf( mesh, fh );
			qMin = std::min( q, qMin );
			qMean += q;

		}
		qMean /= mesh.n_faces();
		return std::make_tuple( qMin, qMean );
	}
}