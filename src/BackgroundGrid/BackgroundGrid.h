#pragma once

#include <string>
#include <assert.h>
#include <vector>
#include <list>
#include <experimental/filesystem>

#include <glog/logging.h>
#include <functional>

// Open Mesh
#include "../RatRace/precision.h"
#include "../RatRace/MeshHeader.h"

class BackgroundGrid
{
protected:
	real xMin_, xMax_, yMin_, yMax_; // AABB of domain

	size_t Nx_, Ny_;	// Number of grid points
	real hx_, hy_;		// Stepsize

	size_t vSize_;		// size of v
	std::vector<real> v_;			// Vector with average triangle size per vertex
	std::vector<bool> isDomain_;

	std::function<real( TriMesh&, const TriMesh::VertexHandle )> valuePerVertex_;

public:
	BackgroundGrid() {};
	BackgroundGrid( const std::experimental::filesystem::path& filename );
	BackgroundGrid( TriMesh& mesh, decltype( valuePerVertex_ ) valuePerVertex, const size_t& Nx = 0, const size_t& Ny = 0, const real scalingfactor = 1.0 );
	~BackgroundGrid();

	// Access functions
	real getSize( real x, real y ) const;		// get size-value at a certain position
	real getSize( const TriMesh::Point& p ) const;
	real getSize( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;
	real getSize( PolyMesh& mesh, const TriMesh::EdgeHandle eh ) const;
	real getSize( TriMesh& mesh, const TriMesh::EdgeHandle eh ) const;
	real getMinSize( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const;
	real getNeighMinSize( TriMesh& mesh, const TriMesh::EdgeHandle eh ) const;
	size_t Nx() const;
	size_t Ny() const;
	real hx() const;
	real hy() const;
	real xMin() const;
	real xMax() const;
	real yMin() const;
	real yMax() const;

	// Print functions
	void toVTK( const std::experimental::filesystem::path& name, const bool& onlyInterior = true, const bool& isBinary = true ) const;
	// DEPRECATED!! Will be deleted soon!
	void compare( const BackgroundGrid& sg, const std::experimental::filesystem::path& name ) const;

	// Access functions
	real& operator()( size_t i, size_t j );
	real operator()( size_t i, size_t j ) const;
	real& operator()( size_t i );
	real operator()( size_t i ) const;

	size_t lexIdx( const size_t &i, const size_t &j ) const;

	void applyToMeshVertices( TriMesh& mesh );

	void writeBinary( const std::experimental::filesystem::path& name ) const;
	void readBinary( const std::experimental::filesystem::path& name );

private:
	std::tuple<real, real, real, real> findAABB( TriMesh& mesh );
	size_t coords2LexInd( real x, real y ); // calculates the grid index that is in x and y value smaller than the values for a given coordinate (bottom left corner of the cell)
	size_t x2col( real x ) const;
	size_t y2row( real y ) const;
	size_t lexInd2col( const size_t& j ) const;
	size_t lexInd2row( const size_t& i ) const;

	real x( const size_t &i ) const;
	real y( const size_t &j ) const;

	void initGrid( TriMesh& mesh, const size_t& Nx, const size_t& Ny, const real scalingfactor );
	void fillGrid( TriMesh& mesh, const std::vector<real>& vertex_avg );
	real findMinEdgeLength( TriMesh& mesh );
	void collectVertexValues( TriMesh mesh, std::vector<real>& vertex_avg );

	void setOuterVals();
	void setGridSize( TriMesh& mesh, real minlength );
	void setGridSpecs( TriMesh& mesh );
};


inline real BackgroundGrid::getSize( real x, real y ) const {
	// this function was created according to the wikipedia article: "Bilinear interpolation"

	size_t i1 = this->x2col( x );
	size_t i2 = i1 + 1;
	if ( i2 == Nx_ ) {
		i1--;
		i2--;
	}
	size_t j1 = this->y2row( y );
	size_t j2 = j1 + 1;
	if ( j2 == Ny_ ) {
		j1--;
		j2--;
	}

	const real x1 = this->x( i1 );
	const real y1 = this->y( j1 );
	const real x2 = this->x( i2 );
	const real y2 = this->y( j2 );

	const real pref = 1.f / ( ( x2 - x1 ) * ( y2 - y1 ) );
	const real a = x2 - x;
	const real b = x - x1;
	const real c = ( *this )( i1, j1 );
	const real d = ( *this )( i1, j2 );
	const real e = ( *this )( i2, j1 );
	const real f = ( *this )( i2, j2 );
	const real g = y2 - y;
	const real h = y - y1;

	return pref * ( a * ( c * g + d * h ) + b * ( e * g + f * h ) );
}
inline real BackgroundGrid::getSize( const TriMesh::Point& p ) const {
	return getSize( p[0], p[1] );
}
inline real BackgroundGrid::getSize( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
	real val = 0;
	const size_t n_steps = 10;
	for ( size_t i = 0; i < n_steps; ++i ) {
		real x = (real)( n_steps - 1 - i ) / (real)( n_steps - 1 ) * p1[0] + (real)i / (real)( n_steps - 1 ) * p2[0];
		real y = (real)( n_steps - 1 - i ) / (real)( n_steps - 1 ) * p1[1] + (real)i / (real)( n_steps - 1 ) * p2[1];

		val += this->getSize( x, y );
	}
	val /= (real)n_steps;

	return val;
}
inline real BackgroundGrid::getSize( PolyMesh& mesh, const TriMesh::EdgeHandle eh ) const {
	const TriMesh::HalfedgeHandle heh = mesh.halfedge_handle( eh, 0 );

	const TriMesh::VertexHandle vh1 = mesh.from_vertex_handle( heh );
	const TriMesh::VertexHandle vh2 = mesh.to_vertex_handle( heh );

	const TriMesh::Point p1 = mesh.point( vh1 );
	const TriMesh::Point p2 = mesh.point( vh2 );

	return this->getSize( p1, p2 );
}
inline real BackgroundGrid::getSize( TriMesh& mesh, const TriMesh::EdgeHandle eh ) const {
	const TriMesh::HalfedgeHandle heh = mesh.halfedge_handle( eh, 0 );

	const TriMesh::VertexHandle vh1 = mesh.from_vertex_handle( heh );
	const TriMesh::VertexHandle vh2 = mesh.to_vertex_handle( heh );

	const TriMesh::Point p1 = mesh.point( vh1 );
	const TriMesh::Point p2 = mesh.point( vh2 );

	return this->getSize( p1, p2 );
}
inline real BackgroundGrid::getMinSize( const TriMesh::Point& p1, const TriMesh::Point& p2 ) const {
	real val = FLT_MAX;
	const size_t n_steps = 10;
	for( size_t i = 0; i < n_steps; ++i ) {
		real x = (real)( n_steps - 1 - i ) / (real)( n_steps - 1 ) * p1[0] + (real)i / (real)( n_steps - 1 ) * p2[0];
		real y = (real)( n_steps - 1 - i ) / (real)( n_steps - 1 ) * p1[1] + (real)i / (real)( n_steps - 1 ) * p2[1];

		val = std::min( val, this->getSize( x, y ) );
	}

	return val;
}
inline real BackgroundGrid::getNeighMinSize( TriMesh& mesh, const TriMesh::EdgeHandle eh ) const {
	auto heh = mesh.halfedge_handle( eh, 0 );
	auto vh1 = mesh.from_vertex_handle( heh );
	auto vh2 = mesh.to_vertex_handle( heh );
	auto p1 = mesh.point( vh1 );
	auto p2 = mesh.point( vh2 );

	auto s = getSize( p1, p2 );

	if( !mesh.is_boundary( heh ) ) {
		auto vh3 = mesh.to_vertex_handle( mesh.next_halfedge_handle( heh ) );
		auto p3 = mesh.point( vh3 );
		s = std::min( s, getSize( p1, p3 ) );
		s = std::min( s, getSize( p2, p3 ) );
	}
	if( !mesh.is_boundary( mesh.opposite_halfedge_handle( heh ) ) ) {
		auto vh4 = mesh.to_vertex_handle( mesh.next_halfedge_handle( mesh.opposite_halfedge_handle( heh ) ) );
		auto p4 = mesh.point( vh4 );
		s = std::min( s, getSize( p1, p4 ) );
		s = std::min( s, getSize( p2, p4 ) );
	}

	return s;
}

inline size_t BackgroundGrid::Nx() const {
	return Nx_;
}
inline size_t BackgroundGrid::Ny() const {
	return Ny_;
}
inline real BackgroundGrid::hx() const {
	return hx_;
}
inline real BackgroundGrid::hy() const {
	return hy_;
}
inline real BackgroundGrid::xMin() const {
	return xMin_;
}
inline real BackgroundGrid::xMax() const {
	return xMax_;
}
inline real BackgroundGrid::yMin() const {
	return yMin_;
}
inline real BackgroundGrid::yMax() const {
	return yMax_;
}

inline real& BackgroundGrid::operator()( size_t i, size_t j ) {
	return v_[lexIdx( i, j )];
}
inline real BackgroundGrid::operator()( size_t i, size_t j ) const {
	return v_[lexIdx( i, j )];
}
inline real& BackgroundGrid::operator()( size_t i ) {
	assert( i < vSize_ );
	return v_[i];
}
inline real BackgroundGrid::operator()( size_t i ) const {
	assert( i < vSize_ );
	return v_[i];
}

inline size_t BackgroundGrid::lexIdx( const size_t &i, const size_t &j ) const {
	assert( i < Nx_ );
	assert( j < Ny_ );
	return i * Ny_ + j;
}

inline void BackgroundGrid::applyToMeshVertices( TriMesh& mesh ) {
	for ( auto vh : mesh.vertices() ) {
		TriMesh::Point p = mesh.point( vh );
		p[2] = getSize( p[0], p[1] );
		//p[2] = 5;
		mesh.set_point( vh, p );
	}
}

inline size_t BackgroundGrid::coords2LexInd( real x, real y ) {
	//assert(x <= xMax_ && x >= xMin_);
	//assert(y <= yMax_ && y >= yMin_);
	if ( x > xMax_ ) x = xMax_;
	if ( x < xMin_ ) x = xMin_;
	if ( y > yMax_ ) y = yMax_;
	if ( y < yMin_ ) y = yMin_;

	size_t i = (size_t)( ( x - xMin_ ) / hx_ );
	size_t j = (size_t)( ( y - yMin_ ) / hy_ );

	return this->lexIdx( i, j );
}
inline size_t BackgroundGrid::x2col( real x ) const {
	//assert(x <= xMax_ && x >= xMin_);
	if ( x > xMax_ ) x = xMax_;
	if ( x < xMin_ ) x = xMin_;
	return (size_t)( ( x - xMin_ ) / hx_ );
}
inline size_t BackgroundGrid::y2row( real y ) const {
	//assert(y <= yMax_ && y >= yMin_);
	if ( y > yMax_ ) y = yMax_;
	if ( y < yMin_ ) y = yMin_;
	return (size_t)( ( y - yMin_ ) / hy_ );
}
inline size_t BackgroundGrid::lexInd2col( const size_t& j ) const {
	assert( j < vSize_ );
	return j % Ny_;
}
inline size_t BackgroundGrid::lexInd2row( const size_t& i ) const {
	assert( i < vSize_ );
	return i / Ny_;
}

inline real BackgroundGrid::x( const size_t &i ) const {
	assert( i < Nx_ );

	return xMin_ + i * hx_;
}
inline real BackgroundGrid::y( const size_t &j ) const {
	assert( j < Ny_ );

	return yMin_ + j * hy_;
}
