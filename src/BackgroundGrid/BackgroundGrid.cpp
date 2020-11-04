#include "BackgroundGrid.h"

#include <algorithm>
#include <fstream>
#include <tuple>
#include <functional>

inline real triangleArea( real x1, real y1, real x2, real y2, real x3, real y3 ) {
	return std::fabs( ( x1 * ( y2 - y3 ) + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 ) ) * 0.5f );
}

inline real signedTriangleArea( real x1, real y1, real x2, real y2, real x3, real y3 ) {
	return ( x1 * ( y2 - y3 ) + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 ) ) * 0.5f;
}

inline real triangleArea( const TriMesh::Point& p1, const TriMesh::Point& p2, const TriMesh::Point& p3 ) {
	return triangleArea( p1[0], p1[1], p2[0], p2[1], p3[0], p3[1] );
}

inline auto barycentricCoordinates( const std::array<TriMesh::Point, 3>& triangle, const TriMesh::Point& p ) {
	const real A = signedTriangleArea( triangle[0][0], triangle[0][1], triangle[1][0], triangle[1][1], triangle[2][0], triangle[2][1] );
	const real Ainv = 1.f / A;
	const real A0 = signedTriangleArea( p[0], p[1], triangle[1][0], triangle[1][1], triangle[2][0], triangle[2][1] );
	const real A1 = signedTriangleArea( triangle[0][0], triangle[0][1], p[0], p[1], triangle[2][0], triangle[2][1] );
	const real A2 = signedTriangleArea( triangle[0][0], triangle[0][1], triangle[1][0], triangle[1][1], p[0], p[1] );

	return std::make_tuple( A0 * Ainv, A1 * Ainv, A2 * Ainv );
}

inline real barycentricInterpolation( const std::vector<TriMesh::Point>& triangle, const real& a, const real& b, const real& c ) {
	return a * triangle[0][2] + b * triangle[1][2] + c * triangle[2][2];
}

// Check if barycentric coordinates are inside the triangle, i.e. very close to zero or negative.
inline bool isInsideTriangle( const real& a, const real& b, const real& c ) {
	const real eps = 1e-5;
	if( a < eps || b < eps || c < eps ) {
		return false;
	} else {
		return true;
	}
}

std::tuple<real, real, real, real> BackgroundGrid::findAABB( TriMesh& mesh ) {
	const auto n_vertices = mesh.n_vertices();
	auto xMax = -std::numeric_limits<real>::max();
	auto yMax = -std::numeric_limits<real>::max();
	auto xMin = std::numeric_limits<real>::max();
	auto yMin = std::numeric_limits<real>::max();
	for ( unsigned int i = 0; i < n_vertices; ++i ) {
		TriMesh::Point point = mesh.point( mesh.vertex_handle( i ) );

		xMax = std::max( xMax, point[0] );
		yMax = std::max( yMax, point[1] );
		xMin = std::min( xMin, point[0] );
		yMin = std::min( yMin, point[1] );
	}
	// to make sure that all points lie in the interior of the grid, enlarge the grid a bit
	// TODO this should not be necessary. Points are mapped to interior anyways!!!
	const real x_range_old = xMax - xMin;
	const real y_range_old = yMax - yMin;
	xMax += 0.1 * x_range_old;
	xMin -= 0.1 * x_range_old;
	yMax += 0.1 * y_range_old;
	yMin -= 0.1 * y_range_old;

	return std::make_tuple( xMax, xMin, yMax, yMin );
}

BackgroundGrid::BackgroundGrid( const std::experimental::filesystem::path& filename ) {
	if ( !std::experimental::filesystem::exists( filename ) ) {
		LOG( ERROR ) << "File with name '" << filename.string() << "' does not exist!";
	}

	if ( filename.extension() == ".bin" ) {
		readBinary( filename );
	} else {
		LOG(ERROR) << "Extension of file '" << filename.string() << "' is not compatible with BackgroundGrid!";
	}
}

BackgroundGrid::BackgroundGrid( TriMesh & mesh, decltype( valuePerVertex_ ) valuePerVertex, const size_t & Nx, const size_t & Ny, const real scalingfactor ) : valuePerVertex_{ valuePerVertex } {
	initGrid( mesh, Nx, Ny, scalingfactor );

	std::vector<real> vertexAvg( mesh.n_vertices() );
	collectVertexValues( mesh, vertexAvg );

	fillGrid( mesh, vertexAvg );
}

BackgroundGrid::~BackgroundGrid() {
}

void BackgroundGrid::toVTK( const std::experimental::filesystem::path& name, const bool& onlyInterior, const bool& isBinary ) const {
	using std::endl;

	auto swapEnd = [](real& var) {
		char* varArray = reinterpret_cast<char*>( &var );
		for ( long i = 0; i < static_cast<long>( sizeof( var ) / 2 ); i++ )
			std::swap( varArray[sizeof( var ) - 1 - i], varArray[i] );
	};
	
	std::ofstream ofs;
	if ( isBinary )
		ofs.open( name, std::ios::binary | std::ios::out );
	else
		ofs.open( name );

	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << name.stem().string() << endl;
	if(isBinary)
		ofs << "BINARY" << endl;
	else
		ofs << "ASCII" << endl;
	ofs << "DATASET STRUCTURED_POINTS" << endl;
	ofs << "DIMENSIONS " << Nx_ << " " << Ny_ << " " << 1 << endl;
	ofs << "ORIGIN " << xMin_ << " " << yMin_ << " " << 0 << endl;
	ofs << "SPACING " << hx_ << " " << hy_ << " " << 1 << endl;
	
	ofs << "POINT_DATA " << vSize_ << endl;
	ofs << "SCALARS value real 1" << endl;
	ofs << "LOOKUP_TABLE data_table" << endl;
	
	for ( size_t j = 0; j < Ny_; ++j ) {
		for ( size_t i = 0; i < Nx_; ++i ) {
			if ( isBinary ) {
				real val;
				if ( onlyInterior && !isDomain_[lexIdx( i, j )] )
					val = 0;
				else
					val = ( *this )( i, j );
				swapEnd( val );
				ofs.write( (char*)&val, sizeof( val ) );
			} 
			else {
				if ( onlyInterior && !isDomain_[lexIdx( i, j )] )
					ofs << "0 ";
				else
					ofs << ( *this )( i, j ) << "\n";
			}
		}
	}
	if ( isBinary )
		ofs << "\n";
	
	ofs.close();
}

void BackgroundGrid::compare( const BackgroundGrid& sg, const std::experimental::filesystem::path& name ) const {
	assert( Nx_ == sg.Nx_ );
	assert( Ny_ == sg.Ny_ );

	// norm both on average edge length
	real v1_avg = 0;		// this size grid
	real v2_avg = 0;		// the other size grid
	for ( size_t i = 0; i < vSize_; ++i ) {
		v1_avg += v_[i];
		v2_avg += sg.v_[i];
	}
	v1_avg /= (real)vSize_;
	v2_avg /= (real)vSize_;

	real* v_comp = new real[vSize_];

	for ( size_t i = 0; i < vSize_; ++i ) {
		v_comp[i] = 0.5 * std::log2f( v_[i] * v2_avg / ( sg.v_[i] * v1_avg ) );
	}

	std::ofstream ofs( name );

	// Create Header
	ofs << "P3\n#Output from BackgroundGrid\n" << Nx_ << " " << Ny_ << "\n255\n" << std::endl;

	//std::ofstream ofs_vtk("../../vtkOutput.vtk");
	//ofs_vtk << "# vtk DataFile Version 1.0" << std::endl;
	//ofs_vtk << "Density comparison for BSG generation" << std::endl;
	//ofs_vtk << "ASCII" << std::endl;
	//ofs_vtk << "DATASET POLYDATA" << std::endl;
	//ofs_vtk << "POINTS " << Nx_ * Ny_ << " real" << std::endl;

	for ( size_t i = 0; i < Ny_; ++i ) {
		for ( size_t j = 0; j < Nx_; ++j ) {
			if ( !isDomain_[lexIdx( j, Ny_ - 1 - i )] ) {
				ofs << "255 255 255 ";
				//ofs_vtk << this->getX(j) << " " << this->getY(i) << 0 << std::endl;
				continue;
			}

			size_t lexInd = ( *this ).lexIdx( j, Ny_ - 1 - i );

			real value = v_comp[lexInd];

			bool signbit = std::signbit( value );

			value = std::abs( value );

			if ( value > 1 ) value = 1;
			value = 1 - value;

			// value between 0 and 360
			//value = 120 * value;

			if ( signbit ) {
				ofs << (int)( value * 255 ) << " " << (int)( value * 255 ) << " " << 255 << " ";
				//ofs_vtk << this->getX(j) << " " << this->getY(i) << -value << std::endl;
			} else {
				ofs << 255 << " " << (int)( value * 255 ) << " " << (int)( value * 255 ) << " ";
				//ofs_vtk << this->getX(j) << " " << this->getY(i) << value << std::endl;
			}

			//rgb v_rgb = this->h2rgb(value);
			//
			//ofs << (int)(value * 255) << " " << (int)(value * 255) << " " << (int)(value * 255) << " ";
		}
		ofs << std::endl;
	}

	ofs.close();
	//ofs_vtk.close();

	delete[] v_comp;
}


void BackgroundGrid::writeBinary( const std::experimental::filesystem::path& name ) const {
	std::ofstream ofs( name, std::ios::binary | std::ios::out );
	if ( !ofs.is_open() ) {
		std::cerr << "Cannot write to file with name '" << name << "'" << std::endl;
		std::cerr << "Line: " << __LINE__ << "\nFile: " << __FILE__ << std::endl;
	}

	auto write = [ &ofs ]( const auto& x ) {
		ofs.write( (char*)&x, sizeof( x ) );
	};

	write( Nx_ );
	write( Ny_ );
	write( xMin_ );
	write( xMax_ );
	write( yMin_ );
	write( yMax_ );
	ofs.write( (char*)&v_[0], v_.size() * sizeof( v_[0] ) );

	// store isDomain_ in a char-vector
	std::vector<char>domCharBuf( vSize_ / 8 + 1, 0 );
	for ( int i = 0; i < domCharBuf.size(); ++i ) {
		for ( int j = 0; j < 8; ++j ) {
			int idx = i * 8 + j;
			if ( idx >= vSize_ - 1 ) break;
			domCharBuf[i] |= ( isDomain_[idx] << j );
		}
	}
	ofs.write( (char*)&domCharBuf[0], domCharBuf.size() * sizeof( domCharBuf[0] ) );

	ofs.close();
}

void BackgroundGrid::readBinary( const std::experimental::filesystem::path & name ) {
	if ( !std::experimental::filesystem::exists( name ) ) {
		std::cerr << "File with name '" << name.string() << "' does not exist!" << std::endl;
		std::cerr << "Line: " << __LINE__ << "\nFile: " << __FILE__ << std::endl;
	}

	std::ifstream ifs( name, std::ios::binary | std::ios::in );
	if ( !ifs.is_open() ) {
		std::cerr << "Cannot read file with name '" << name << "'" << std::endl;
		std::cerr << "Line: " << __LINE__ << "\nFile: " << __FILE__ << std::endl;
	}

	auto read = [ &ifs ]( const auto& x ) {
		ifs.read( (char*)&x, sizeof( x ) );
	};

	read( Nx_ );
	read( Ny_ );
	read( xMin_ );
	read( xMax_ );
	read( yMin_ );
	read( yMax_ );

	vSize_ = Nx_ * Ny_;
	v_.resize( vSize_ );

	hx_ = ( xMax_ - xMin_ ) / ( Nx_ - 1 );
	hy_ = ( yMax_ - yMin_ ) / ( Ny_ - 1 );

	ifs.read( (char*)&v_[0], v_.size() * sizeof( v_[0] ) );

	std::vector<char>domCharBuf( vSize_ / 8 + 1, 0 );
	ifs.read( (char*)&domCharBuf[0], domCharBuf.size() * sizeof( domCharBuf[0] ) );

	std::vector<bool>domBoolBuf( vSize_, false );
	for ( int i = 0; i < domCharBuf.size(); ++i ) {
		for ( int j = 0; j < 8; ++j ) {
			int idx = i * 8 + j;
			if ( idx >= vSize_ - 1 ) break;

			domBoolBuf[idx] = domCharBuf[i] & ( 1 << j );
		}
	}

	isDomain_ = domBoolBuf;


	ifs.close();
}

void BackgroundGrid::initGrid( TriMesh & mesh, const size_t & Nx, const size_t & Ny, const real scalingfactor ) {
	Nx_ = Nx;
	Ny_ = Ny;
	// max & min values
	std::tie( xMax_, xMin_, yMax_, yMin_ ) = findAABB( mesh );

	if ( Nx_ == 0 || Ny_ == 0 ) {
		LOG( WARNING ) << "\nSize Grid dimensions not stated or zero. Use adaptive dimension settings";
		setGridSize( mesh, findMinEdgeLength( mesh ) );
	}

	vSize_ = Nx_ * Ny_;
	v_.resize( vSize_ );
	for ( auto& e : v_ ) e = -FLT_MAX;
	isDomain_ = std::vector<bool>( Nx_ * Ny_ );

	// Grid spacing
	hx_ = ( xMax_ - xMin_ ) / ( Nx_ - 1 );
	hy_ = ( yMax_ - yMin_ ) / ( Ny_ - 1 );
}

void BackgroundGrid::fillGrid( TriMesh& mesh, const std::vector<real>& vertex_avg ) {
	//// iterate over all faces
	const auto nFaces = mesh.n_faces();
	#pragma omp parallel for
	for ( int f{ 0 }; f < nFaces; ++f ) {
		TriMesh::FaceHandle fh = mesh.face_handle( f );
		std::vector<TriMesh::Point> triangle;	// points of the triangle

		for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); ++fv_it ) {
			TriMesh::Point p{ mesh.point( *fv_it )[0], mesh.point( *fv_it )[1], vertex_avg[fv_it->idx()] };
			triangle.push_back( p );
		}


		// get the rectangle of grid points that surrounds the triangle
		real x_t_max = std::max( std::max( triangle[0][0], triangle[1][0] ), triangle[2][0] );
		real x_t_min = std::min( std::min( triangle[0][0], triangle[1][0] ), triangle[2][0] );
		real y_t_max = std::max( std::max( triangle[0][1], triangle[1][1] ), triangle[2][1] );
		real y_t_min = std::min( std::min( triangle[0][1], triangle[1][1] ), triangle[2][1] );
		size_t i_min = x2col( x_t_min );
		size_t i_max = x2col( x_t_max ) + 1;
		size_t j_min = y2row( y_t_min );
		size_t j_max = y2row( y_t_max ) + 1;

		// check for all inner points of rectangle if they are inside the triangle. If yes, calculate the depth
		for ( size_t i = i_min + 1; i < i_max; ++i ) {
			const real xv = x( i );

			for ( size_t j = j_min + 1; j < j_max; ++j ) {
				const real yv = y( j );
				const auto[a, b, c] = barycentricCoordinates( { triangle[0], triangle[1], triangle[2] }, { xv,yv,0 } );

				if ( isInsideTriangle( a, b, c ) ) {
					// calculate average edge length
					( *this )( i, j ) = barycentricInterpolation( triangle, a, b, c );
					isDomain_[lexIdx( i, j )] = true;
				}
			}
		}

	}

	setOuterVals();
}

real BackgroundGrid::findMinEdgeLength( TriMesh& mesh ) {
	real minlength_squared = std::numeric_limits<real>::max();

	for ( TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it ) {
		const TriMesh::HalfedgeHandle  heh = mesh.halfedge_handle( *e_it, 0 );
		const TriMesh::VertexHandle vh1 = mesh.to_vertex_handle( heh );
		const TriMesh::VertexHandle vh2 = mesh.from_vertex_handle( heh );

		const real x_dist = mesh.point( vh2 )[0] - mesh.point( vh1 )[0];
		const real y_dist = mesh.point( vh2 )[1] - mesh.point( vh1 )[1];
		const real d_squared = x_dist * x_dist + y_dist * y_dist;

		minlength_squared = std::min( minlength_squared, d_squared );
	}

	return sqrt( minlength_squared );
}

void BackgroundGrid::collectVertexValues( TriMesh mesh, std::vector<real>& vertex_avg ) {
	// go over all vertices
	const auto nVertices = mesh.n_vertices();

	#pragma omp parallel for
	for ( auto i = 0; i < nVertices; ++i ) {
		TriMesh::VertexHandle vh = mesh.vertex_handle( i );
		vertex_avg[i] = valuePerVertex_(mesh, vh);
	}
}


inline void gaussSeidelStep( BackgroundGrid& sg, int i, int j, real* v, real* vNew ) {
	int n = 0;
	real val = 0;
	if ( i > 0 && v[sg.lexIdx( i - 1, j )] != -FLT_MAX ) {
		val += v[sg.lexIdx( i - 1, j )];
		++n;
	}
	if ( i < sg.Nx() - 1 && v[sg.lexIdx( i + 1, j )] != -FLT_MAX ) {
		val += v[sg.lexIdx( i + 1, j )];
		++n;
	}
	if ( j > 0 && v[sg.lexIdx( i, j - 1 )] != -FLT_MAX ) {
		val += v[sg.lexIdx( i, j - 1 )];
		++n;
	}
	if ( j < sg.Ny() - 1 && v[sg.lexIdx( i, j + 1 )] != -FLT_MAX ) {
		val += v[sg.lexIdx( i, j + 1 )];
		++n;
	}

	if ( n != 0 )
		vNew[sg.lexIdx( i, j )] = val / n;
}

inline void gaussSeidelStepInterior( BackgroundGrid& sg, int i, int j, real* v, real* vNew ) {
	vNew[sg.lexIdx( i, j )] = 0.25f * (
		v[sg.lexIdx( i - 1, j )] +
		v[sg.lexIdx( i + 1, j )] +
		v[sg.lexIdx( i, j - 1 )] +
		v[sg.lexIdx( i, j + 1 )] );
}

void BackgroundGrid::setOuterVals() {

	auto v2 = v_;
	auto vDiff = decltype( v_ )( vSize_, 0 );
	real* pv1 = v_.data();
	real* pv2 = v2.data();

	const auto& v1Size = vSize_;

	std::vector<std::array<int, 2>> warmUpRed;
	std::vector<std::array<int, 2>> warmUpBlack;

	for ( int i = 0; i < Nx_; ++i ) {
		for ( int j = 0; j < Ny_; ++j ) {
			if ( isDomain_[lexIdx( i, j )] ) continue;
			if ( ( i + j ) % 2 == 0 )
				warmUpRed.push_back( { i,j } );
			else
				warmUpBlack.push_back( { i,j } );
		}
	}

	const auto warmUpRedSize = warmUpRed.size();
	const auto warmUpBlackSize = warmUpBlack.size();

	// warm-up phase (spread reasonable values fast)
	for ( long k = 0; k < 1e6; ++k ) {

		#pragma omp parallel for
		for ( int idx = 0; idx < warmUpRedSize; ++idx ) {
			const int i = warmUpRed[idx][0];
			const int j = warmUpRed[idx][1];
			gaussSeidelStep( *this, i, j, pv1, pv2 );
		}
		#pragma omp parallel for
		for ( int idx = 0; idx < warmUpBlackSize; ++idx ) {
			const int i = warmUpBlack[idx][0];
			const int j = warmUpBlack[idx][1];
			gaussSeidelStep( *this, i, j, pv2, pv2 );
		}

		std::swap( pv1, pv2 );
		auto minIt = std::min_element( v_.begin(), v_.end() );

		if ( *minIt != -FLT_MAX ) {
			// every node has a reasonable value
			//std::cout << "nIterations = " << k << "   min value = " << *minIt << std::endl;
			break;
		}
	}

	std::vector<std::array<int, 2>> interiorRed;
	std::vector<std::array<int, 2>> interiorBlack;

	for ( int i = 1; i < Nx_ - 1; ++i ) {
		for ( int j = 1; j < Ny_ - 1; ++j ) {
			if ( isDomain_[lexIdx( i, j )] ) continue;
			if ( ( i + j ) % 2 == 0 )
				interiorRed.push_back( { i,j } );
			else
				interiorBlack.push_back( { i,j } );
		}
	}

	const auto interiorRedSize = interiorRed.size();
	const auto interiorBlackSize = interiorBlack.size();

	// main loop (after warm-up)
	for ( long k = 0; k < 1e6; ++k ) {

		#pragma omp parallel for
		for ( int idx = 0; idx < interiorRedSize; ++idx ) {
			const auto i = interiorRed[idx][0];
			const auto j = interiorRed[idx][1];
			gaussSeidelStepInterior( *this, i, j, pv1, pv2 );
		}
		// edges red
		for ( int i = 2; i < Nx_ - 1; i += 2 ) {
			const auto j = 0;
			pv2[lexIdx( i, j )] = (
				pv1[lexIdx( i - 1, j )] +
				pv1[lexIdx( i + 1, j )] +
				pv1[lexIdx( i, j + 1 )] ) / 3;
		}
		for ( int j = 2; j < Ny_ - 1; j += 2 ) {
			const auto i = 0;
			pv2[lexIdx( i, j )] = (
				pv1[lexIdx( i + 1, j )] +
				pv1[lexIdx( i, j - 1 )] +
				pv1[lexIdx( i, j + 1 )] ) / 3;
		}
		for ( int i = Ny_ % 2 == 0 ? 1 : 2; i < Nx_ - 1; i += 2 ) {
			const auto j = Ny_ - 1;
			pv2[lexIdx( i, j )] = (
				pv1[lexIdx( i - 1, j )] +
				pv1[lexIdx( i + 1, j )] +
				pv1[lexIdx( i, j - 1 )] ) / 3;
		}
		for ( int j = Nx_ % 2 == 0 ? 1 : 2; j < Ny_ - 1; j += 2 ) {
			const auto i = Nx_ - 1;
			pv2[lexIdx( i, j )] = (
				pv1[lexIdx( i - 1, j )] +
				pv1[lexIdx( i, j - 1 )] +
				pv1[lexIdx( i, j + 1 )] ) / 3;
		}
		// corners red
		pv2[lexIdx( 0, 0 )] = ( pv1[lexIdx( 0 + 1, 0 )] + pv1[lexIdx( 0, 0 + 1 )] ) / 2;
		if ( Ny_ % 2 == 1 )
			pv2[lexIdx( Nx_ - 1, 0 )] = ( pv1[lexIdx( Nx_ - 1 - 1, 0 )] + pv1[lexIdx( Nx_ - 1, 0 + 1 )] ) / 2;
		if ( Nx_ % 2 == 1 )
			pv2[lexIdx( 0, Ny_ - 1 )] = ( pv1[lexIdx( 0 + 1, Ny_ - 1 )] + pv1[lexIdx( 0, Ny_ - 1 - 1 )] ) / 2;
		if ( ( Ny_ + Nx_ ) % 2 == 0 )
			pv2[lexIdx( Nx_ - 1, Ny_ - 1 )] = ( pv1[lexIdx( Nx_ - 1 - 1, Ny_ - 1 )] + pv1[lexIdx( Nx_ - 1, Ny_ - 1 - 1 )] ) / 2;

		#pragma omp parallel for
		for ( int idx = 0; idx < interiorBlackSize; ++idx ) {
			const int i = interiorBlack[idx][0];
			const int j = interiorBlack[idx][1];
			gaussSeidelStepInterior( *this, i, j, pv2, pv2 );
		}
		// edges black
		for ( int i = 1; i < Nx_ - 1; i += 2 ) {
			const int j = 0;
			pv2[lexIdx( i, j )] = (
				pv2[lexIdx( i - 1, j )] +
				pv2[lexIdx( i + 1, j )] +
				pv2[lexIdx( i, j + 1 )] ) / 3;
		}
		for ( int j = 1; j < Ny_ - 1; j += 2 ) {
			const int i = 0;
			pv2[lexIdx( i, j )] = (
				pv2[lexIdx( i + 1, j )] +
				pv2[lexIdx( i, j - 1 )] +
				pv2[lexIdx( i, j + 1 )] ) / 3;
		}
		for ( int i = Ny_ % 2 == 1 ? 1 : 2; i < Nx_ - 1; i += 2 ) {
			const int j = Ny_ - 1;
			pv2[lexIdx( i, j )] = (
				pv2[lexIdx( i - 1, j )] +
				pv2[lexIdx( i + 1, j )] +
				pv2[lexIdx( i, j - 1 )] ) / 3;
		}
		for ( int j = Nx_ % 2 == 1 ? 1 : 2; j < Ny_ - 1; j += 2 ) {
			const int i = Nx_ - 1;
			pv2[lexIdx( i, j )] = (
				pv2[lexIdx( i - 1, j )] +
				pv2[lexIdx( i, j - 1 )] +
				pv2[lexIdx( i, j + 1 )] ) / 3;
		}
		// corners black
		if ( Ny_ % 2 == 0 )
			pv2[lexIdx( Nx_ - 1, 0 )] = ( pv2[lexIdx( Nx_ - 1 - 1, 0 )] + pv2[lexIdx( Nx_ - 1, 0 + 1 )] ) / 2;
		if ( Nx_ % 2 == 0 )
			pv2[lexIdx( 0, Ny_ - 1 )] = ( pv2[lexIdx( 0 + 1, Ny_ - 1 )] + pv2[lexIdx( 0, Ny_ - 1 - 1 )] ) / 2;
		if ( ( Ny_ + Nx_ ) % 2 == 1 )
			pv2[lexIdx( Nx_ - 1, Ny_ - 1 )] = ( pv2[lexIdx( Nx_ - 1 - 1, Ny_ - 1 )] + pv2[lexIdx( Nx_ - 1, Ny_ - 1 - 1 )] ) / 2;

		#pragma omp parallel for
		for ( int i = 0; i < v1Size; ++i ) {
			vDiff[i] = std::fabs( pv2[i] - pv1[i] );
		}
		auto maxIt = std::max_element( vDiff.begin(), vDiff.end() );
		std::swap( pv1, pv2 );

		if ( *maxIt <= 1e-5 * pv1[std::distance( vDiff.begin(), maxIt )] ) {
			//std::cout << *maxIt << " at " << std::distance( vDiff.begin(), maxIt ) << " value1 = " << v1[std::distance( vDiff.begin(), maxIt )] << " value2 = " << v2[std::distance( vDiff.begin(), maxIt )] << "\n";
			//std::cout << "nIterations = " << k << std::endl;
			break;
		}
	}

}

void BackgroundGrid::setGridSize( TriMesh& mesh, real minlength ) {
	Nx_ = (size_t)( ( xMax_ - xMin_ ) / minlength );
	Ny_ = (size_t)( ( yMax_ - yMin_ ) / minlength );
	LOG( WARNING ) << "SizeGrid dimensions:    Nx: " << Nx_ << "  |  Ny: " << Ny_;
}

void BackgroundGrid::setGridSpecs( TriMesh& mesh ) {
	vSize_ = Nx_ * Ny_;
	v_.resize( vSize_ );
	for ( auto& e : v_ ) e = -1;
	isDomain_ = std::vector<bool>( Nx_ * Ny_ );

	// Grid spacing
	hx_ = ( xMax_ - xMin_ ) / ( Nx_ - 1 );
	hy_ = ( yMax_ - yMin_ ) / ( Ny_ - 1 );

	std::vector<real> vertex_avg( mesh.n_vertices() );
	collectVertexValues( mesh, vertex_avg );

	fillGrid( mesh, vertex_avg );

	setOuterVals();
}
