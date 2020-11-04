#pragma once

#include "BackgroundGrid.h"
#include <typeinfo>

class SizeGrid : public BackgroundGrid
{
public:
	SizeGrid() : BackgroundGrid() {}
	SizeGrid( TriMesh& mesh, const size_t& Nx = 0, const size_t& Ny = 0, const real scalingfactor = 1.0 )
		: BackgroundGrid{ mesh,
			[]( TriMesh& mesh, const TriMesh::VertexHandle vh ) {
				real avg{ 0.f };
				auto nEdges{ 0 };

				for ( auto voh_it = mesh.voh_iter( vh ); voh_it.is_valid(); ++voh_it ) {
					++nEdges;
					auto p1 = mesh.point( mesh.from_vertex_handle( *voh_it ) );
					auto p2 = mesh.point( mesh.to_vertex_handle( *voh_it ) );
					p1[2] = 0;
					p2[2] = 0;
					avg += ( p1 - p2 ).length();
				}

				return avg / (real)nEdges;
			},
			Nx, Ny, scalingfactor } { }
	SizeGrid( const std::experimental::filesystem::path& filename ) : BackgroundGrid{ filename } {}
};

// load size grid if it already exists in cache, otherwise generate and store it.
inline SizeGrid loadSizeGrid( TriMesh& mesh, const std::experimental::filesystem::path& cacheFolder, const std::experimental::filesystem::path& meshFile, const size_t& sizeGridSizeX, const size_t& sizeGridSizeY ) {
	namespace fs = std::experimental::filesystem;

	std::string typeName = typeid( real ).name();

	fs::path sizeGridFile = cacheFolder / ( meshFile.stem().string() + "_SizeGrid_" + std::to_string( sizeGridSizeX ) + "_" + std::to_string( sizeGridSizeY ) + "_" + typeName + ".bin" );

	if ( !fs::exists( sizeGridFile.parent_path() ) ) {
		LOG( INFO ) << "Directory '" << sizeGridFile.parent_path() << "' does not exist yet and is therefore now created";
		if ( !fs::create_directories( sizeGridFile.parent_path() ) ) {
			LOG( ERROR ) << "Directory '" << sizeGridFile.parent_path() << "' cannot be created!";
		}
	}
	if ( fs::exists( sizeGridFile ) ) {
		LOG( INFO ) << "Read SizeGrid from cache file '" << fs::canonical( sizeGridFile ) << "'";
		return SizeGrid( sizeGridFile.string() );
	} else {
		LOG( INFO ) << "Generate SizeGrid";
		SizeGrid sg( mesh, sizeGridSizeX, sizeGridSizeY );
		LOG( INFO ) << "Store in cache file '" << sizeGridFile << "'";
		sg.writeBinary( sizeGridFile.string() );
		return sg;
	}
}
