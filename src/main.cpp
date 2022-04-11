#include <experimental/filesystem>
// glog
#include <glog/logging.h>

#include "workingDirectory.h"

#include "RatRace.h"
#include "MeshOptimization.h"

#include "QualityMetrics.h"
#include "Stopwatch/Stopwatch.h"
#include "BackgroundGrid/SizeGrid.h"

namespace fs = std::experimental::filesystem;

//template<int = 0> void discreteMeshOptimization( PolyMesh& mesh, const float grid_scale = 0.5f, int n_iter = 100 );

void printQuality( PolyMesh& mesh, const fs::path& outputPath ) {
	// print quality information
	auto nVal2 = 0;
	for( auto vh : mesh.vertices() ) {
		auto fh = *( mesh.vf_begin( vh ) );
		std::vector<PolyMesh::Point> points;
		points.reserve( 4 );
		for( auto fv : mesh.fv_range( fh ) ) {
			points.push_back( mesh.point( fv ) );
		}
		auto q = QualityMetrics::conditionSurf( points );
		if( mesh.valence( vh ) == 2 && q < 0.01 ) {
			++nVal2;
			LOG( INFO ) << "Low quality vertex: " << vh.idx() << " | " << mesh.point( vh );
			if( !vh.is_boundary() )
				LOG( INFO ) << "   Vertex is not on boundary!!!";
		}
	}
	LOG( INFO ) << "Number of valence 2 vertices: " << nVal2;
	std::vector<int> qVec( 10, 0 );
	real qMin = std::numeric_limits<real>::max(), qMean = 0;
	for( auto fh : mesh.faces() ) {
		std::vector<PolyMesh::Point> points;
		points.reserve( 4 );
		for( auto fv : mesh.fv_range( fh ) ) {
			points.push_back( mesh.point( fv ) );
		}
		auto q = QualityMetrics::conditionSurf( points );

		qMin = std::min( q, qMin );
		qMean += q;

		if( q < 0.1 )
			LOG( INFO ) << "Low quality at " << points[0];

		if( q < 0 ) {
			++qVec[0];
		} else
			++qVec[( q * qVec.size() ) - 0.0001];
	}
	qMean /= mesh.n_faces();
	LOG( INFO ) << "qMin = " << qMin;
	LOG( INFO ) << "qMean = " << qMean;
	for( int i = 0; i < qVec.size(); ++i ) {
		LOG( INFO ) << "Quality " << (double)( i ) / qVec.size() << " - " << (double)( i + 1 ) / qVec.size() << "\t|\t" << qVec[i];
	}
	
	//mesh.request_face_colors();
	//for( auto fh : mesh.faces() ) {
	//	auto q = QualityMetrics::conditionSurf( mesh, fh );
	//	PolyMesh::Point c = { 0,0,0 };
	//	if( q > 0.1 )
	//		c = PolyMesh::Point( 255 ) * q;
	//	else
	//		//c[0] = 200 * -q + 55;
	//		c[0] = 255;
	//	rr.mesh().set_color( fh, { c[0],c[1],c[2] } );
	//}
	//OpenMesh::IO::write_mesh( mesh, ( outputPath / "quality.off" ).string(), OpenMesh::IO::Options::FaceColor, 10 );

	// count valence
	int nVal3 = 0, nVal4 = 0, nVal5 = 0;
	for( auto vh : mesh.vertices() ) {
		if( vh.is_boundary() )
			continue;
		auto val = vh.valence();
		if( val < 4 )
			++nVal3;
		else if( val == 4 )
			++nVal4;
		else
			++nVal5;
	}
	LOG( INFO ) << "Valences: " << nVal3 << " < 4 | " << nVal4 << " == 4 | " << nVal5 << " > 4";
	LOG( INFO ) << "singular vertices [%]: " << 100 * (float)( nVal3 + nVal5 ) / (float)( nVal3 + nVal5 + nVal4 );
}

int main( int argc, char* argv[] ) 
{
	//////////////////////////////////////////////////
	//                                              //
	//                 *** Init ***                 //
	//                                              //
	//----------------------------------------------//

	// --- init logging ---
	google::InitGoogleLogging( argv[0] );
	//FLAGS_timestamp_in_logfile_name = false;
	FLAGS_logtostderr = true;
		
	// --- paths ---
	fs::path basePath = WORKING_DIRECTORY;
	LOG( INFO ) << "basePath = " << fs::canonical( basePath );
	LOG_ASSERT( fs::exists( basePath ) );
	std::string meshName = "transition4to1";
	fs::path inputPath = basePath / "meshes";
	fs::path outputPath = basePath / "output";
	fs::path cachePath = basePath / "cache";

	fs::create_directories( outputPath );

	// uncomment for logging to output path
	//FLAGS_logtostderr = false;
	//FLAGS_alsologtostderr = true;
	//FLAGS_log_dir = outputPath.string();

	fs::path meshFile = inputPath / (meshName + ".off");

	LOG( INFO ) << "Mesh: " << meshFile;


	RatRace rr;
	TriMesh initMesh;
	OpenMesh::IO::read_mesh( initMesh, meshFile.string() );
	OpenMesh::IO::write_mesh( initMesh, ( outputPath / "triangleMesh.off" ).string() );
	const auto nQuads = initMesh.n_faces() / 2;
	LOG( INFO ) << "nQuads = " << nQuads;


	//////////////////////////////////////////////////
	//                                              //
	//             *** Tri to Quad ***              //
	//                                              //
	//----------------------------------------------//

	//Stopwatch sw;
	//sw.start();

	//rr.initBlossomQuad( initMesh );			// Blossom-Quad
	//rr.initTriangleMerging( initMesh );		// Triangle Merging
	rr.initHybrid( initMesh );				// Hybrid Approach

	OpenMesh::IO::write_mesh( rr.mesh(), ( outputPath / "hybridMesh.off" ).string() );
	LOG( INFO ) << "Start RatRace";

	// --- merge remaining triangles ---
	rr.run();

	//sw.stop();
	//LOG( INFO ) << "Runtime: " << sw.runtimeStr<Stopwatch::Milliseconds>();
	LOG( INFO ) << "Finish RatRace";
	OpenMesh::IO::write_mesh( rr.mesh(), ( outputPath / "quadMesh.off" ).string() );

	//////////////////////////////////////////////////
	//                                              //
	//           *** post-processing ***            //
	//                                              //
	//----------------------------------------------//

	SizeGrid sizegrid = loadSizeGrid( initMesh, cachePath, meshFile, 1000, 1000 );
	rr.postProcessing( nQuads, &sizegrid );
	// rr.postProcessing( nQuads, &sizegrid, outputPath ); // use this for more intermediate results

	//////////////////////////////////////////////////
	//                                              //
	//          *** print quality info ***          //
	//                                              //
	//----------------------------------------------//
	
	printQuality( rr.mesh(), outputPath );
	OpenMesh::IO::write_mesh( rr.mesh(), ( outputPath / "finalMesh.off" ).string(), OpenMesh::IO::Options::Default, 10 );

	LOG( INFO ) << "RatRace finished";
}