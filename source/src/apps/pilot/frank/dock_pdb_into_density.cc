#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/import_pose/import_pose.hh>

// minimize pose into density
#include <protocols/electron_density/DockFragmentsIntoDensityMover.hh>
#include <protocols/jd2/JobDistributor.hh>

//#include <core/io/silent/ProteinSilentStruct.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>


#include <string>

#include <utility/file/FileName.hh>

OPT_KEY( Integer, bw )
OPT_KEY( Integer, n_to_search )
OPT_KEY( Integer, n_filtered )
OPT_KEY( Integer, n_output )
OPT_KEY( Integer, n_rotpertrans )
OPT_KEY( Integer, searchsep )
OPT_KEY( Integer, movestep )
OPT_KEY( Integer, ncyc )
OPT_KEY( Real, clust_radius )
OPT_KEY( Real, frag_dens )
OPT_KEY( Boolean, min_bb )
OPT_KEY( Boolean, min )
OPT_KEY( Real, point_radius ) // min distance between points
OPT_KEY( Boolean, convolute_single_residue) // convolute pose on middle ca, results in fine point selection
OPT_KEY( Real, laplacian_offset ) // min distance between points

using namespace ObjexxFCL::format;

// main
int main(int argc, char* argv[]) {
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT( bw, "spharm bandwidth", 32 );
		NEW_OPT( n_to_search, "how many translations to search", 100 );
		NEW_OPT( searchsep, "min distance between search points", 3 );
		NEW_OPT( n_filtered,  "how many solutions to take to refinement", 100 );
		NEW_OPT( n_output, "how many solutions to output", 10 );
		NEW_OPT( n_rotpertrans, "how many solutions to output", 10 );
		NEW_OPT( movestep, "grid spacing over which to search", 1 );
		NEW_OPT( ncyc, "Min cycles", 1 );
		NEW_OPT( clust_radius, "Cluster radius", 3.0 );
		NEW_OPT( frag_dens, "Fragment density", 0.9 );
		NEW_OPT( min_bb, "minimize backbone?", false );
		NEW_OPT( min, "rb min?", true );
		NEW_OPT( point_radius, "minimum translation point to point radius", 5 );
		NEW_OPT( convolute_single_residue, "Convolute only on middle reside", false);
		NEW_OPT( laplacian_offset, "Activates laplacian scoring, and sets laplacian filter offset distance ", 0 );

		devel::init(argc, argv);

		// force some options
		option[ out::nooutput ].value(true);

		protocols::electron_density::DockFragmentsIntoDensityMoverOP dock( new protocols::electron_density::DockFragmentsIntoDensityMover );
		dock->setB( option[ bw ] );
		dock->setTopN( option[ n_to_search ] , option[ n_filtered ] , option[ n_output ] );
		dock->setGridStep( option[ movestep ] );
		dock->setMinDist( option[ searchsep ] );
		dock->setNCyc(option[ ncyc ]());
		dock->setClusterRadius(option[ clust_radius ]());
		dock->setFragDens(option[ frag_dens ]());
		dock->setMinBackbone(option[ min_bb ]());
		dock->setDoRefine(option[ min ]());
		dock->setMaxRotPerTrans( option[ n_rotpertrans ]() );
		dock->setPointRadius(option[ point_radius ]());
		dock->setConvoluteSingleR( option[ convolute_single_residue ]());
		dock->setLaplacianOffset( option[ laplacian_offset ]());


		if ( option[ in::file::native ].user() ) {
			core::pose::PoseOP native_pose( new core::pose::Pose() );
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]().name() , core::import_pose::PDB_file);
			dock->setNative( native_pose );
		}

		if ( option[ out::file::silent ].user() ) {
			std::string silent_fn = option[ out::file::silent ]();
			dock->setOutputSilent( silent_fn );
		}

		protocols::jd2::JobDistributor::get_instance()->go( dock );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
