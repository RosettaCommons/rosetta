#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>  //make_pose_from_sequence
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/util/SwitchResidueTypeSet.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

// minimize pose into density
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/electron_density/DockIntoDensityMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>


#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>

//#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <iostream>
#include <string>
#include <list>
#include <algorithm>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

OPT_KEY( Integer, bw )
OPT_KEY( Integer, n_to_search )
OPT_KEY( Integer, n_filtered )
OPT_KEY( Integer, n_output )
OPT_KEY( Integer, movestep )
OPT_KEY( Integer, ncyc )
OPT_KEY( Real, clust_radius )
OPT_KEY( Real, frag_dens )
OPT_KEY( Boolean, min_bb )
OPT_KEY( Boolean, min )

using namespace ObjexxFCL::format;

// main
int main(int argc, char* argv[]) {
try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( bw, "spharm bandwidth", 16 );
	NEW_OPT( n_to_search, "how many translations to search", 1000 );
	NEW_OPT( n_filtered,  "how many solutions to take to refinement", 100 );
	NEW_OPT( n_output, "how many solutions to output", 10 );
	NEW_OPT( movestep, "grid spacing over which to search", 2 );
	NEW_OPT( ncyc, "Min cycles", 1 );
	NEW_OPT( clust_radius, "Cluster radius", 3.0 );
	NEW_OPT( frag_dens, "Fragment density", 0.9 );
	NEW_OPT( min_bb, "minimize backbone?", false );
	NEW_OPT( min, "rb min?", true );

	devel::init(argc, argv);

	// force some options
	option[ out::nooutput ].value(true);

	protocols::electron_density::DockIntoDensityMoverOP dock( new protocols::electron_density::DockIntoDensityMover );
	dock->setB( option[ bw ] );
	dock->setTopN( option[ n_to_search ] , option[ n_filtered ] , option[ n_output ] );
	dock->setGridStep( option[ movestep ] );
	dock->setNCyc(option[ ncyc ]());
	dock->setClusterRadius(option[ clust_radius ]());
	dock->setFragDens(option[ frag_dens ]());
	dock->setMinBackbone(option[ min_bb ]());
	dock->setDoRefine(option[ min ]());

	if( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]().name() );
		dock->setNative( native_pose );
	}

	if (option[ out::file::silent ].user()) {
		std::string silent_fn = option[ out::file::silent ]();
		dock->setOutputSilent( silent_fn );
	}

	protocols::jd2::JobDistributor::get_instance()->go( dock );

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
  return 0;
}
