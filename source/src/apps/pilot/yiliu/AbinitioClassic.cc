// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// keep these headers first for compilation with Visual Studio C++
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>

// Unit Headers
#include <protocols/abinitio/AbrelaxApplication.hh>


// Package Headers
#include <core/kinematics/util.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/JumpingFoldConstraints.hh>
#include <protocols/abinitio/Templates.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/StrandConstraints.hh>
#include <protocols/abinitio/FragmentMover.hh>

#include <protocols/Protocol.hh>
#include <protocols/relax_protocols.hh>

#include <protocols/jumping/SheetBuilder.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/ResiduePairJumpSetup.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <protocols/jumping/StrandPairing.hh>
#include <protocols/jumping/util.hh>

// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/MetricValue.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/chemical/ChemicalManager.hh>


#include <core/conformation/util.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

#include <core/sequence/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
#include <protocols/evaluation/TimeEvaluator.hh>
#include <protocols/evaluation/PCA.hh>
#include <protocols/evaluation/PoseMetricEvaluator.hh>
#include <protocols/evaluation/ConstraintEvaluator.hh>

#include <protocols/loops/SlidingWindowLoopClosure.hh>
#include <protocols/loops/ShortLoopClosure.hh>
#include <protocols/loops/LoopClosure.hh>
#include <protocols/loops/LoopClass.hh>
#include <protocols/loops/util.hh>

#include <protocols/loops/looprelax_protocols.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/RGFilter.hh>
#include <protocols/simple_filters/COFilter.hh>
#include <protocols/simple_filters/SheetFilter.hh>

#include <protocols/abinitio/KinematicTaskControl.hh>

//#include <protocols/simple_moves/MinMover.hh>

//numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/util.hh>
#include <utility/excn/Exceptions.hh>
// C++ headers
#include <cstdlib>
#include <string>
#include <vector>



//==#ifdef BOINC
//==#include <utility/boinc/boinc_util.hh>
//==#include <protocols/boinc/boinc.hh>
//==#include "boinc_zip.h"
//==#endif // BOINC

/// Must have this after BOINC stuff to avoid windows build error
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/rbsegment_relaxRelax_main.hh>
#include <protocols/loop_build/LoopBuild.hh>
#include <protocols/relax_protocols.hh>
#include <protocols/ligand_docking/ligand_dock_impl.hh>>
#include <protocols/abinitio/vs_test.hh>
#include <protocols/protein_interface_design/conditional_jobdist_mover.hh>
#include <protocols/protein_interface_design/setup_dockdesign_mover.hh>

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/MetricValue.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

int

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



main( int argc, char * argv [] )
{
	try {

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using std::string;
	using utility::vector1;

	//YL, move the options register functions out of the boinc section
	// has to be called before devel::init. Which is really stupid.
	protocols::abinitio::ClassicAbinitio::register_options();
	// options, random initialization
	devel::init( argc, argv );
	if ( option[ run::checkpoint ] || option[ run::checkpoint_interval ].user() ) {
		protocols::checkpoint::checkpoint_with_interval( option[ run::checkpoint_interval ] );
  }
	//declare fragments file
	std::string frag_large_file, frag_small_file;

	if ( option[ in::file::fragA ].user() ) {
		frag_large_file  = option[ in::file::fragA ]();
	} else {
		frag_large_file  = option[ in::file::frag9 ]();
	}

	if ( option[ in::file::fragB ].user() ) {
		frag_small_file  = option[ in::file::fragB ]();
	} else {
		frag_small_file  = option[ in::file::frag3 ]();
	}

  using namespace core::fragment;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	// declare 9mer fragments
	core::fragment::FragSetOP fragset_large;
	fragset_large = FragmentIO().read(frag_large_file);

	// declare 3mer fragments
	core::fragment::FragSetOP fragset_small;
	fragset_small = FragmentIO().read( frag_small_file );

	kinematics::MoveMapOP movemap = new kinematics::MoveMap;

	//yl, prepare for classic abinitio protocols
	protocols::abinitio::ClassicAbinitio abinitio(
		fragset_small,
		fragset_large,
		movemap
	);

	//core::pose::PoseOP native_pose;
	//if ( option[ in::file::native ].user() ) {
	//	native_pose = new pose::Pose;
	//	core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
	//	pose::set_ss_from_phipsi( *native_pose );
	//}

	abinitio.init(*native_pose);
	int const nstruct = std::max( 1, option [ out::nstruct ]() );
	while(nstruct){
		abinitio.apply(*native_pose);
	}
	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
