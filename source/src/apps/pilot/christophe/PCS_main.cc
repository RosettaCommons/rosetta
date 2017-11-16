// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file PCS_main.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


// Unit headers
#include <apps/pilot/christophe/PCS_main.hh>

// Package headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftTensor.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergy.hh>
#include <core/util/Tracer.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>


////
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>


#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/JumpClaimer.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>
#include <protocols/topology_broker/RigidChunkClaimer.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
#include <protocols/topology_broker/LoopFragmentClaimer.hh>
#include <protocols/topology_broker/CutBiasClaimer.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <protocols/abinitio/FragmentSampler.hh>
#include <protocols/abinitio/ConstraintFragmentSampler.hh>
#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/abinitio/BrokerMain.hh>


// keep these headers first for compilation with Visual Studio C++
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>

// Unit Headers


// Package Headers
#include <core/kinematics/util.hh>

#include <protocols/abinitio/FragmentMover.hh>
#include <protocols/abinitio/SmoothFragmentMover.hh>
#include <protocols/abinitio/GunnCost.hh>


//#include <protocols/Protocol.hh>
#include <protocols/relax/util.hh>
#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/jumping/SecondaryStructure.hh>
// Project Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <basic/MetricValue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>

#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/chemical/ChemicalManager.hh>


#include <core/conformation/util.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
//#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>


#include <core/sequence/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
#include <protocols/evaluation/TimeEvaluator.hh>
#include <protocols/evaluation/PCA.hh>
#include <protocols/evaluation/PoseMetricEvaluator.hh>
#include <protocols/evaluation/ConstraintEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/evaluation/ChemicalShiftEvaluator.hh>

#include <protocols/loops/SlidingWindowLoopClosure.hh>
#include <protocols/loops/WidthFirstSlidingWindowLoopClosure.hh>
//#include <protocols/loops/ShortLoopClosure.hh>
//#include <protocols/loops/LoopClosure.hh>
#include <protocols/loops/LoopClass.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/util.hh>

#include <protocols/idealize/idealize.hh>

#include <protocols/viewer/viewers.hh>


//#include <protocols/loops/looprelax_protocols.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/RGFilter.hh>
#include <protocols/simple_filters/COFilter.hh>
#include <protocols/simple_filters/SheetFilter.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>

#include <protocols/moves/MoverStatus.hh>

//numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/util.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <protocols/topology_broker/JumpClaimer.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>
#include <protocols/topology_broker/RigidChunkClaimer.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
#include <protocols/topology_broker/LoopFragmentClaimer.hh>
#include <protocols/topology_broker/CutBiasClaimer.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <protocols/abinitio/FragmentSampler.hh>
#include <protocols/abinitio/ConstraintFragmentSampler.hh>
#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/abinitio/BrokerMain.hh>


////

#include <protocols/abinitio/BrokerMain.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


// Numeric headers

// Objexx headers

// C++ headers


static basic::Tracer TR_PCS_main( "apps.pilot.christophe.PCS_main" );

void run() {

	using namespace core;
	using namespace protocols;
	using namespace abinitio;
	using namespace fragment;
	using namespace jumping;
	using namespace evaluation;
	using namespace basic::options;
	//using namespace basic::options::OptionKeys;
	if ( !option[ OptionKeys::in::path::database ].user() ) {
		option[ OptionKeys::in::path::database ].def( "/home/christophe/minirosetta_database_svn_PCS");
	}

	if ( option[ OptionKeys::run::checkpoint ] || option[ OptionKeys::run::checkpoint_interval ].user() ) {
		protocols::checkpoint::checkpoint_with_interval( option[ OptionKeys::run::checkpoint_interval ] );
	}

	protocols::abinitio::Broker_main();

	return;
}

int
main( int argc, char* argv [] )
{
	try {

		using namespace core;
		using namespace scoring;
		using namespace methods;
		using namespace pcs;
		using namespace pose;
		using namespace basic::options;


		pose::Pose pdb;
		utility::vector1<std::string> vec_filename;
		utility::vector1<PCS_tensor> vec_tensor;
		utility::vector1<core::Real> vec_score;
		core::Real pcs_score_total;
		core::Size i;


		protocols::abinitio::register_options_broker();

		devel::init( argc, argv );

		std::string const native( option[ OptionKeys::in::file::native ]() );
		core::import_pose::pose_from_file( pdb, native , core::import_pose::PDB_file);

		PCS_Energy PCS_e = PCS_Energy();

		for ( i = 1; i <= 1; ++i ) {
			pcs_score_total = PCS_e.calculate_pcs_score(pdb, true);
			TR_PCS_main << "Score for iteration " << i << " " << pcs_score_total << std::endl;
		}


		//utility_exit();

		TR_PCS_main << "CALLING BROKER_MAIN" << std::endl;

		//protocols::abinitio::Broker_main();

		run();

		TR_PCS_main << "It is going to crash for unknown reason" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return(1);
}
