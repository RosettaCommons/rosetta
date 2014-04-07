// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Application level code for relax-type protocols
/// @detailed

// keep these headers first for compilation with Visual Studio C++
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jd2/Job.hh>
// Unit Headers


// Package Headers
#include <core/kinematics/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
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

#include <core/io/pdb/pose_io.hh>
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
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/util.hh>

#include <protocols/idealize.hh>


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
#include <protocols/topology_broker/weights/all.hh>
#include <protocols/topology_broker/SequenceClaimer.hh>
#include <protocols/topology_broker/RigidChunkClaimer.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
#include <protocols/topology_broker/LoopFragmentClaimer.hh>
#include <protocols/topology_broker/CutBiasClaimer.hh>
#include <protocols/jumping/SecondaryStructure.hh>
#include <protocols/abinitio/FragmentSampler.hh>
#include <protocols/abinitio/ConstraintFragmentSampler.hh>
#include <protocols/abinitio/AbrelaxMover.hh>


////////////////////////////////////////////////////////////////////////////////////////////////////


using core::Size;
using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;
using namespace jumping;
using namespace evaluation;
using namespace basic::options;
//using namespace basic::options::OptionKeys;



class TestJD2Mover : public moves::Mover
{
public:
	/// @brief
	/// 	empty constructor fills values with the values
	///		read in from the commandline
	TestJD2Mover() :
		Mover()
	{
		Mover::type( "TestJD2Mover" );
	}


	virtual void apply( core::pose::Pose & pose ){

		protocols::jd2::JobDistributor *jd = protocols::jd2::JobDistributor::get_instance();

		std::cout << jd->current_job()->input_tag()
		          << jd->current_job()->nstruct_index()
							<< std::endl;


	};
	virtual void test_move( core::pose::Pose & pose ){apply(pose);}

private:  //  Data
};

typedef utility::pointer::owning_ptr< TestJD2Mover > TestJD2MoverOP;
typedef utility::pointer::owning_ptr< TestJD2Mover const > TestJD2MoverCOP;





void run() {
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/mtyka/minirosetta_database");
	}


	TestJD2MoverOP test = new TestJD2Mover;
	protocols::jd2::JobDistributor::get_instance()->go( test );
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		devel::init( argc, argv );
		run();

	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}

	return 0;
}


