// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//kdrew: based on code from mini/src/apps/public/scenarios/doug_dock_design_min_mod2_cal_cal.cc
//   and https://svn.rosettacommons.org/source/branches/releases/rosetta-3.1/manual/advanced/example_protocol.cc

//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/oop/OopRandomPuckMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/ncbb/util.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::ncbb;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::oop;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


// tracer - used to replace cout
static thread_local basic::Tracer TR( "OopDesign" );

// application specific options
namespace oop_design {
// pert options

IntegerOptionKey const pert_num( "oop_design::pert_num" );
IntegerOptionKey const design_loop_num( "oop_design::design_loop_num" );

IntegerVectorOptionKey const oop_design_positions( "oop_design::oop_design_positions" );

}

class OopDesignMover : public Mover {

public:

	//default ctor
	OopDesignMover(): Mover("OopDesignMover"){}

	//default dtor
	virtual ~OopDesignMover(){}

	//methods
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "OopDesignMover"; }

};

typedef utility::pointer::shared_ptr< OopDesignMover > OopDesignMoverOP;
typedef utility::pointer::shared_ptr< OopDesignMover const > OopDesignMoverCOP;


int
main( int argc, char* argv[] )
{
	try {


		/*********************************************************************************************************************
		Common Setup
		**********************************************************************************************************************/

		// add application specific options to options system
		// There are far more options here than you will realistically need for a program of this complexity - but this gives you an idea of how to fine-grain option-control everything

		option.add( oop_design::pert_num, "Number of iterations of perturbation loop per design" ).def(10);
		option.add( oop_design::design_loop_num, "Number of iterations of pertubation and design" ).def(10);

		utility::vector1< core::Size > empty_vector(0);
		option.add( oop_design::oop_design_positions, "Positions of oop to design" ).def( empty_vector );

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		//create mover instance
		OopDesignMoverOP OD_mover( new OopDesignMover() );

		setup_filter_stats();

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( OD_mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}//main

void
OopDesignMover::apply(
	core::pose::Pose & pose
)
{
	// create score function
	scoring::ScoreFunctionOP score_fxn( get_score_function() );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	// get a fold tree suitable for docking (local helper function)
	setup_pert_foldtree( pose );

	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn, 1.0 ) );

	// jump, rot_mag, trans_mag
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, 1.0, 0.5 ) );

	/*********************************************************
	Oop Setup
	**********************************************************/

	// get oop start and end positions
	Size pep_start( pose.conformation().chain_begin( 2 ) ); Size pep_end( pose.total_residue() );
	TR << "oop_start: " << pep_start << " oop_end: " << pep_end << std::endl;

	// create movemap for oop
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb_true_range(pep_start, pep_end);

	//kdrew: automatically find oop positions
	utility::vector1< core::Size > oop_seq_positions = core::pose::ncbb::initialize_oops(pose);

	for ( Size i = 1; i <= oop_seq_positions.size(); i++  ) {
		pert_pep_mm->set_bb( oop_seq_positions[i], false );

		if ( score_fxn->has_zero_weight( core::scoring::atom_pair_constraint ) ) {
			score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		}

	}

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.8, 1 ) );

	// create random mover
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, 20 ) ); //repeat 20 times

	/*********************************************************
	OopPuck Setup
	**********************************************************/

	oop::OopRandomSmallMoverOP opm_small( new oop::OopRandomSmallMover( oop_seq_positions, 2.0 ) );
	oop::OopRandomPuckMoverOP opm_puck( new oop::OopRandomPuckMover( oop_seq_positions ) );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a random mover to hold the docking, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_repeat, 0.5 );
	pert_random->add_mover( opm_small, 0.5 );
	pert_random->add_mover( opm_puck, 0.1 );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_random );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, mc ) );

	/*********************************************************
	Design Setup
	**********************************************************/

	using core::pack::task::operation::TaskOperationCOP;
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

	/*********************************************************
	Minimize Setup
	**********************************************************/

	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );
	//kdrew: set backbone of target false and backbone of oop true, decide whether to do this or not
	desn_mm->set_bb( false );
	desn_mm->set_bb_true_range( pep_start, pep_end );
	desn_mm->set_chi( true );
	desn_mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min( new TaskAwareMinMover( desn_min, desn_tf ) );

	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/

	TR << "Main loop..." << std::endl;

	ncbb_design_main_loop( Size( option[ oop_design::design_loop_num ].value() ),
		Size( option[ oop_design::pert_num ].value() ),
		pose,
		pert_trial,
		option[ oop_design::oop_design_positions ].value(),
		pep_start,
		pep_end,
		desn_ta_min,
		score_fxn,
		mc
	);

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	curr_job->add_string_real_pair( "ENERGY_FINAL ", (*score_fxn)(pose) );

	TR << "Ending main loop..." << std::endl;
	TR << "Checking pose energy..." << std::endl;
	calculate_statistics( curr_job, pose, score_fxn );

}
