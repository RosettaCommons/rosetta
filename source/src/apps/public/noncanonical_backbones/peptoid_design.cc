// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// based on code from mini/src/apps/public/scenarios/doug_dock_design_min_mod2_cal_cal.cc
//   and https://svn.rosettacommons.org/source/branches/releases/rosetta-3.1/manual/advanced/example_protocol.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

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
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/RandomOmegaFlipMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>
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
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <sstream>

// DEBUG DEBUG
#include <protocols/moves/PyMolMover.hh>


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


// tracer
static thread_local basic::Tracer TR( "PeptoidDesign" );

// application specific options
namespace peptoid_design {
IntegerOptionKey const pert_num( "peptoid_design::pert_num" );
IntegerOptionKey const design_loop_num( "peptoid_design::design_loop_num" );
BooleanOptionKey const cyclic( "peptoid_design::cyclic" );
IntegerVectorOptionKey const peptoid_design_positions( "peptoid_design::peptoid_design_positions" );
}

class PeptoidDesignMover : public Mover {

public:

	//default ctor
	PeptoidDesignMover(): Mover("PeptoidDesignMover"){}

	//default dtor
	virtual ~PeptoidDesignMover(){}

	//methods
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "PeptoidDesignMover"; }

};

typedef utility::pointer::shared_ptr< PeptoidDesignMover > PeptoidDesignMoverOP;
typedef utility::pointer::shared_ptr< PeptoidDesignMover const > PeptoidDesignMoverCOP;


int
main( int argc, char* argv[] )
{
	try {
		/*********************************************************************************************************************
		Common Setup
		**********************************************************************************************************************/

		// add application specific options to options system
		option.add( peptoid_design::pert_num, "Number of iterations of perturbation loop per design" ).def(10);
		option.add( peptoid_design::design_loop_num, "Number of iterations of pertubation and design" ).def(10);
		option.add( peptoid_design::cyclic, "Is the second chain a cycle" ).def(false);

		utility::vector1< core::Size > empty_vector(0);
		option.add( peptoid_design::peptoid_design_positions, "Positions of peptoid to design" ).def( empty_vector );

		// init command line options
		devel::init(argc, argv);

		//create mover instance
		PeptoidDesignMoverOP OD_mover( new PeptoidDesignMover() );

		setup_filter_stats();

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( OD_mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}//main

void
PeptoidDesignMover::apply(
	core::pose::Pose & pose
)
{
	// create score function
	scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS ) );
	score_fxn->set_weight( atom_pair_constraint, 10 );
	score_fxn->set_weight( angle_constraint, 10 );
	score_fxn->set_weight( dihedral_constraint, 10 );

	// get a fold tree suitable for docking (local helper function)
	setup_pert_foldtree( pose );

	// setup the cyclization mover if the pose is cyclic
	if ( option[ peptoid_design::cyclic ].value() == true ) {
		simple_moves::CyclizationMoverOP init_cyclization( new CyclizationMover( 2, true, false, 0 ) );
		init_cyclization->apply( pose );
	}

	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn, 1.0 ) );

	// DEBUG
	//moves::PyMolMoverOP pmm( new protocols::moves::PyMolMover() );
	//pmm->keep_history( true );

	/*********************************************************
	Peptoid Pertubation Phase Setup
	**********************************************************/

	// jump, rot_mag, trans_mag
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover( 1, 0.1, 0.1 ) );

	// get peptoid start and end positions
	Size pep_start( pose.conformation().chain_begin( 2 ) ); Size pep_end( pose.total_residue() );
	TR << "peptoid_start: " << pep_start << " peptoid_end: " << pep_end << std::endl;

	// create movemap for peptoid
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb_true_range(pep_start, pep_end);
	pert_pep_mm->set_jump( 1, true );

	// setup the cyclization mover
	simple_moves::CyclizationMoverOP pert_cyclization( new CyclizationMover( 2, false, true, 1, score_fxn, pert_pep_mm ) );

	// create random torsion mover
	simple_moves::RandomTorsionMoverOP pert_pep_rand_tor( new simple_moves::RandomTorsionMover( pert_pep_mm, 0.5, 10 ) );

	// create an omega flip mover
	simple_moves::RandomOmegaFlipMoverOP pert_pep_omg_flip( new simple_moves::RandomOmegaFlipMover( pert_pep_mm ) );

	// create a random mover to hold the docking, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_rand_tor, 1 );
	//pert_random->add_mover( pert_pep_omg_flip, 0.01 );

	// create a repeat mover
	moves::RepeatMoverOP pert_repeat( new moves::RepeatMover( pert_random, 100 ) );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_repeat );
	//pert_sequence->add_mover( pmm );
	if ( option[ peptoid_design::cyclic ].value() == true ) {
		pert_sequence->add_mover( pert_cyclization );
	}
	//pert_sequence->add_mover( pmm );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, mc ) );

	/*********************************************************
	Design Setup
	**********************************************************/

	using core::pack::task::operation::TaskOperationCOP;
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

	// maybe a restrict to interface operation

	/*********************************************************
	Minimize Setup
	**********************************************************/

	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );
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

	// create a list of peptoid sidechains (this is inefficient)
	chemical::ResidueTypeCOPs const & rt_caps( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->residue_types_DO_NOT_USE() );
	std::set< std::string > peptoid_name_set;
	for ( Size i(1); i <= rt_caps.size(); ++i ) {
		if ( rt_caps[i]->is_peptoid() ) {
			peptoid_name_set.insert( rt_caps[i]->name3() );
		}
	}

	peptoid_name_set.erase( peptoid_name_set.find( "004" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "006" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "010" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "013" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "111" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "206" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "307" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "313" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "402" ) );
	peptoid_name_set.erase( peptoid_name_set.find( "405" ) );

	peptoid_name_set.erase( peptoid_name_set.find( "623" ) ); // wrong order of chi angle atoms, not sure why we have been able to pack with this before?
	peptoid_name_set.erase( peptoid_name_set.find( "624" ) ); // wrong order of chi angle atoms, not sure why we have been able to pack with this before?

	peptoid_name_set.erase( peptoid_name_set.find( "701" ) ); // not sure what is up but don't want to design this anyway
	peptoid_name_set.erase( peptoid_name_set.find( "702" ) ); // not sure what is up but don't want to design this anyway
	peptoid_name_set.erase( peptoid_name_set.find( "703" ) ); // not sure what is up but don't want to design this anyway
	peptoid_name_set.erase( peptoid_name_set.find( "704" ) ); // not sure what is up but don't want to design this anyway

	TR << "Allowed sidechains: ";
	for ( std::set< std::string >::const_iterator j(peptoid_name_set.begin()); j != peptoid_name_set.end(); ++j ) {
		TR << *j << " ";
	}
	TR << std::endl;


	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/

	TR << "Main loop..." << std::endl;

	ncbb_design_main_loop( Size( option[ peptoid_design::design_loop_num ].value() ),
		Size( option[ peptoid_design::pert_num ].value() ),
		pose,
		pert_trial,
		option[ peptoid_design::peptoid_design_positions ].value(),
		pep_start,
		pep_end,
		desn_ta_min,
		score_fxn,
		mc
	);

	mc->recover_low( pose );

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	curr_job->add_string_real_pair( "ENERGY_FINAL ", (*score_fxn)(pose) );

	TR << "Ending main loop..." << std::endl;
	TR << "Checking pose energy..." << std::endl;
	calculate_statistics( curr_job, pose, score_fxn );
}


