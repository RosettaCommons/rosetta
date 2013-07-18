// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//kdrew: based on code from mini/src/apps/public/scenarios/doug_dock_design_min_mod2_cal_cal.cc
//			and https://svn.rosettacommons.org/source/branches/releases/rosetta-3.1/manual/advanced/example_protocol.cc

//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.


// Project Headers
#include <core/pose/Pose.hh>
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

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
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

// Here is a sketch of the basic flow of the program...
//
// Pertubation Phase
//   +-Monte Carlo Mover---------------------------------------+
//   | +-Random Mover-( 1 / 2 / 1 / 1 )--------------------+ | |
//   | | +-Docking Mover-----------------------------------+ | |
//   | | | small rigid body movements between the peptide	 | | |
//   | | | and protein for conformational diversity        | | |
//   | | +-------------------------------------------------+ | |
//   | | +-Peptide Modeling--------------------------------+ | |
//   | | | move peptide with small/shear moves to generate | | |
//   | | | conformational diversity                        | | |
//   | | +-------------------------------------------------+ | |
//   | +-----------------------------------------------------+ |
//   | +-Rotamer Trials Mover--------------------------------+ |
//   | | quick sidechain packing to find optimal rotamers    | |
//   | | before the next cycle                               | |
//   | +-----------------------------------------------------+ |
//   +---------------------------------------------------------+
//
// Design Minimization Phase
//   | +-Pack Rotamers Mover---------------------------------+ |
//   | | repack and design rotamers to explore sequence 	   | |
//   | | space  	                                           | |
//   | +-----------------------------------------------------+ |
//   | +-Minimization Mover----------------------------------+ |
//   | | energy minimize the current conformation            | |
//   | +-----------------------------------------------------+ |

//

// tracer - used to replace cout
static basic::Tracer TR("ODDM");

// application specific options
namespace oddm {
	// pert options
	RealOptionKey const mc_temp( "oddm::mc_temp" );
	RealOptionKey const pert_mc_temp( "oddm::pert_mc_temp" );
	RealOptionKey const pert_dock_rot_mag( "oddm::pert_dock_rot_mag" );
	RealOptionKey const pert_dock_trans_mag( "oddm::pert_dock_trans_mag" );
	RealOptionKey const pert_pep_small_temp( "oddm::pert_pep_small_temp" );
	RealOptionKey const pert_pep_small_H( "oddm::pert_pep_small_H" );
	RealOptionKey const pert_pep_small_L( "oddm::pert_pep_small_L" );
	RealOptionKey const pert_pep_small_E( "oddm::pert_pep_small_E" );
	RealOptionKey const pert_pep_shear_temp( "oddm::pert_pep_shear_temp" );
	RealOptionKey const pert_pep_shear_H( "oddm::pert_pep_shear_H" );
	RealOptionKey const pert_pep_shear_L( "oddm::pert_pep_shear_L" );
	RealOptionKey const pert_pep_shear_E( "oddm::pert_pep_shear_E" );

	IntegerOptionKey const pert_pep_num_rep( "oddm::pert_pep_num_rep" );
	IntegerOptionKey const pert_num( "oddm::pert_num" );
	IntegerOptionKey const dock_design_loop_num( "oddm::dock_design_loop_num" );

	BooleanOptionKey const final_design_min( "oddm::final_design_min" );
	BooleanOptionKey const use_soft_rep( "oddm::use_soft_rep" );
	BooleanOptionKey const mc_initial_pose( "oddm::mc_initial_pose" );
	BooleanOptionKey const oop_design_first( "oddm::oop_design_first" );

	BooleanOptionKey const pymol( "oddm::pymol" );
	BooleanOptionKey const keep_history( "oddm::keep_history" );

	// design options
	RealOptionKey const desn_mc_temp( "oddm::desn_mc_temp" );

	//IntegerVectorOptionKey const oop_positions( "oddm::oop_positions" );

}

class OopDockDesignMinimizeMover : public Mover {

	public:

		//default ctor
		OopDockDesignMinimizeMover(): Mover("OopDockDesignMinimizeMover"){}

		//default dtor
		virtual ~OopDockDesignMinimizeMover(){}

		//methods
		void setup_pert_foldtree( core::pose::Pose & pose);
		void setup_filter_stats();
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "OopDockDesignMinimizeMover"; }

};

typedef utility::pointer::owning_ptr< OopDockDesignMinimizeMover > OopDockDesignMinimizeMoverOP;
typedef utility::pointer::owning_ptr< OopDockDesignMinimizeMover const > OopDockDesignMinimizeMoverCOP;


int
main( int argc, char* argv[] )
{
    try {
	/*********************************************************************************************************************
	Common Setup
	**********************************************************************************************************************/

	// add application specific options to options system
	// There are far more options here than you will realistically need for a program of this complexity - but this gives you an idea of how to fine-grain option-control everything
	option.add( oddm::mc_temp, "The temperature to use for the outer loop of the ODDM protocol. Defaults to 1.0." ).def( 1.0 );
	option.add( oddm::pert_mc_temp, "The temperature to use for the pertubation phase of the ODDM protocol. Defaults to 0.8." ).def( 0.8 );
	option.add( oddm::pert_dock_rot_mag, "The rotation magnitude for the ridged body pertubation in the pertubation phase of the ODDM protocol. Defaults to 1.0." ).def( 1 );
	option.add( oddm::pert_dock_trans_mag, "The translation magnitude for the ridged body pertubation in the pertubation phase of the ODDM protocol. Defaults to 0.5." ).def( 0.5 );
	option.add( oddm::pert_pep_small_temp, "" ).def( 0.8 );
	option.add( oddm::pert_pep_shear_temp, "" ).def( 0.8 );

	option.add( oddm::pert_pep_small_H, "" ).def( 2.0 );
	option.add( oddm::pert_pep_small_L, "" ).def( 2.0 );
	option.add( oddm::pert_pep_small_E, "" ).def( 2.0 );
	option.add( oddm::pert_pep_shear_H, "" ).def( 2.0 );
	option.add( oddm::pert_pep_shear_L, "" ).def( 2.0 );
	option.add( oddm::pert_pep_shear_E, "" ).def( 2.0 );

	option.add( oddm::pert_pep_num_rep, "Number of small and shear iterations for the peptide" ).def( 100 );
	option.add( oddm::pert_num, "Number of iterations of perturbation loop per design" ).def(10);
	option.add( oddm::dock_design_loop_num, "Number of iterations of pertubation and design" ).def(10);

	option.add( oddm::final_design_min, "Do a final repack/design and minimization. Default true" ).def(true);
	option.add( oddm::use_soft_rep, "Use soft repulsion for pertubation and initial design. Default false" ).def(false);
	option.add( oddm::mc_initial_pose, "Allow initial pose to be considered as lowest energy pose. Default false" ).def(false);
	option.add( oddm::oop_design_first, "Design before pertubation (want when initial struct is aligned to hotspot)  Default false" ).def(false);

	option.add( oddm::pymol, "Set up pymol mover. Default false" ).def(false);
	option.add( oddm::keep_history, "Keep history in pymol. Requires oddm::pymol set to true. Default false" ).def(false);

	option.add( oddm::desn_mc_temp, "The temperature to use for the design/minimization phase of the ODDM protocol. Defaults to 0.8." ).def( 0.8 );

	//utility::vector1< core::Size > empty_vector(0);
	//option.add( oddm::oop_positions, "The positions of the first residues of oop rings" ).def( empty_vector );

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);

	//create mover instance
	OopDockDesignMinimizeMoverOP ODDM_mover( new OopDockDesignMinimizeMover() );

	ODDM_mover->setup_filter_stats();

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( ODDM_mover );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}//main

void
OopDockDesignMinimizeMover::apply(
	core::pose::Pose & pose
)
{
	// create score function
	//kdrew: old standard scoring function, using MM scoring function now because of NCAAs
	//scoring::ScoreFunctionOP score_fxn( getScoreFunction() );
	//scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::ScoreFunctionOP score_fxn = getScoreFunction();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	//kdrew: soft_rep score function
	//scoring::ScoreFunctionOP soft_score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::ScoreFunctionOP soft_score_fxn  = getScoreFunction();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*soft_score_fxn);
	soft_score_fxn->set_etable( FA_STANDARD_SOFT );

	scoring::ScoreFunctionOP pert_score_fxn;
	if( option[ oddm::use_soft_rep ].value() )
	{
    	pert_score_fxn = soft_score_fxn;
	}
	else
	{
    	pert_score_fxn = score_fxn;
	}

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	// get a fold tree suitable for docking (local helper function)
	setup_pert_foldtree( pose );


	//pose.conformation().show_residue_connections();

	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn, option[ oddm::mc_temp ].value() ) );

	/*********************************************************************************************************************
	Pertubation Phase
	**********************************************************************************************************************/

	// create a monte carlo object for the pertubation phase
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *pert_score_fxn, option[ oddm::pert_mc_temp ].value() ) );

	/*********************************************************
	Docking Setup
	**********************************************************/

	// create a rigid body mover to move the peptide around in the pocket
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, option[ oddm::pert_dock_rot_mag].value(),  option[ oddm::pert_dock_trans_mag].value()) );

	/*********************************************************
	Peptide Setup
	**********************************************************/

	// get peptide start and end positions
	Size pep_start( pose.conformation().chain_begin( 2 ) ); Size pep_end( pose.total_residue() );
	TR << "pep_start: " << pep_start << " pep_end: " << pep_end << std::endl;

	// create movemap for peptide
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb_true_range(pep_start, pep_end);

	////kdrew: automatically find oop positions
	utility::vector1< core::Size > oop_seq_positions = core::pose::ncbb::initialize_oops(pose);

	for( Size i = 1; i <= oop_seq_positions.size(); i++  )
	{
		pert_pep_mm->set_bb( oop_seq_positions[i], false );

		if( score_fxn->has_zero_weight( core::scoring::atom_pair_constraint ) )
		{
			score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		}

	}

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, option[ oddm::pert_pep_small_temp ].value(), 1 ) );
	pert_pep_small->angle_max( 'H', option[ oddm::pert_pep_small_H ].value() );
	pert_pep_small->angle_max( 'L', option[ oddm::pert_pep_small_L ].value() );
	pert_pep_small->angle_max( 'E', option[ oddm::pert_pep_small_E ].value() );

	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_mm, option[ oddm::pert_pep_shear_temp ].value(), 1 ) );
	pert_pep_shear->angle_max( 'H', option[ oddm::pert_pep_shear_H ].value() );
	pert_pep_shear->angle_max( 'L', option[ oddm::pert_pep_shear_L ].value() );
	pert_pep_shear->angle_max( 'E', option[ oddm::pert_pep_shear_E ].value() );

	// create random mover
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );
	//pert_pep_random->add_mover( pert_pep_shear, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, option[ oddm::pert_pep_num_rep ].value() ) );

	/*********************************************************
	OopPuck Setup
	**********************************************************/

	//oop::OopRandomSmallMoverOP opm_small( new oop::OopRandomSmallMover( option[ oddm::oop_positions ].value(), 2.0 ) );
	//oop::OopRandomPuckMoverOP opm_puck( new oop::OopRandomPuckMover( option[ oddm::oop_positions ].value() ) );
	oop::OopRandomSmallMoverOP opm_small( new oop::OopRandomSmallMover( oop_seq_positions, 2.0 ) );
	oop::OopRandomPuckMoverOP opm_puck( new oop::OopRandomPuckMover( oop_seq_positions ) );

	/******************************************************************************
	Rotamer Trials Setup
	*******************************************************************************/

	// create a task factory and task operations
	TaskFactoryOP pert_tf(new TaskFactory());
	pert_tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

	operation::ReadResfileOP pert_rrop( new operation::ReadResfile() );
	pert_rrop->default_filename();
	pert_tf->push_back( pert_rrop );

	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

    //operation::RestrictToInterfaceOP pert_rtio( new operation::RestrictToInterface(1, 2) ); //magic numbers: assume chains 1 and 2
	//pert_tf->push_back( pert_rtio );

	// create a rotamer trials mover
	simple_moves::RotamerTrialsMoverOP pert_rt(new simple_moves::EnergyCutRotamerTrialsMover( pert_score_fxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

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
	pert_sequence->add_mover( pert_rt );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	/*********************************************************************************************************************
	Design Min Phase
	**********************************************************************************************************************/
	/*********************************************************
	Design Setup
	**********************************************************/

	// create a task factory and task operations
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

	operation::ReadResfileOP desn_rrop( new operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );

	/*
	desn_tf->push_back( pert_rtio ); //not bothering to construct it a second time
	*/

	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover() );
	desn_pr->task_factory( desn_tf );
	desn_pr->score_function( pert_score_fxn );

	/*********************************************************
	Minimize Setup
	**********************************************************/

	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );
	//kdrew: set backbone of target false and backbone of oop true, decide whether to do this or not
	desn_mm->set_bb( false );
	desn_mm->set_bb_true_range( pep_start, pep_end );
	//desn_mm->set_bb( true );
	desn_mm->set_chi( true );
	desn_mm->set_jump( 1, true );

	//kdrew: loop through oop positions and set VZP and VYP torsions to minimize
	/*
	for( Size i = 1; i <= oop_seq_positions.size(); i++  )
	{
   		AtomID aidVZP( pose.residue(oop_seq_positions[i]).atom_index("VZP"), oop_seq_positions[i] );
		AtomID aidCYP( pose.residue(oop_seq_positions[i]).atom_index("CYP"), oop_seq_positions[i] );
		AtomID aidN( pose.residue(oop_seq_positions[i]).atom_index("N"), oop_seq_positions[i] );
		AtomID aidCA( pose.residue(oop_seq_positions[i]).atom_index("CA"), oop_seq_positions[i] );
		DOF_ID const vzp_tor_dof_id( pose.atom_tree().torsion_angle_dof_id( aidVZP, aidCYP, aidN, aidCA ) );
		desn_mm->set( vzp_tor_dof_id, true );

   		AtomID aidVYP( pose.residue(oop_seq_positions[i]+1).atom_index("VYP"), oop_seq_positions[i]+1 );
		AtomID aidCZP( pose.residue(oop_seq_positions[i]+1).atom_index("CZP"), oop_seq_positions[i]+1 );
		AtomID aidNpost( pose.residue(oop_seq_positions[i]+1).atom_index("N"), oop_seq_positions[i]+1 );
		AtomID aidCApost( pose.residue(oop_seq_positions[i]+1).atom_index("CA"), oop_seq_positions[i]+1 );
		DOF_ID const vyp_tor_dof_id( pose.atom_tree().torsion_angle_dof_id( aidVYP, aidCZP, aidNpost, aidCApost ) );
		desn_mm->set( vyp_tor_dof_id, true );
	}
	*/

	// create minimization mover
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min = new TaskAwareMinMover( desn_min, desn_tf );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a sequence mover to hold pack rotamers and minimization movers
	moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
	desn_sequence->add_mover( desn_pr );
	desn_sequence->add_mover( desn_ta_min );

	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/

	/*
	//Note you probably do not want a main loop like this in your protocol - it would be better to have most of the code preceding this as a mover's apply() function, so here you would have the mover/monte carlo statements but not the outer for loop and pose IO
	core::pose::Pose const copy_pose(pose);
	for ( Size i = 1; i <= Size( option[ OptionKeys::out::nstruct ].value() ) ; ++i ) {
		pose = copy_pose;
		pert_mc->reset(pose);
		std::stringstream filename_pdb;
		filename_pdb << "design_" << i << ".pdb";


	}
	*/
	TR << "Main loop..." << std::endl;

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

//kdrew: only turn on pymol observer in debug mode
//#ifndef NDEBUG
if( option[ oddm::pymol ].value() )
{
	protocols::moves::PyMolObserverOP pymover = protocols::moves::AddPyMolObserver(pose, option[ oddm::keep_history ].value() );
}
//#endif

	//pose.dump_pdb("pre_main_loop.pdb");
	for ( Size k = 1; k <= Size( option[ oddm::dock_design_loop_num ].value() ); ++k ) {


		pert_mc->reset(pose);

		//kdrew: a quick design/repack prior to pertubation, often the initial structure given is aligned to hotspot Ca Cb vector
		//kdrew: and do not want to perturb away until designed in hotspot residue
		if( k == 1 && option[ oddm::oop_design_first ].value() )
		{
			desn_sequence->apply( pose );
		}

		// pert loop
		for( Size j = 1; j <= Size( option[ oddm::pert_num ].value() ); ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
			curr_job->add_string_real_pair( "ENERGY_PERT (pert score)", (*pert_score_fxn)(pose) );
		}
		pert_mc->recover_low( pose );
		curr_job->add_string_real_pair( "ENERGY_PERT (pert score) recovered low", (*pert_score_fxn)(pose) );

		// design
		TR << "DESIGN: " << k << std::endl;
		desn_sequence->apply( pose );
		curr_job->add_string_real_pair( "ENERGY_DESN (hard score)", (*score_fxn)(pose) );

		//kdrew: reset mc after first cycle if not considering initial pose
		if( !option[ oddm::mc_initial_pose ].value() && k == 1 )
		{
			mc->reset(pose);
			TR<< "after mc->reset" << std::endl;
			mc->show_state();
		}

		TR<< "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		TR<< "post mc->boltzmann" << std::endl;
		mc->show_state();

	}//dock_design for loop

	mc->recover_low( pose );

	curr_job->add_string_real_pair( "ENERGY_FINAL (pert score) ", (*pert_score_fxn)(pose) );
	curr_job->add_string_real_pair( "ENERGY_FINAL (hard score) ", (*score_fxn)(pose) );

	TR << "Ending main loop..." << std::endl;

	TR << "Checking pose energy..." << std::endl;

		// create  MetricValues
	basic::MetricValue< core::Real > mv_sasa_complex;
	basic::MetricValue< core::Real > mv_sasa_seperated;
	basic::MetricValue< utility::vector1< core::Size > > mv_unsat_res_complex;
	basic::MetricValue< utility::vector1< core::Size > > mv_unsat_res_seperated;
	basic::MetricValue< core::Real > mv_pack_complex;
	basic::MetricValue< core::Real > mv_pack_seperated;

	basic::MetricValue< core::Real > mv_repack_sasa_seperated;
	basic::MetricValue< utility::vector1< core::Size > > mv_repack_unsat_res_seperated;
	basic::MetricValue< core::Real > mv_repack_pack_seperated;
	core::Real repack_energy_seperated;
	core::Real repack_hbond_ener_sum_seperated;

	core::Real energy_complex;
	core::Real energy_seperated;
	core::Real hbond_ener_sum_complex;
	core::Real hbond_ener_sum_seperated;

	// calc energy
	energy_complex = (*score_fxn)(pose);

	TR << "Energy less than cutoff, doing final design and running filters..." << std::endl;

	if ( option[ oddm::final_design_min].value() )
	{
		// get packer task from task factory
		PackerTaskOP final_desn_pt( *(desn_tf->create_task_and_apply_taskoperations( pose )) );

		// add extra chi and extra chi cut off to pt
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			final_desn_pt->nonconst_residue_task( i ).or_ex1( true );
			final_desn_pt->nonconst_residue_task( i ).or_ex2( true );
			final_desn_pt->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
		}

		// create a pack rotamers mover for the final design
		simple_moves::PackRotamersMoverOP final_desn_pr( new simple_moves::PackRotamersMover(score_fxn, final_desn_pt, 10 ) );
		//final_desn_pr->packer_task( final_desn_pt );
		//final_desn_pr->score_function( score_fxn );
		//final_desn_pr->nloop( 10 );

		// design with final pr mover
		final_desn_pr->apply( pose );

		// create move map for minimization
		kinematics::MoveMapOP final_min_mm( new kinematics::MoveMap() );
		final_min_mm->set_bb( true );
		final_min_mm->set_chi( true );
		final_min_mm->set_jump( 1, true );

		// create minimization mover
		simple_moves::MinMoverOP final_min( new simple_moves::MinMover( final_min_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );
		// final min (okay to use ta min here)
		final_min->apply( pose );
	}

	// make copy of pose to calc stats
	Pose stats_pose( pose );

	// complex stats
	energy_complex = (*score_fxn)(stats_pose);
	stats_pose.metric("sasa","total_sasa",mv_sasa_complex);
	stats_pose.metric("unsat", "residue_bur_unsat_polars", mv_unsat_res_complex);
	utility::vector1< core::Size > const unsat_res_complex(mv_unsat_res_complex.value());
	stats_pose.metric( "pack", "total_packstat", mv_pack_complex );
	scoring::EnergyMap complex_emap( stats_pose.energies().total_energies() );
	hbond_ener_sum_complex = complex_emap[ hbond_sr_bb ] + complex_emap[ hbond_lr_bb ] + complex_emap[ hbond_bb_sc ] + complex_emap[ hbond_sc ];

	// seperate designed chain from other chains
	protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( pose, 1 ) ); // HARDCODED JUMP NUMBER
	translate->step_size( 1000.0 );
	translate->apply( stats_pose );
	//stats_pose.dump_pdb("stats_trans1000.pdb");

	Pose repack_stats_pose( stats_pose );

	//kdrew: probably should repack and minimize here after separation
	TaskFactoryOP tf(new TaskFactory());
	tf->push_back( new core::pack::task::operation::InitializeFromCommandline );
	//kdrew: do not do design, makes NATAA if res file is not specified
	operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
	tf->push_back( rtrp );
	simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( score_fxn );
	packer->apply( repack_stats_pose );

	// create move map for minimization
	kinematics::MoveMapOP separate_min_mm( new kinematics::MoveMap() );
	separate_min_mm->set_bb( true );
	separate_min_mm->set_chi( true );
	separate_min_mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP separate_min( new simple_moves::MinMover( separate_min_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );
	// final min (okay to use ta min here)
	separate_min->apply( repack_stats_pose );

	// seperate stats
	energy_seperated = (*score_fxn)(stats_pose);
	repack_energy_seperated = (*score_fxn)(repack_stats_pose);
	stats_pose.metric("sasa","total_sasa",mv_sasa_seperated);
	repack_stats_pose.metric("sasa","total_sasa",mv_repack_sasa_seperated);
	stats_pose.metric("unsat", "residue_bur_unsat_polars", mv_unsat_res_seperated);
	repack_stats_pose.metric("unsat", "residue_bur_unsat_polars", mv_repack_unsat_res_seperated);
	utility::vector1< core::Size > const unsat_res_seperated(mv_unsat_res_seperated.value());
	stats_pose.metric( "pack", "total_packstat", mv_pack_seperated );
	repack_stats_pose.metric( "pack", "total_packstat", mv_repack_pack_seperated );
	scoring::EnergyMap seperated_emap( stats_pose.energies().total_energies() );
	hbond_ener_sum_seperated = seperated_emap[ hbond_sr_bb ] + seperated_emap[ hbond_lr_bb ] + seperated_emap[ hbond_bb_sc ] + seperated_emap[ hbond_sc ];
	scoring::EnergyMap repack_seperated_emap( repack_stats_pose.energies().total_energies() );
	repack_hbond_ener_sum_seperated = repack_seperated_emap[ hbond_sr_bb ] + repack_seperated_emap[ hbond_lr_bb ] + repack_seperated_emap[ hbond_bb_sc ] + repack_seperated_emap[ hbond_sc ];

	// add values to job so that they will be output in the pdb
	curr_job->add_string_real_pair( "ENERGY_COMPLEX:\t\t", energy_complex );
	curr_job->add_string_real_pair( "ENERGY_SEPERATE:\t\t", energy_seperated );
	curr_job->add_string_real_pair( "ENERGY_DIFF:\t\t", energy_complex - energy_seperated );
	curr_job->add_string_real_pair( "REPACK_ENERGY_SEPERATE:\t\t", repack_energy_seperated );
	curr_job->add_string_real_pair( "REPACK_ENERGY_DIFF:\t\t", energy_complex - repack_energy_seperated );

	curr_job->add_string_real_pair( "SASA_COMPLEX:\t\t", mv_sasa_complex.value() );
	curr_job->add_string_real_pair( "SASA_SEPERATE:\t\t", mv_sasa_seperated.value() );
	curr_job->add_string_real_pair( "SASA_DIFF:\t\t", mv_sasa_complex.value() - mv_sasa_seperated.value() );
	curr_job->add_string_real_pair( "REPACK_SASA_SEPERATE:\t\t", mv_repack_sasa_seperated.value() );
	curr_job->add_string_real_pair( "REPACK_SASA_DIFF:\t\t", mv_sasa_complex.value() - mv_repack_sasa_seperated.value() );

	curr_job->add_string_real_pair( "HB_ENER_COMPLEX:\t\t", hbond_ener_sum_complex );
	curr_job->add_string_real_pair( "HB_ENER_SEPERATE:\t\t", hbond_ener_sum_seperated );
	curr_job->add_string_real_pair( "HB_ENER_DIFF:\t\t", hbond_ener_sum_complex - hbond_ener_sum_seperated );
	curr_job->add_string_real_pair( "REPACK_HB_ENER_SEPERATE:\t\t", repack_hbond_ener_sum_seperated );
	curr_job->add_string_real_pair( "REPACK_HB_ENER_DIFF:\t\t", hbond_ener_sum_complex - repack_hbond_ener_sum_seperated );

	curr_job->add_string_real_pair( "PACK_COMPLEX:\t\t", mv_pack_complex.value() );
	curr_job->add_string_real_pair( "PACK_SEPERATE:\t\t", mv_pack_seperated.value() );
	curr_job->add_string_real_pair( "PACK_DIFF:\t\t", mv_pack_complex.value() - mv_pack_seperated.value() );
	curr_job->add_string_real_pair( "REPACK_PACK_SEPERATE:\t\t", mv_repack_pack_seperated.value() );
	curr_job->add_string_real_pair( "REPACK_PACK_DIFF:\t\t", mv_pack_complex.value() - mv_repack_pack_seperated.value() );

}
// this only works for two chains and assumes the protein is first and the peptide is second
// inspired by protocols/docking/DockingProtocol.cc
void
OopDockDesignMinimizeMover::setup_pert_foldtree(
	core::pose::Pose & pose
)
{
	using namespace kinematics;

	// get current fold tree
	FoldTree f( pose.fold_tree() );
	f.clear();

	// get the start and end for both chains
	Size pro_start( pose.conformation().chain_begin( 1 ) );
	Size pro_end( pose.conformation().chain_end( 1 ) );
	Size pep_start( pose.conformation().chain_begin( 2 ) );
	Size pep_end( pose.conformation().chain_end( 2 ) );

	// get jump positions based on the center of mass of the chains
	Size dock_jump_pos_pro( geometry::residue_center_of_mass( pose, pro_start, pro_end ) );
	Size dock_jump_pos_pep( geometry::residue_center_of_mass( pose, pep_start, pep_end ) );

	// build fold tree
	Size jump_index( f.num_jump() + 1 );
	//-1 is a magic number for PEPTIDE EDGE.  There is a constant defined with the fold tree that should have been used here.
	f.add_edge( pro_start, dock_jump_pos_pro, -1 );
	f.add_edge( dock_jump_pos_pro, pro_end, -1 );
	f.add_edge( pep_start, dock_jump_pos_pep, -1);
	f.add_edge( dock_jump_pos_pep, pep_end, -1 );
	f.add_edge( dock_jump_pos_pro, dock_jump_pos_pep, jump_index );

	// set pose foldtree to foldtree we just created
	f.reorder(1);
	f.check_fold_tree();
	assert( f.check_fold_tree() );

	std::cout << "AFTER: " << f << std::endl;

	pose.fold_tree( f );
}

void
OopDockDesignMinimizeMover::setup_filter_stats()
{
	/*********************************************************************************************************************
	Filter / Stats Setup
	*********************************************************************************************************************/

	// create and register sasa calculator
	pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculator() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

	// create and register hb calculator
	pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator( new pose_metric_calculators::NumberHBondsCalculator() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );

	// create and register unsat calculator
	pose::metrics::PoseMetricCalculatorOP unsat_calculator( new pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds") ) ;
	pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );

	// create and register packstat calculator
	pose::metrics::PoseMetricCalculatorOP pack_calcculator( new pose_metric_calculators::PackstatCalculator() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "pack", pack_calcculator );

}


