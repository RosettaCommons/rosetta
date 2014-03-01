// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/moves/Mover.hh>
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
#include <protocols/simple_moves/hbs/HbsMover.hh>
#include <protocols/simple_moves/hbs/HbsRandomSmallMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/ncbb/NcbbDockDesignProtocol.hh>
#include <protocols/ncbb/NcbbDockDesignProtocolCreator.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
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
using namespace protocols::simple_moves::hbs;
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
static basic::Tracer TR("NDDP");

namespace protocols {
namespace ncbb {

NcbbDockDesignProtocol::NcbbDockDesignProtocol():
		mc_temp_ ( 1.0),
		pert_mc_temp_ (0.8),
		pert_dock_rot_mag_ (1.0),
		pert_dock_trans_mag_ (0.5),
		pert_pep_small_temp_ (0.8),
		pert_pep_small_H_ (2.0),
		pert_pep_small_L_ (2.0),
		pert_pep_small_E_ (2.0),
		pert_pep_shear_temp_ (0.8),
		pert_pep_shear_H_ (2.0),
		pert_pep_shear_L_ (2.0),
		pert_pep_shear_E_ (2.0),

		pert_pep_num_rep_ (100),
		pert_num_ (10),
		dock_design_loop_num_ (10),

		no_design_(false),
		final_design_min_ (true),
		use_soft_rep_ (false),
		mc_initial_pose_ (false),
		ncbb_design_first_ (false),

		pymol_ (false),
		keep_history_ (false)
{ 
		Mover::type("NcbbDockDesignProtocol");
	
		score_fxn_ = getScoreFunction(); 
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
}

NcbbDockDesignProtocol::NcbbDockDesignProtocol(
	scoring::ScoreFunctionOP score_function,
	core::Real const mc_temp,
	core::Real const pert_mc_temp,
	core::Real const pert_dock_rot_mag,
	core::Real const pert_dock_trans_mag,
	core::Real const pert_pep_small_temp,
	core::Real const pert_pep_small_H,
	core::Real const pert_pep_small_L,
	core::Real const pert_pep_small_E,
	core::Real const pert_pep_shear_temp,
	core::Real const pert_pep_shear_H,
	core::Real const pert_pep_shear_L,
	core::Real const pert_pep_shear_E,

	core::Size const pert_pep_num_rep,
	core::Size const pert_num,
	core::Size const dock_design_loop_num,

	bool const no_design,
	bool const final_design_min,
	bool const use_soft_rep,
	bool const mc_initial_pose,
	bool const ncbb_design_first,

	bool const pymol,
	bool const keep_history

): Mover("NcbbDockDesignProtocol"),
	score_fxn_(score_function),
	mc_temp_ (mc_temp),
	pert_mc_temp_(pert_mc_temp),
	pert_dock_rot_mag_(pert_dock_rot_mag),
	pert_dock_trans_mag_(pert_dock_trans_mag),
	pert_pep_small_temp_(pert_pep_small_temp),
	pert_pep_small_H_(pert_pep_small_H),
	pert_pep_small_L_(pert_pep_small_L),
	pert_pep_small_E_(pert_pep_small_E),
	pert_pep_shear_temp_(pert_pep_shear_temp),
	pert_pep_shear_H_(pert_pep_shear_H),
	pert_pep_shear_L_(pert_pep_shear_L),
	pert_pep_shear_E_(pert_pep_shear_E),
	pert_pep_num_rep_(pert_pep_num_rep),
	pert_num_(pert_num),
	dock_design_loop_num_(dock_design_loop_num),
	no_design_(no_design),
	final_design_min_(final_design_min),
	use_soft_rep_(use_soft_rep),
	mc_initial_pose_(mc_initial_pose),
	ncbb_design_first_(ncbb_design_first),
	pymol_(pymol),
	keep_history_(keep_history) 
{	
}

NcbbDockDesignProtocol::NcbbDockDesignProtocol(
	scoring::ScoreFunctionOP score_function, 
	core::Real const mc_temp,
	core::Real const pert_dock_rot_mag,
	core::Real const pert_dock_trans_mag,
	core::Size const dock_design_loop_num,
	bool const no_design,
	bool const final_design_min,
	bool const pymol,
	bool const keep_history
): Mover("NcbbDockDesignProtocol"),
	score_fxn_(score_function),
	mc_temp_ (mc_temp),
	pert_mc_temp_ (0.8),
	pert_dock_rot_mag_ (pert_dock_rot_mag),
	pert_dock_trans_mag_ (pert_dock_trans_mag),
	pert_pep_small_temp_ (0.8),
	pert_pep_small_H_ (2.0),
	pert_pep_small_L_ (2.0),
	pert_pep_small_E_ (2.0),
	pert_pep_shear_temp_ (0.8),
	pert_pep_shear_H_ (2.0),
	pert_pep_shear_L_ (2.0),
	pert_pep_shear_E_ (2.0),

	pert_pep_num_rep_ (100),
	pert_num_ (10),
	dock_design_loop_num_(dock_design_loop_num),

	no_design_(no_design),
	final_design_min_(final_design_min),
	use_soft_rep_ (false),
	mc_initial_pose_ (false),
	ncbb_design_first_ (false),

	pymol_(pymol),
	keep_history_(keep_history)
{}

void
NcbbDockDesignProtocol::apply(
	core::pose::Pose & pose
)
{
	scoring::ScoreFunctionOP soft_score_fxn  = getScoreFunction();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*soft_score_fxn);
	soft_score_fxn->set_etable( FA_STANDARD_SOFT );

	scoring::ScoreFunctionOP pert_score_fxn;

	if ( use_soft_rep_ )
		pert_score_fxn = soft_score_fxn;
	else
		pert_score_fxn = score_fxn_;

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	// get a fold tree suitable for docking (local helper function)
	setup_pert_foldtree( pose );

	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn_, mc_temp_ ) );

	/*********************************************************************************************************************
	Pertubation Phase
	**********************************************************************************************************************/
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *pert_score_fxn, pert_mc_temp_ ) );

	/*********************************************************
	Docking Setup
	**********************************************************/
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, pert_dock_rot_mag_,  pert_dock_trans_mag_) );

	/*********************************************************
	Peptide Setup
	**********************************************************/
	Size pep_start( pose.conformation().chain_begin( 2 ) ); Size pep_end( pose.total_residue() );
	TR << "pep_start: " << pep_start << " pep_end: " << pep_end << std::endl;

	// create movemap for peptide
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb_true_range(pep_start, pep_end);

	////kdrew: automatically find ncbb positions
	utility::vector1< core::Size > ncbb_seq_positions = core::pose::ncbb::initialize_ncbbs(pose);

	////awatkins: initialize specific vectors for each supported patch type, too
	utility::vector1< core::Size > oop_seq_positions = core::pose::ncbb::initialize_oops(pose);
	utility::vector1< core::Size > hbs_seq_positions = core::pose::ncbb::initialize_hbs(pose);
	
	for( Size i = 1; i <= ncbb_seq_positions.size(); ++i  )
	{
		pert_pep_mm->set_bb( ncbb_seq_positions[i], false );

		if ( score_fxn_->has_zero_weight( core::scoring::atom_pair_constraint ) )
			score_fxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	}
	

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, pert_pep_small_temp_, 1 ) );
	pert_pep_small->angle_max( 'H', pert_pep_small_H_ );
	pert_pep_small->angle_max( 'L', pert_pep_small_L_ );
	pert_pep_small->angle_max( 'E', pert_pep_small_E_ );

	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_mm, pert_pep_shear_temp_ , 1 ) );
	pert_pep_shear->angle_max( 'H', pert_pep_shear_H_ );
	pert_pep_shear->angle_max( 'L', pert_pep_shear_L_ );
	pert_pep_shear->angle_max( 'E', pert_pep_shear_E_ );

	// create random mover
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );
	//pert_pep_random->add_mover( pert_pep_shear, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, pert_pep_num_rep_ ) );

	//simple_moves::hbs::HbsRandomSmallMoverOP hpm_small( new simple_moves::hbs::HbsRandomSmallMover( hbs_seq_positions, 2.0 ) );
	simple_moves::oop::OopRandomSmallMoverOP opm_small( new simple_moves::oop::OopRandomSmallMover( oop_seq_positions, 2.0 ) );
	simple_moves::oop::OopRandomPuckMoverOP opm_puck( new simple_moves::oop::OopRandomPuckMover( oop_seq_positions ) );

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

	// create a rotamer trials mover
	simple_moves::RotamerTrialsMoverOP pert_rt(new simple_moves::EnergyCutRotamerTrialsMover( pert_score_fxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

	/*********************************************************
	Common Setup
	**********************************************************/

	// create a random mover to hold the docking, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_repeat, 0.5 );
	//pert_random->add_mover( hpm_small, 0.5 );
	if (oop_seq_positions.size() > 0) {
		pert_random->add_mover( opm_small, 0.5 );
		pert_random->add_mover( opm_puck, 0.1 );
	}

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_random );
	pert_sequence->add_mover( pert_rt );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	/*********************************************************
	Design Setup
	**********************************************************/
	// create a task factory and task operations
	TaskFactoryOP desn_tf( new TaskFactory() );
	desn_tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

	operation::ReadResfileOP desn_rrop( new operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );

	if( no_design_ )
	{
		core::pack::task::operation::RestrictToRepackingOP rtrp( new core::pack::task::operation::RestrictToRepacking() );
		desn_tf->push_back( rtrp );
	}

	// create a pack rotamers mover
	simple_moves::PackRotamersMoverOP desn_pack_rotamers( new simple_moves::PackRotamersMover() );
	desn_pack_rotamers->task_factory( desn_tf );
	desn_pack_rotamers->score_function( pert_score_fxn );

	/*********************************************************
	Minimize Setup
	**********************************************************/
	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );
	//kdrew: set backbone of target false and backbone of ncbb true, decide whether to do this or not
	desn_mm->set_bb( false );
	desn_mm->set_bb_true_range( pep_start, pep_end );
	desn_mm->set_chi( true );
	desn_mm->set_jump( 1, true );
	
	// create minimization mover
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	// definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min = new TaskAwareMinMover( desn_min, desn_tf );

	/*********************************************************
	Common Setup
	**********************************************************/
	moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
	desn_sequence->add_mover( desn_pack_rotamers );
	desn_sequence->add_mover( desn_ta_min );

	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/
	TR << "Main loop..." << std::endl;

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	if( pymol_ )
		protocols::moves::PyMolObserverOP pymover = protocols::moves::AddPyMolObserver(pose, keep_history_ );

	//pose.dump_pdb("pre_main_loop.pdb");
	for ( Size k = 1; k <= Size( dock_design_loop_num_ ); ++k ) {
		pert_mc->reset(pose);

		//kdrew: a quick design/repack prior to pertubation, often the initial structure given is aligned to hotspot Ca Cb vector
		//kdrew: and do not want to perturb away until designed in hotspot residue
		if( k == 1 && ncbb_design_first_ )
			desn_sequence->apply( pose );

		// Perturbation phase - loop
		for( Size j = 1; j <= Size( pert_num_ ); ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
			curr_job->add_string_real_pair( "ENERGY_PERT (pert score)", (*pert_score_fxn)(pose) );
		}
		pert_mc->recover_low( pose );
		curr_job->add_string_real_pair( "ENERGY_PERT (pert score) recovered low", (*pert_score_fxn)(pose) );

		// Design phase
		TR << "DESIGN: " << k << std::endl;
		desn_sequence->apply( pose );
		curr_job->add_string_real_pair( "ENERGY_DESN (hard score)", (*score_fxn_)(pose) );

		//kdrew: reset mc after first cycle if not considering initial pose
		if( !mc_initial_pose_  && k == 1 )
		{
			mc->reset(pose);
			TR << "after mc->reset" << std::endl;
			mc->show_state();
		}

		TR << "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		TR << "post mc->boltzmann" << std::endl;
		mc->show_state();
	}//dock_design for loop

	mc->recover_low( pose );
	curr_job->add_string_real_pair( "ENERGY_FINAL (pert score) ", (*pert_score_fxn)(pose) );
	curr_job->add_string_real_pair( "ENERGY_FINAL (hard score) ", (*score_fxn_)(pose) );

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
	energy_complex = (*score_fxn_)(pose);

	TR << "Energy less than cutoff, doing final design and running filters..." << std::endl;

	if ( final_design_min_ )
	{
		// get packer task from task factory
		PackerTaskOP final_desn_packer_task( *(desn_tf->create_task_and_apply_taskoperations( pose )) );

		// add extra chi and extra chi cut off to pt
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			final_desn_packer_task->nonconst_residue_task( i ).or_ex1( true );
			final_desn_packer_task->nonconst_residue_task( i ).or_ex2( true );
			final_desn_packer_task->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
		}

		// create a pack rotamers mover for the final design
		simple_moves::PackRotamersMoverOP final_desn_pack_rotamers( new simple_moves::PackRotamersMover(score_fxn_, final_desn_packer_task, 10 ) );

		// design with final pr mover
		final_desn_pack_rotamers->apply( pose );

		// create move map for minimization
		kinematics::MoveMapOP final_min_mm( new kinematics::MoveMap() );
		final_min_mm->set_bb( true );
		final_min_mm->set_chi( true );
		final_min_mm->set_jump( 1, true );

		// create minimization mover
		simple_moves::MinMoverOP final_min( new simple_moves::MinMover( final_min_mm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );
		// final min (okay to use ta min here)
		final_min->apply( pose );
	}

	// make copy of pose to calc stats
	Pose stats_pose( pose );

	// complex stats
	energy_complex = (*score_fxn_)(stats_pose);
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
	packer->score_function( score_fxn_ );
	packer->apply( repack_stats_pose );

	// create move map for minimization
	kinematics::MoveMapOP separate_min_mm( new kinematics::MoveMap() );
	separate_min_mm->set_bb( true );
	separate_min_mm->set_chi( true );
	separate_min_mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP separate_min( new simple_moves::MinMover( separate_min_mm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );
	// final min (okay to use ta min here)
	separate_min->apply( repack_stats_pose );

	// seperate stats
	energy_seperated = (*score_fxn_)(stats_pose);
	repack_energy_seperated = (*score_fxn_)(repack_stats_pose);
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
NcbbDockDesignProtocol::setup_pert_foldtree(
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
	f.add_edge( pro_start, dock_jump_pos_pro, Edge::PEPTIDE );
	f.add_edge( dock_jump_pos_pro, pro_end, Edge::PEPTIDE );
	f.add_edge( pep_start, dock_jump_pos_pep, Edge::PEPTIDE );
	f.add_edge( dock_jump_pos_pep, pep_end, Edge::PEPTIDE );
	f.add_edge( dock_jump_pos_pro, dock_jump_pos_pep, jump_index );

	// set pose foldtree to foldtree we just created
	f.reorder(1);
	f.check_fold_tree();
	assert( f.check_fold_tree() );

	std::cout << "AFTER: " << f << std::endl;

	pose.fold_tree( f );
}

void
NcbbDockDesignProtocol::setup_filter_stats()
{
	// create and register sasa calculator
	pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
	if (!pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ))
		pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

	// create and register hb calculator
	pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator( new pose_metric_calculators::NumberHBondsCalculator() );
	if (!pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "num_hbonds" ))
		pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );

	// create and register unsat calculator
	pose::metrics::PoseMetricCalculatorOP unsat_calculator( new pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds") ) ;
	if (!pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "unsat" ))
		pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );

	// create and register packstat calculator
	pose::metrics::PoseMetricCalculatorOP pack_calculator( new pose_metric_calculators::PackstatCalculator() );
	if (!pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "pack" ))
		pose::metrics::CalculatorFactory::Instance().register_calculator( "pack", pack_calculator );
}


protocols::moves::MoverOP
NcbbDockDesignProtocol::clone() const
{
	return new NcbbDockDesignProtocol (
		score_fxn_,
  		mc_temp_,
  		pert_mc_temp_,
  		pert_dock_rot_mag_,
  		pert_dock_trans_mag_,
  		pert_pep_small_temp_,
  		pert_pep_small_H_,
  		pert_pep_small_L_,
  		pert_pep_small_E_,
  		pert_pep_shear_temp_,
  		pert_pep_shear_H_,
  		pert_pep_shear_L_,
  		pert_pep_shear_E_,
  		pert_pep_num_rep_,
  		pert_num_,
  		dock_design_loop_num_,
  		no_design_,
  		final_design_min_,
  		use_soft_rep_,
  		mc_initial_pose_,
  		ncbb_design_first_,
  		pymol_,
		keep_history_
		);
}

void
NcbbDockDesignProtocol::parse_my_tag
( 
	utility::tag::TagCOP tag, 
	basic::datacache::DataMap &data, 
	protocols::filters::Filters_map const &, 
	protocols::moves::Movers_map const &, 
	core::pose::Pose const &
) {

	if(tag->hasOption( "scorefxn"))
	{
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn" ) );
		if ( ! data.has( "scorefxn", scorefxn_key ) )
			throw utility::excn::EXCN_RosettaScriptsOption("ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		score_fxn_ = data.get< core::scoring::ScoreFunction* >( "scorefxns", scorefxn_key );
	}

	if(tag->hasOption( "mc_temp"))
		this->mc_temp_ = tag->getOption<core::Real>("mc_temp", mc_temp_);
	else
		mc_temp_ = 1.0;
 
	if(tag->hasOption( "pert_mc_temp"))
		pert_mc_temp_ = tag->getOption<core::Real>("pert_mc_temp", pert_mc_temp_);
	else
		pert_mc_temp_ = 0.8;

	if(tag->hasOption( "pert_dock_rot_mag"))
		pert_dock_rot_mag_ = tag->getOption<core::Real>("pert_dock_rot_mag", pert_dock_rot_mag_);
	else
		pert_dock_rot_mag_ = 1.0;
 
	if(tag->hasOption( "pert_dock_trans_mag"))
		pert_dock_trans_mag_ = tag->getOption<core::Real>("pert_dock_trans_mag", pert_dock_trans_mag_);
 	else
		pert_dock_trans_mag_ = 0.5;
 
	if(tag->hasOption( "pert_pep_small_temp"))
		pert_pep_small_temp_ = tag->getOption<core::Real>("pert_pep_small_temp", pert_pep_small_temp_);
 	else
		pert_pep_small_temp_ = 0.8;
 
	if(tag->hasOption( "pert_pep_small_H"))
		pert_pep_small_H_ = tag->getOption<core::Real>("pert_pep_small_H", pert_pep_small_H_);
	else
		pert_pep_small_H_ = 2.0;
 
	if(tag->hasOption( "pert_pep_small_L"))
		pert_pep_small_L_ = tag->getOption<core::Real>("pert_pep_small_L", pert_pep_small_L_);
 	else
		pert_pep_small_L_ = 2.0;
 
	if(tag->hasOption( "pert_pep_small_E"))
		pert_pep_small_E_ = tag->getOption<core::Real>("pert_pep_small_E", pert_pep_small_E_);
 	else
		pert_pep_small_E_ = 2.0;
 
	if(tag->hasOption( "pert_pep_shear_temp"))
		pert_pep_shear_temp_ = tag->getOption<core::Real>("pert_pep_shear_temp", pert_pep_shear_temp_);
 	else
		pert_pep_shear_temp_ = 0.8;

	if(tag->hasOption( "pert_pep_shear_H"))
		pert_pep_shear_H_ = tag->getOption<core::Real>("pert_pep_shear_H", pert_pep_shear_H_);
 	else
		pert_pep_shear_H_ = 2.0;
 
	if(tag->hasOption( "pert_pep_shear_L"))
		pert_pep_shear_L_ = tag->getOption<core::Real>("pert_pep_shear_L", pert_pep_shear_L_);
 	else
		pert_pep_shear_L_ = 2.0;

	if(tag->hasOption( "pert_pep_shear_E"))
		pert_pep_shear_E_ = tag->getOption<core::Real>("pert_pep_shear_E", pert_pep_shear_E_);
 	else
		pert_pep_shear_E_ = 2.0;
 
	if(tag->hasOption( "pert_pep_num_rep"))
		pert_pep_num_rep_ = tag->getOption<core::Size>("pert_pep_num_rep", pert_pep_num_rep_);
 	else
		pert_pep_num_rep_ = 100;

	if(tag->hasOption( "pert_num"))
		pert_num_ = tag->getOption<core::Size>("pert_num", pert_num_);
 	else
		pert_num_ = 10;

	if(tag->hasOption( "dock_design_loop_num"))
		dock_design_loop_num_ = tag->getOption<core::Size>("dock_design_loop_num", dock_design_loop_num_);
 	else
		dock_design_loop_num_ = 10;

	if(tag->hasOption( "no_design"))
		no_design_ = tag->getOption<bool>("no_design", no_design_);
 	else
		no_design_ = false;

	if(tag->hasOption( "final_design_min"))
		final_design_min_ = tag->getOption<bool>("final_design_min", final_design_min_);
 	else
		final_design_min_ = true;
 
	if(tag->hasOption( "use_soft_rep"))
		use_soft_rep_ = tag->getOption<bool>("use_soft_rep", use_soft_rep_);
 	else
		use_soft_rep_ = false;
 
	if(tag->hasOption( "mc_initial_pose"))
		mc_initial_pose_ = tag->getOption<bool>("mc_initial_pose", mc_initial_pose_);
 	else
		mc_initial_pose_ = false;

	if(tag->hasOption( "ncbb_design_first"))
		ncbb_design_first_ = tag->getOption<bool>("ncbb_design_first", ncbb_design_first_);
 	else
		ncbb_design_first_ = false;
 
	if(tag->hasOption( "pymol"))
		pymol_ = tag->getOption<bool>("pymol", pymol_);
 	else
		pymol_ = false;
 
	if(tag->hasOption( "keep_history"))
		keep_history_ = tag->getOption<bool>("keep_history", keep_history_);
 	else
		keep_history_ = false;
}

// MoverCreator
moves::MoverOP NcbbDockDesignProtocolCreator::create_mover() const {
        return new NcbbDockDesignProtocol();
}

std::string NcbbDockDesignProtocolCreator::keyname() const {
        return NcbbDockDesignProtocolCreator::mover_name();
}

std::string NcbbDockDesignProtocolCreator::mover_name(){
        return "NcbbDockDesign";
}


}
}
