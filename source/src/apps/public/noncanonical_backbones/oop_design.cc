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
	//RealOptionKey const mc_temp( "oop_design::mc_temp" );
	//RealOptionKey const pert_dock_rot_mag( "oop_design::pert_dock_rot_mag" );
	//RealOptionKey const pert_dock_trans_mag( "oop_design::pert_dock_trans_mag" );

	//IntegerOptionKey const pert_pep_num_rep( "oop_design::pert_pep_num_rep" );
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
		void setup_pert_foldtree( core::pose::Pose & pose);
		void setup_filter_stats();
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
	/*
	option.add( oop_design::mc_temp, "The temperature to use for the outer loop of the ODDM protocol. Defaults to 1.0." ).def( 1.0 );
	option.add( oop_design::pert_dock_rot_mag, "The rotation magnitude for the ridged body pertubation in the pertubation phase of the ODDM protocol. Defaults to 1.0." ).def( 1 );
	option.add( oop_design::pert_dock_trans_mag, "The translation magnitude for the ridged body pertubation in the pertubation phase of the ODDM protocol. Defaults to 0.5." ).def( 0.5 );
	*/

	//option.add( oop_design::pert_pep_num_rep, "Number of small and shear iterations for the peptide" ).def( 100 );
	option.add( oop_design::pert_num, "Number of iterations of perturbation loop per design" ).def(10);
	option.add( oop_design::design_loop_num, "Number of iterations of pertubation and design" ).def(10);

	utility::vector1< core::Size > empty_vector(0);
	option.add( oop_design::oop_design_positions, "Positions of oop to design" ).def( empty_vector );

	// init command line options
	//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
	devel::init(argc, argv);

	//create mover instance
	OopDesignMoverOP OD_mover( new OopDesignMover() );

	OD_mover->setup_filter_stats();

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
	//scoring::ScoreFunctionOP score_fxn = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	// get a fold tree suitable for docking (local helper function)
	setup_pert_foldtree( pose );

	// create a monte carlo object for the full cycle
	//moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn, option[ oop_design::mc_temp ].value() ) );
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

	for( Size i = 1; i <= oop_seq_positions.size(); i++  )
	{
		pert_pep_mm->set_bb( oop_seq_positions[i], false );

		if( score_fxn->has_zero_weight( core::scoring::atom_pair_constraint ) )
		{
			score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		}

	}

	//utility::vector1< core::Size > oop_seq_positions;
	//for ( Size i = 1; i <= pose.total_residue(); ++i )
	//{
	//	if( pose.residue(i).has_variant_type(chemical::OOP_PRE) == 1)
	//	{
	//		using basic::options::option;
	//		using namespace basic::options::OptionKeys;
	//		oop_seq_positions.push_back( i );
	//		//kdrew: set up constraints, if not detected at readin
	//		if ( ! option[ in::detect_oops ].user() )
	//		{
	//			core::pose::ncbb::add_oop_constraint(pose, i);
	//		}
	//		//kdrew: do not use small/shear mover on oop positions, use oop mover instead
	//		pert_pep_mm->set_bb( i, false );

	//		if( score_fxn->has_zero_weight( core::scoring::atom_pair_constraint ) )
	//		{
	//			score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	//		}
	//	}
	//}

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
	//pert_sequence->add_mover( pert_rt);

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
	//desn_mm->set_bb( true );
	desn_mm->set_chi( true );
	desn_mm->set_jump( 1, true );

	// create minimization mover
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min( new TaskAwareMinMover( desn_min, desn_tf ) );

	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/

	TR << "Main loop..." << std::endl;

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	//pose.dump_pdb("pre_main_loop.pdb");
	for ( Size k = 1; k <= Size( option[ oop_design::design_loop_num ].value() ); ++k ) {

		mc->reset(pose);

		// pert loop
		for( Size j = 1; j <= Size( option[ oop_design::pert_num ].value() ); ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
		}
		mc->recover_low( pose );

		// design
		TR << "DESIGN: " << k << std::endl;
		//kdrew: treating packer task as throw away object because it becomes invalid after design substitutions.
		PackerTaskOP task( TaskFactory::create_packer_task( pose ));
		//PackerTaskOP task = desn_tf->create_packer_task( pose ) ;

		utility::vector1<bool> allowed_aas(20, true);
		allowed_aas[aa_cys] = false;
		allowed_aas[aa_gly] = false;
		allowed_aas[aa_pro] = false;

		for (Size i=1; i<=pep_start-1; i++)
		{
			//TR << "  not designed" << std::endl;
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_from_command_line();
		}

		utility::vector1<Size> oop_designable_positions = option[ oop_design::oop_design_positions ].value();
		//kdrew: internal indexing for oop chain
		Size oop_pos = 0;
		// Set which residues can be designed
		for (Size i=pep_start; i<=pep_end; i++) {
			oop_pos++;
			TR << "oop postion " << oop_pos << std::endl;
			if ( oop_designable_positions.end() == find(oop_designable_positions.begin(), oop_designable_positions.end(), oop_pos ))
			{
				TR << "  not designed" << std::endl;
				task->nonconst_residue_task(i).restrict_to_repacking();
				task->nonconst_residue_task(i).initialize_from_command_line();
			}
			else {
				TR << "  designed" << std::endl;
				bool temp = allowed_aas[pose.residue(i).aa()];
				allowed_aas[pose.residue(i).aa()] = true;
				task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
				task->nonconst_residue_task(i).or_include_current(true);
				task->nonconst_residue_task(i).initialize_from_command_line();
				allowed_aas[pose.residue(i).aa()] = temp;
			}
		}

		//desn_tf->modify_task(pose, task);

		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover(score_fxn, task) );

		// create a sequence mover to hold pack rotamers and minimization movers
		moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
		desn_sequence->add_mover( desn_pr );
		desn_sequence->add_mover( desn_ta_min );
		desn_sequence->apply( pose );

		TR<< "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		TR<< "post mc->boltzmann" << std::endl;
		mc->show_state();

	}//dock_design for loop

	mc->recover_low( pose );

	curr_job->add_string_real_pair( "ENERGY_FINAL ", (*score_fxn)(pose) );

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
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
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
OopDesignMover::setup_pert_foldtree(
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
	Size dock_jump_pos_pro( core::pose::residue_center_of_mass( pose, pro_start, pro_end ) );
	Size dock_jump_pos_pep( core::pose::residue_center_of_mass( pose, pep_start, pep_end ) );

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
OopDesignMover::setup_filter_stats()
{
	/*********************************************************************************************************************
	Filter / Stats Setup
	*********************************************************************************************************************/

	// create and register sasa calculator
	pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
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


