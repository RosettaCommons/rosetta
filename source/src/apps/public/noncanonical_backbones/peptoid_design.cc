// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// based on code from mini/src/apps/public/scenarios/doug_dock_design_min_mod2_cal_cal.cc
//			and https://svn.rosettacommons.org/source/branches/releases/rosetta-3.1/manual/advanced/example_protocol.cc

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
		void setup_pert_foldtree( core::pose::Pose & pose);
		void setup_filter_stats();
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "PeptoidDesignMover"; }

};

typedef utility::pointer::owning_ptr< PeptoidDesignMover > PeptoidDesignMoverOP;
typedef utility::pointer::owning_ptr< PeptoidDesignMover const > PeptoidDesignMoverCOP;


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

		OD_mover->setup_filter_stats();

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
	score_fxn->set_weight(	atom_pair_constraint, 10 );
	score_fxn->set_weight(	angle_constraint, 10 );
	score_fxn->set_weight(	dihedral_constraint, 10 );

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
	desn_tf->push_back( new core::pack::task::operation::InitializeFromCommandline );

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
	simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( desn_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	TaskAwareMinMoverOP desn_ta_min = new TaskAwareMinMover( desn_min, desn_tf );

	// create a list of peptoid sidechains (this is inefficient)
	chemical::ResidueTypeCOPs const & rt_caps( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->residue_types() );
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
	for ( std::set< std::string >::const_iterator j(peptoid_name_set.begin()); j != peptoid_name_set.end(); ++j) {
		TR << *j << " ";
	}
	TR << std::endl;


	/*********************************************************************************************************************
	Main Loop
	**********************************************************************************************************************/

	TR << "Main loop..." << std::endl;

	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	//pose.dump_pdb("pre_main_loop.pdb");
	for ( Size k = 1; k <= Size( option[ peptoid_design::design_loop_num ].value() ); ++k ) {

		mc->reset(pose);

		// pert loop
		for( Size j = 1; j <= Size( option[ peptoid_design::pert_num ].value() ); ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
		}
		mc->recover_low( pose );

		// design
		TR << "DESIGN: " << k << std::endl;
		//kdrew: treating packer task as throw away object because it becomes invalid after design substitutions.
		PackerTaskOP task( TaskFactory::create_packer_task( pose ) );
		//PackerTaskOP task = desn_tf->create_packer_task( pose ) ;


		// set all residues not in chain 2 to repack
		for (Size i=1; i<=pep_start-1; i++) {
			//TR << "  not designed" << std::endl;
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_from_command_line();
		}

		utility::vector1<Size> peptoid_designable_positions = option[ peptoid_design::peptoid_design_positions ].value();
		//kdrew: internal indexing for peptoid chain
		Size peptoid_pos = 0;
		// Set which residues can be designed
		for (Size i=pep_start; i<=pep_end; i++) {
			peptoid_pos++;
			TR << "peptoid postion " << peptoid_pos << std::endl;
			if ( peptoid_designable_positions.end() == find(peptoid_designable_positions.begin(), peptoid_designable_positions.end(), peptoid_pos )) {
				TR << "  not designed" << std::endl;
				task->nonconst_residue_task(i).restrict_to_repacking();
				task->nonconst_residue_task(i).initialize_from_command_line();
			}	else {
				TR << "  designed" << std::endl;
				for ( std::set< std::string >::const_iterator j(peptoid_name_set.begin()); j != peptoid_name_set.end(); ++j) {
					task->nonconst_residue_task(i).allow_noncanonical_aa( *j );
				}
				task->nonconst_residue_task(i).or_include_current(true);
				task->nonconst_residue_task(i).initialize_from_command_line();
			}
		}

		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP desn_pr( new simple_moves::PackRotamersMover(score_fxn, task) );

		// create a cyclization mover to run before design
		simple_moves::CyclizationMoverOP desn_cyclization( new CyclizationMover( 2, false, true, 3, score_fxn, desn_mm ) );

		// create a sequence mover to hold pack rotamers and minimization movers
		moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
		if ( option[ peptoid_design::cyclic ].value() == true ) {
			desn_sequence->add_mover( desn_cyclization );
		}
		//desn_sequence->add_mover( pmm );
		desn_sequence->add_mover( desn_pr );
		//desn_sequence->add_mover( pmm );
		desn_sequence->add_mover( desn_ta_min );
		//desn_sequence->add_mover( pmm );
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
PeptoidDesignMover::setup_pert_foldtree(
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
PeptoidDesignMover::setup_filter_stats()
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


