// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Utilities used by various NCBB based design and dock-design applications.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/jd2/util.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
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
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


// tracer - used to replace cout
static basic::Tracer TR( "NCBB" );

namespace protocols {
namespace ncbb {

Size
give_dihedral_index(
	Size n,
	utility::vector1< char > uniqs,
	std::string dihedral_pattern,
	std::string alpha_beta_pattern
) {
	// for all characters from uniqs that are earlier than this residue's
	// add 2 or 3 depending on if those uniqs are for alphas or betas

	Size index = 1;

	for ( Size i = 1; i < n; ++i ) {
		for ( Size j = 0; j < dihedral_pattern.length(); ++j ) {
			if ( dihedral_pattern[j] == uniqs[i] ) {
				if ( alpha_beta_pattern[j] == 'A' || alpha_beta_pattern[j]=='P' ) {
					// it's an alpha
					index += 2;
				} else if ( alpha_beta_pattern[j] == 'B' ) {
					index += 3;
				}
				j = dihedral_pattern.length();
			}
		}
	}
	return index;
}

Size
get_number_dihedrals (
	utility::vector1< char > uniqs,
	std::string const & dihedral_pattern,
	std::string const & alpha_beta_pattern
) {
	Size number = 0;
	for ( Size i = 1; i <= uniqs.size(); ++i ) {
		// look through dihedral pattern until you find uniqs
		for ( Size j = 0; j < dihedral_pattern.length(); ++j ) {
			if ( dihedral_pattern[j] == uniqs[i] ) {
				if ( alpha_beta_pattern[j] == 'A' || alpha_beta_pattern[j] == 'P' ) {
					// it's an alpha or a peptoid
					// TODO: figure out an appropriate way to decide if you want to do omega sampling
					// for peptoids; I don't suggest this yet--at worst you can minimize it later?
					number += 2;
				} else if ( alpha_beta_pattern[j] == 'B' ) {
					number += 3;
				}
				j = dihedral_pattern.length();
			}
		}
	}
	return number;
}

void
count_uniq_char( std::string pattern, Size & num, utility::vector1<char> & uniqs ) {
	num = 1;
	uniqs.push_back( pattern[ 0 ] );
	for ( Size i = 1; i < pattern.length(); ++i ) {
		// assume i unique
		Size increment = 1;
		for ( Size j = 0; j < i; ++j ) {
			if ( pattern[j]  == pattern[ i ] ) {
				// repeat
				increment = 0;
			} else {
				// new!
				increment *= 1;
			}
		}
		if ( increment == 1 ) uniqs.push_back( pattern[ i ] );
		num += increment;
	}
}

void
ncbb_design_main_loop( Size loop_num, Size pert_num, core::pose::Pose pose, TrialMoverOP pert_trial, utility::vector1<Size> designable_positions, Size pep_start, Size pep_end, TaskAwareMinMoverOP desn_ta_min, core::scoring::ScoreFunctionOP score_fxn, MonteCarloOP mc ) {

	for ( Size k = 1; k <= loop_num; ++k ) {

		mc->reset(pose);

		// pert loop
		for ( Size j = 1; j <= pert_num; ++j ) {
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

		for ( Size i=1; i<=pep_start-1; i++ ) {
			//TR << "  not designed" << std::endl;
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_from_command_line();
		}

		//kdrew: internal indexing for chain
		Size pos = 0;
		// Set which residues can be designed
		for ( Size i=pep_start; i<=pep_end; i++ ) {
			pos++;
			TR << "position " << pos << std::endl;
			if ( designable_positions.end() == find(designable_positions.begin(), designable_positions.end(), pos ) ) {
				TR << "  not designed" << std::endl;
				task->nonconst_residue_task(i).restrict_to_repacking();
				task->nonconst_residue_task(i).initialize_from_command_line();
			} else {
				TR << "  designed" << std::endl;
				bool temp = allowed_aas[pose.residue(i).aa()];
				allowed_aas[pose.residue(i).aa()] = true;
				task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
				task->nonconst_residue_task(i).or_include_current(true);
				task->nonconst_residue_task(i).initialize_from_command_line();
				allowed_aas[pose.residue(i).aa()] = temp;
			}
		}

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
}

void
final_design_min( core::pose::Pose & pose, ScoreFunctionOP score_fxn_, core::pack::task::TaskFactoryOP desn_tf ) {
	// get packer task from task factory
	PackerTaskOP final_desn_pt( desn_tf->create_task_and_apply_taskoperations( pose ) );

	// add extra chi and extra chi cut off to pt
	for ( Size i = 1; i <= pose.size(); ++i ) {
		final_desn_pt->nonconst_residue_task( i ).or_ex1( true );
		final_desn_pt->nonconst_residue_task( i ).or_ex2( true );
		final_desn_pt->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
	}

	// create a pack rotamers mover for the final design
	simple_moves::PackRotamersMoverOP final_desn_pr( new simple_moves::PackRotamersMover(score_fxn_, final_desn_pt, 10 ) );
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
	simple_moves::MinMoverOP final_min( new simple_moves::MinMover( final_min_mm, score_fxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );
	// final min (okay to use ta min here)
	final_min->apply( pose );
}

void
calculate_statistics( core::pose::Pose const & pose, core::scoring::ScoreFunctionOP score_fxn ) {

	// create  MetricValues
	basic::MetricValue<  Real  > mv_sasa_complex;
	basic::MetricValue<  Real  > mv_sasa_seperated;
	basic::MetricValue< utility::vector1< core::Size > > mv_unsat_res_complex;
	basic::MetricValue< utility::vector1< core::Size > > mv_unsat_res_seperated;
	basic::MetricValue<  Real  > mv_pack_complex;
	basic::MetricValue<  Real  > mv_pack_seperated;

	basic::MetricValue<  Real  > mv_repack_sasa_seperated;
	basic::MetricValue< utility::vector1< core::Size > > mv_repack_unsat_res_seperated;
	basic::MetricValue<  Real  > mv_repack_pack_seperated;
	Real  repack_energy_seperated;
	Real  repack_hbond_ener_sum_seperated;

	Real  energy_complex;
	Real  energy_seperated;
	Real  hbond_ener_sum_complex;
	Real  hbond_ener_sum_seperated;

	// calc energy
	//energy_complex = (*score_fxn)(pose);

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
	tf->push_back( operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
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
	simple_moves::MinMoverOP separate_min( new simple_moves::MinMover( separate_min_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );
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
	protocols::jd2::add_string_real_pair_to_current_job( "ENERGY_COMPLEX:\t\t", energy_complex );
	protocols::jd2::add_string_real_pair_to_current_job( "ENERGY_SEPERATE:\t\t", energy_seperated );
	protocols::jd2::add_string_real_pair_to_current_job( "ENERGY_DIFF:\t\t", energy_complex - energy_seperated );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_ENERGY_SEPERATE:\t\t", repack_energy_seperated );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_ENERGY_DIFF:\t\t", energy_complex - repack_energy_seperated );

	protocols::jd2::add_string_real_pair_to_current_job( "SASA_COMPLEX:\t\t", mv_sasa_complex.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "SASA_SEPERATE:\t\t", mv_sasa_seperated.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "SASA_DIFF:\t\t", mv_sasa_complex.value() - mv_sasa_seperated.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_SASA_SEPERATE:\t\t", mv_repack_sasa_seperated.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_SASA_DIFF:\t\t", mv_sasa_complex.value() - mv_repack_sasa_seperated.value() );

	protocols::jd2::add_string_real_pair_to_current_job( "HB_ENER_COMPLEX:\t\t", hbond_ener_sum_complex );
	protocols::jd2::add_string_real_pair_to_current_job( "HB_ENER_SEPERATE:\t\t", hbond_ener_sum_seperated );
	protocols::jd2::add_string_real_pair_to_current_job( "HB_ENER_DIFF:\t\t", hbond_ener_sum_complex - hbond_ener_sum_seperated );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_HB_ENER_SEPERATE:\t\t", repack_hbond_ener_sum_seperated );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_HB_ENER_DIFF:\t\t", hbond_ener_sum_complex - repack_hbond_ener_sum_seperated );

	protocols::jd2::add_string_real_pair_to_current_job( "PACK_COMPLEX:\t\t", mv_pack_complex.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "PACK_SEPERATE:\t\t", mv_pack_seperated.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "PACK_DIFF:\t\t", mv_pack_complex.value() - mv_pack_seperated.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_PACK_SEPERATE:\t\t", mv_repack_pack_seperated.value() );
	protocols::jd2::add_string_real_pair_to_current_job( "REPACK_PACK_DIFF:\t\t", mv_pack_complex.value() - mv_repack_pack_seperated.value() );
}

// this only works for two chains and assumes the protein is first and the peptide is second
// inspired by protocols/docking/DockingProtocol.cc
void
setup_pert_foldtree(
	core::pose::Pose & pose
) {
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
	debug_assert( f.check_fold_tree() );

	std::cout << "AFTER: " << f << std::endl;

	pose.fold_tree( f );
}

void
setup_filter_stats() {
	/*********************************************************************************************************************
	Filter / Stats Setup
	*********************************************************************************************************************/

	// create and register sasa calculator
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		Warning() << "In setup_filter_stats, calculator " << "sasa"
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
		pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}

	// create and register hb calculator
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "num_hbonds" ) ) {
		Warning() << "In setup_filters_stats, calculator " << "num_hbonds"
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator( new pose_metric_calculators::NumberHBondsCalculator() );
		pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );
	}

	// create and register unsat calculator
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "unsat" ) ) {
		Warning() << "In setup_filters_stats, calculator " << "unsat"
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		pose::metrics::PoseMetricCalculatorOP unsat_calculator( new pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds") ) ;
		pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );
	}

	// create and register packstat calculator
	pose::metrics::PoseMetricCalculatorOP pack_calcculator( new pose_metric_calculators::PackstatCalculator() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "pack", pack_calcculator );

}

void
init_common_options( utility::tag::TagCOP tag, basic::datacache::DataMap &data, ScoreFunctionOP score_fxn_, Real & mc_temp_, Real & pert_mc_temp_, Real & pert_dock_rot_mag_, Real & pert_dock_trans_mag_, Real & pert_pep_small_temp_, Real & pert_pep_small_H_, Real & pert_pep_small_L_, Real & pert_pep_small_E_, Real & pert_pep_shear_temp_, Real & pert_pep_shear_H_, Real & pert_pep_shear_L_, Real & pert_pep_shear_E_, Size & pert_pep_num_rep_, Size & pert_num_, Size & dock_design_loop_num_, bool & no_design_, bool & final_design_min_, bool & use_soft_rep_, bool & mc_initial_pose_, bool & pymol_, bool & keep_history_ ) {

	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_key( tag->getOption< std::string >( "scorefxn" ) );
		if ( ! data.has( "scorefxns", scorefxn_key ) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
		score_fxn_ = data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_key );
	}

	if ( tag->hasOption( "mc_temp" ) ) {
		mc_temp_ = tag->getOption< Real >( "mc_temp", mc_temp_ );
	} else {
		mc_temp_ = 1.0;
	}

	if ( tag->hasOption( "pert_mc_temp" ) ) {
		pert_mc_temp_ = tag->getOption< Real >( "pert_mc_temp", pert_mc_temp_ );
	} else {
		pert_mc_temp_ = 0.8;
	}

	if ( tag->hasOption( "pert_dock_rot_mag" ) ) {
		pert_dock_rot_mag_ = tag->getOption< Real >( "pert_dock_rot_mag", pert_dock_rot_mag_ );
	} else {
		pert_dock_rot_mag_ = 1.0;
	}

	if ( tag->hasOption( "pert_dock_trans_mag" ) ) {
		pert_dock_trans_mag_ = tag->getOption< Real >( "pert_dock_trans_mag", pert_dock_trans_mag_ );
	} else {
		pert_dock_trans_mag_ = 0.5;
	}

	if ( tag->hasOption( "pert_pep_small_temp" ) ) {
		pert_pep_small_temp_ = tag->getOption< Real >( "pert_pep_small_temp", pert_pep_small_temp_ );
	} else {
		pert_pep_small_temp_ = 0.8;
	}

	if ( tag->hasOption( "pert_pep_small_H" ) ) {
		pert_pep_small_H_ = tag->getOption< Real >( "pert_pep_small_H", pert_pep_small_H_ );
	} else {
		pert_pep_small_H_ = 2.0;
	}

	if ( tag->hasOption( "pert_pep_small_L" ) ) {
		pert_pep_small_L_ = tag->getOption< Real >( "pert_pep_small_L", pert_pep_small_L_ );
	} else {
		pert_pep_small_L_ = 2.0;
	}

	if ( tag->hasOption( "pert_pep_small_E" ) ) {
		pert_pep_small_E_ = tag->getOption< Real >( "pert_pep_small_E", pert_pep_small_E_ );
	} else {
		pert_pep_small_E_ = 2.0;
	}

	if ( tag->hasOption( "pert_pep_shear_temp" ) ) {
		pert_pep_shear_temp_ = tag->getOption< Real >( "pert_pep_shear_temp", pert_pep_shear_temp_ );
	} else {
		pert_pep_shear_temp_ = 0.8;
	}

	if ( tag->hasOption( "pert_pep_shear_H" ) ) {
		pert_pep_shear_H_ = tag->getOption< Real >( "pert_pep_shear_H", pert_pep_shear_H_ );
	} else {
		pert_pep_shear_H_ = 2.0;
	}

	if ( tag->hasOption( "pert_pep_shear_L" ) ) {
		pert_pep_shear_L_ = tag->getOption< Real >( "pert_pep_shear_L", pert_pep_shear_L_ );
	} else {
		pert_pep_shear_L_ = 2.0;
	}

	if ( tag->hasOption( "pert_pep_shear_E" ) ) {
		pert_pep_shear_E_ = tag->getOption< Real >( "pert_pep_shear_E", pert_pep_shear_E_ );
	} else {
		pert_pep_shear_E_ = 2.0;
	}

	if ( tag->hasOption( "pert_pep_num_rep" ) ) {
		pert_pep_num_rep_ = tag->getOption< Size >("pert_pep_num_rep", pert_pep_num_rep_);
	} else {
		pert_pep_num_rep_ = 100;
	}

	if ( tag->hasOption( "pert_num" ) ) {
		pert_num_ = tag->getOption< Size >( "pert_num", pert_num_ );
	} else {
		pert_num_ = 10;
	}

	if ( tag->hasOption( "dock_design_loop_num" ) ) {
		dock_design_loop_num_ = tag->getOption< Size >( "dock_design_loop_num", dock_design_loop_num_ );
	} else {
		dock_design_loop_num_ = 10;
	}

	if ( tag->hasOption( "no_design" ) ) {
		no_design_ = tag->getOption< bool >( "no_design", no_design_ );
	} else {
		no_design_ = false;
	}

	if ( tag->hasOption( "final_design_min" ) ) {
		final_design_min_ = tag->getOption< bool >( "final_design_min", final_design_min_ );
	} else {
		final_design_min_ = true;
	}

	if ( tag->hasOption( "use_soft_rep" ) ) {
		use_soft_rep_ = tag->getOption< bool >("use_soft_rep", use_soft_rep_);
	} else {
		use_soft_rep_ = false;
	}

	if ( tag->hasOption( "mc_initial_pose" ) ) {
		mc_initial_pose_ = tag->getOption< bool >( "mc_initial_pose", mc_initial_pose_ );
	} else {
		mc_initial_pose_ = false;
	}

	if ( tag->hasOption( "pymol" ) ) {
		pymol_ = tag->getOption< bool >( "pymol", pymol_ );
	} else {
		pymol_ = false;
	}

	if ( tag->hasOption( "keep_history" ) ) {
		keep_history_ = tag->getOption< bool >( "keep_history", keep_history_ );
	} else {
		keep_history_ = false;
	}
}

using namespace utility::tag;
void add_attributes_for_common_options( AttributeList & attlist ) {
	attlist + XMLSchemaAttribute::attribute_w_default( "mc_temp", xsct_real, "Temperature of Monte Carlo simulation", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_mc_temp", xsct_real, "Temperature of Monte Carlo simulation for the perturbation stage", "0.8" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_dock_rot_mag", xsct_real, "Size of rotations for docking", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_dock_trans_mag", xsct_real, "Size of translations for docking", "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_small_temp", xsct_real, "Temperature of Monte Carlo simulation for peptide small moves", "0.8" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_small_H", xsct_real, "Size of small moves on helices", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_small_L", xsct_real, "Size of small moves on loops", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_small_E", xsct_real, "Size of small moves on strands", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_shear_temp", xsct_real, "Temperature of Monte Carlo simulation for peptide shear moves", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_shear_H", xsct_real, "Size of small moves on helices", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_shear_L", xsct_real, "Size of small moves on loops", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_shear_E", xsct_real, "Size of small moves on strands", "2.0" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_pep_num_rep", xsct_non_negative_integer, "Number of peptide moves per cycle", "100" )
		+ XMLSchemaAttribute::attribute_w_default( "pert_num", xsct_non_negative_integer, "Rounds of perturbation", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "no_design", xsct_rosetta_bool, "Turn off design", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "final_design_min", xsct_rosetta_bool, "Minimize after design", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_soft_rep", xsct_rosetta_bool, "Use a soft-repulsive scoring function", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "mc_initial_pose", xsct_rosetta_bool, "Run a simple MC simulation on the initial pose", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "pymol", xsct_rosetta_bool, "Use the pymol mover", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "keep_history", xsct_rosetta_bool, "Keep the history of the simulation", "false" );
}

}
}
