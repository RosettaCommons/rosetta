// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/coupled_moves/CoupledMovesProtocol.cc
/// @brief Mover that implements the CoupledMovesProtocol
/// @author Noah <nollikai@gmail.com>, refactored by Steven Lewis smlewi@gmail.com
/// @author Anum Glasgow
/// @author Amanda Loshbaugh <aloshbau@gmail.com>

#include <protocols/coupled_moves/CoupledMovesProtocol.hh>
#include <protocols/coupled_moves/CoupledMovesProtocolCreator.hh>

#include <basic/Tracer.hh>

//RosettaScripts
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>

//
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinPackMover.hh>
#include <protocols/simple_moves/CoupledMover.hh>
#include <protocols/minimization_packing/MinPackMover.hh>
#include <protocols/minimization_packing/BoltzmannRotamerMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/jd2/util.hh>
#include <protocols/toolbox/IGEdgeReweighters.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelector.hh>
#include <core/pack/task/operation/ClashBasedRepackShell.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/xyz.functions.hh>
#include <ctime>
#include <iostream>
#include <string>
#include <sstream>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/coupled_moves.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.CoupledMovesProtocol" );

namespace protocols {
namespace coupled_moves {

/*options to deal with with the JD3 thing - leaving in a comment for ease of
notation later
OPT_1GRP_KEY(Integer, coupled_moves, ntrials)
OPT_1GRP_KEY(Integer, coupled_moves, number_ligands)
OPT_1GRP_KEY(Real, coupled_moves, mc_kt)
OPT_1GRP_KEY(Real, coupled_moves, boltzmann_kt)
OPT_1GRP_KEY(Real, coupled_moves, mm_bend_weight)
OPT_1GRP_KEY(Boolean, coupled_moves, trajectory)
OPT_1GRP_KEY(Boolean, coupled_moves, trajectory_gz)
OPT_1GRP_KEY(Integer, coupled_moves, trajectory_stride)
OPT_1GRP_KEY(String, coupled_moves, trajectory_file)
OPT_1GRP_KEY(String, coupled_moves, output_fasta)
OPT_1GRP_KEY(String, coupled_moves, output_stats)
OPT_1GRP_KEY(Boolean, coupled_moves, ligand_mode)
OPT_1GRP_KEY(Boolean, coupled_moves, initial_repack)
OPT_1GRP_KEY(Boolean, coupled_moves, min_pack)
OPT_1GRP_KEY(Boolean, coupled_moves, save_sequences)
OPT_1GRP_KEY(Boolean, coupled_moves, save_structures)
OPT_1GRP_KEY(Real, coupled_moves, ligand_prob)
OPT_1GRP_KEY(Boolean, coupled_moves, fix_backbone)
OPT_1GRP_KEY(Boolean, coupled_moves, uniform_backrub)
OPT_1GRP_KEY(Boolean, coupled_moves, bias_sampling)
OPT_1GRP_KEY(Boolean, coupled_moves, bump_check)
OPT_1GRP_KEY(Real, coupled_moves, ligand_weight)
OPT_1GRP_KEY(String, coupled_moves, output_prefix)

OPT(in::path::database);
OPT(in::ignore_unrecognized_res);
OPT(out::nstruct);
OPT(packing::resfile);
OPT(in::file::native);
OPT(constraints::cst_fa_weight);
OPT(constraints::cst_fa_file);
OPT(out::pdb_gz);
NEW_OPT(coupled_moves::ntrials, "number of Monte Carlo trials to run", 1000);
NEW_OPT(coupled_moves::number_ligands, "number of ligands in the pose", 1);
NEW_OPT(coupled_moves::mc_kt, "value of kT for Monte Carlo", 0.6);
NEW_OPT(coupled_moves::boltzmann_kt, "value of kT for Boltzmann weighted moves", 0.6);
NEW_OPT(coupled_moves::mm_bend_weight, "weight of mm_bend bond angle energy term", 1.0);
NEW_OPT(coupled_moves::trajectory, "record a trajectory", false);
NEW_OPT(coupled_moves::trajectory_gz, "gzip the trajectory", false);
NEW_OPT(coupled_moves::trajectory_stride, "write out a trajectory frame every N steps", 100);
NEW_OPT(coupled_moves::trajectory_file, "name of trajectory file", "traj.pdb");
NEW_OPT(coupled_moves::output_fasta, "name of FASTA output file", "sequences.fasta");
NEW_OPT(coupled_moves::output_stats, "name of stats output file", "sequences.stats");
NEW_OPT(coupled_moves::ligand_mode, "if true, model protein ligand interaction", false);
NEW_OPT(coupled_moves::initial_repack, "start simulation with repack and design step", true);
NEW_OPT(coupled_moves::min_pack, "use min_pack for initial repack and design step", false);
NEW_OPT(coupled_moves::save_sequences, "save all unique sequences", true);
NEW_OPT(coupled_moves::save_structures, "save structures for all unique sequences", false);
NEW_OPT(coupled_moves::ligand_prob, "probability of making a ligand move", 0.1);
NEW_OPT(coupled_moves::fix_backbone, "do not make any backbone moves", false);
NEW_OPT(coupled_moves::uniform_backrub, "select backrub rotation angle from uniform distribution", false);
NEW_OPT(coupled_moves::bias_sampling, "if true, bias rotamer selection based on energy", true);
NEW_OPT(coupled_moves::bump_check, "if true, use bump check in generating rotamers", true);
NEW_OPT(coupled_moves::ligand_weight, "weight for residue - ligand interactions", 1.0);
NEW_OPT(coupled_moves::output_prefix, "prefix for output files", "");
*/

CoupledMovesProtocol::CoupledMovesProtocol(): Mover(),
	score_fxn_(utility::pointer::make_shared< core::scoring::ScoreFunction >()),
	main_task_factory_(utility::pointer::make_shared< core::pack::task::TaskFactory >())
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	main_task_factory_->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
	main_task_factory_->push_back( utility::pointer::make_shared< IncludeCurrent >() );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( utility::pointer::make_shared< operation::ReadResfile >() );
	} else {
		operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
		main_task_factory_->push_back( rtrop );
	}
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	//main_task_factory_->push_back( new operation::PreserveCBeta );
	score_fxn_ = core::scoring::get_score_function();
	configure_score_fxn();
}

CoupledMovesProtocol::CoupledMovesProtocol( CoupledMovesProtocol const & /*cmp*/ ) = default;

core::Real CoupledMovesProtocol::compute_ligand_score_bonus(
	core::pose::PoseOP pose
) //utility::vector1<core::Size> ligand_resnums )
{
	core::scoring::EnergyMap weights = pose->energies().weights();
	core::scoring::EnergyGraph const & energy_graph( pose->energies().energy_graph() );
	core::scoring::EnergyMap ligand_two_body_energies;
	for ( core::Size i = 1; i <= pose->size(); i++ ) {
		for ( utility::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {
			const auto * edge( static_cast< const core::scoring::EnergyEdge *> (*iru) );
			core::Size const j( edge->get_first_node_ind() );
			core::Size const k( edge->get_second_node_ind() );
			for ( core::Size index = 1; index <= ligand_resnums_.size(); index++ ) {
				if ( j == ligand_resnums_[index] || k == ligand_resnums_[index] ) {
					ligand_two_body_energies += edge->fill_energy_map();
				}
			}
		}
	}
	return (ligand_two_body_energies.dot(weights) * (ligand_weight_ - 1.0));
}

void CoupledMovesProtocol::configure_score_fxn() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//SML moved from ctor without close inspection
	// set up the score function and add the bond angle energy term
	score_fxn_->set_weight(core::scoring::mm_bend, option[ OptionKeys::coupled_moves::mm_bend_weight ]);
	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	core::scoring::methods::EnergyMethodOptions energymethodoptions(score_fxn_->energy_method_options());
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(option[ OptionKeys::backrub::pivot_atoms ]);
	score_fxn_->set_energy_method_options(energymethodoptions);
	return;
}

void CoupledMovesProtocol::initialize_from_options() {
	TR << "CoupledMovesProtocol is examining the options system to set itself up.  If parse_my_tag and an active command line option disagree on a setting, the option is going to win." << std::endl;
	TR << "This applies only for ligand options at this time." << std::endl;
	initialize_ligand_from_options();
	return;
}

void CoupledMovesProtocol::initialize_ligand_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys::coupled_moves;
	if ( option[ligand_mode].user() ) { set_ligand_mode( option[ligand_mode].value() ); }
	if ( option[number_ligands].user() ) { set_number_ligands( option[number_ligands].value() ); }
	if ( option[ligand_weight].user() ) { set_ligand_weight( option[ligand_weight].value() ); }
	if ( option[ligand_prob].user() ) { set_ligand_prob( option[ligand_prob].value() ); }
	return;
}

void CoupledMovesProtocol::exclude_nonclashing_positions( core::pose::PoseOP pose ) {
	core::pack::task::operation::ClashBasedRepackShellOP cbrs( new core::pack::task::operation::ClashBasedRepackShell() );
	core::pack::task::PackerTaskOP task_for_cbss( main_task_factory_->create_task_and_apply_taskoperations( *pose ) );
	core::pack::task::residue_selector::ClashBasedShellSelectorOP rs( new core::pack::task::residue_selector::ClashBasedShellSelector( task_for_cbss, true ) ); // second argument sets focus_on_designable to true, so that ClashBasedShellSelector will find clashes for positions set to designable (NOT repackable) in task. This maintains the original behavior from (Ollikainen 2015).
	cbrs->selector( rs );
	main_task_factory_->push_back( cbrs );
	return;
} //exclude_nonclashing_positions()

void CoupledMovesProtocol::setup_move_and_design_positions( core::pack::task::PackerTaskCOP task, core::pose::PoseCOP pose_copy ) {
	move_positions_.clear();
	design_positions_.clear();
	for ( core::Size i = 1; i <= pose_copy->size(); ++i ) {
		//std::cout << "res " << i << "\tdesign " << task->design_residue(i) << "\tpack\t" << task->pack_residue(i) << std::endl;
		if ( ( task->design_residue(i) ) ) {
			move_positions_.push_back(i);
			design_positions_.push_back(i);
		} else if ( task->pack_residue(i) ) {
			move_positions_.push_back(i);
		}
	}
	bool const have_move(!move_positions_.empty());
	runtime_assert_string_msg(have_move, "move_positions empty in CoupledMovesProtocol. Probably you set everything to unrepackable.");
	bool const have_design(!design_positions_.empty());
	runtime_assert_string_msg(have_design, "design_positions empty in CoupledMovesProtocol. Probably you specified no design positions.");
	return;
} //setup_move_and_design_positions()

protocols::simple_moves::CoupledMoverOP CoupledMovesProtocol::setup_coupled_mover( core::pack::task::PackerTaskOP task, core::pose::PoseOP pose_copy, protocols::simple_moves::CoupledMoverOP coupled_mover ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( ligand_mode_ ) {
		coupled_mover = setup_coupled_mover_for_ligand( pose_copy, task, coupled_mover ) ;
	} else {
		coupled_mover = utility::pointer::make_shared< protocols::simple_moves::CoupledMover >(pose_copy, score_fxn_, task);
	}
	if ( backbone_mover_ == "backrub" ) {
		coupled_mover->set_backbone_mover( backbone_mover_ );
	} else if ( backbone_mover_ == "kic" ) {
		coupled_mover->set_backbone_mover( backbone_mover_ );
		coupled_mover->set_perturber( perturber_ );
	}
	coupled_mover->set_loop_size( loop_size_ );
	coupled_mover->set_fix_backbone( option[OptionKeys::coupled_moves::fix_backbone] );
	coupled_mover->set_bias_sampling( option[OptionKeys::coupled_moves::bias_sampling] );
	coupled_mover->set_temperature( option[OptionKeys::coupled_moves::boltzmann_kt] );
	coupled_mover->set_bump_check( option[OptionKeys::coupled_moves::bump_check] );
	coupled_mover->set_uniform_backrub( option[OptionKeys::coupled_moves::uniform_backrub] );
	return coupled_mover;
}

protocols::simple_moves::CoupledMoverOP CoupledMovesProtocol::setup_coupled_mover_for_ligand( core::pose::PoseOP pose_copy, core::pack::task::PackerTaskOP task, protocols::simple_moves::CoupledMoverOP coupled_mover ) {
	runtime_assert_string_msg((number_ligands_ > 0), "in CoupledMovesProtocol, ligand_mode was requested but with 0 number_ligands - don't do that");
	// ASSUMPTION: ligand(s) allowed sampling by CoupledMover (sampling includes changing conformation and translating in cartesian space) are the last residue(s) in the given PDB
	bool any_ligand(false);
	core::Size const nres(pose_copy->size());
	for ( core::Size res(1); res <= nres; ++res ) {
		if ( ( pose_copy->residue_type(res).is_ligand() ) &&  task->design_residue(res) ) {
			any_ligand = true;
			continue;
		}
		//if we've hit a designable ligand and later hit a nonligand:
		if ( any_ligand && (!pose_copy->residue_type(res).is_ligand()) ) {
			utility_exit_with_message("CoupledMovesProtocol requires that ligands allowed sampling by CoupledMover (sampling includes changing conformation and translating in cartesian space) be at the end of the Pose. ASSUMPTION: the ligand(s) are the last residue(s) in the given PDB.");
		}
	}
	for ( core::Size i(0); i < number_ligands_ ; ++i ) {
		core::Size const this_res(nres - i);
		ligand_resnums_.push_back( this_res );
		runtime_assert_string_msg(pose_copy->residue_type(this_res).is_ligand(), "in CoupledMovesProtocol, more ligands were requested than were present at the end of the Pose.  Note all designable ligands must be at the end of the input PDB.");
	}
	runtime_assert_string_msg((ligand_resnums_.size() > 0), "in CoupledMovesProtocol, ligand_resnums vector is empty - check that the number_ligands was set properly");
	coupled_mover = protocols::simple_moves::CoupledMoverOP(new protocols::simple_moves::CoupledMover(pose_copy, score_fxn_, task, ligand_resnums_[1]));
	coupled_mover->set_ligand_resnum( ligand_resnums_[1], pose_copy );
	coupled_mover->set_ligand_weight( ligand_weight_ );
	core::pack::task::IGEdgeReweighterOP reweight_ligand( new protocols::toolbox::IGLigandDesignEdgeUpweighter( ligand_weight_ ) );
	task->set_IGEdgeReweights()->add_reweighter( reweight_ligand );
	return coupled_mover;
} //setup_coupled_mover_for_ligand()

void CoupledMovesProtocol::setup_move_type( core::Real const move_prob ) {
	if ( ligand_mode_ && move_prob < ligand_prob_ ) {
		move_type_ = "LIGAND";
	} else {
		move_type_ = "RESIDUE";
	}
	return;
} // setup_move_type

void CoupledMovesProtocol::setup_resnum( core::Real const move_prob, protocols::simple_moves::CoupledMoverOP coupled_mover, core::pose::PoseOP pose_copy ) {
	core::Size const random = numeric::random::random_range(1, move_positions_.size());
	core::Size resnum = move_positions_[random]; //can't be const: might be reset by ligand mode
	if ( ligand_mode_ && move_prob < ligand_prob_ ) {
		resnum = ligand_resnums_[numeric::random::random_range(1, ligand_resnums_.size())];
		coupled_mover->set_ligand_resnum( resnum, pose_copy );
	}
	coupled_mover->set_resnum(resnum);
	return;
} // setup_resnum

void CoupledMovesProtocol::write_info_to_log( core::pose::PoseOP pose_copy ) {
	TR << "Design Positions: ";
	for ( core::Size i = 1; i <= design_positions_.size(); i++ ) {
		TR << pose_copy->pdb_info()->number(design_positions_[i]) << " ";
	}
	std::string initial_sequence = "";
	for ( core::Size index = 1; index <= design_positions_.size(); index++ ) {
		initial_sequence += pose_copy->residue(design_positions_[index]).name1();
	}
	TR << "Starting Sequence: " << initial_sequence << std::endl;
	TR << "Starting Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();
	TR << "Running " << /*basic::options::option[ basic::options::OptionKeys::coupled_moves::ntrials ]*/ ntrials_ << " trials..." << std::endl;
}

void CoupledMovesProtocol::apply( core::pose::Pose & pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	initialize_from_options();

	TR << "Initial Score:" << std::endl;
	score_fxn_->show(TR, pose);
	TR.flush();

	protocols::moves::MonteCarlo mc(pose, *score_fxn_, option[ OptionKeys::coupled_moves::mc_kt ]);
	protocols::viewer::add_monte_carlo_viewer(mc, "CoupledMoves", 600, 600);

	std::string output_tag(protocols::jd2::current_output_name());
	output_tag += option[OptionKeys::coupled_moves::output_prefix];

	// start with a fresh copy of the optimized pose
	core::pose::PoseOP pose_copy( new core::pose::Pose(pose) );
	// add constraints if supplied by via constraints::cst_file option
	core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose(*pose_copy);

	bool exclude_nonclashing_positions_option = option[ OptionKeys::coupled_moves::exclude_nonclashing_positions ];
	if ( exclude_nonclashing_positions_option ) {
		exclude_nonclashing_positions( pose_copy );
	}

	core::pack::task::PackerTaskOP task( main_task_factory_->create_task_and_apply_taskoperations( *pose_copy ) );
	setup_move_and_design_positions( task, pose_copy );
	protocols::simple_moves::CoupledMoverOP coupled_mover;
	coupled_mover = setup_coupled_mover( task, pose_copy, coupled_mover );
	if ( option[OptionKeys::coupled_moves::initial_repack] ) {
		if ( option[OptionKeys::coupled_moves::min_pack] ) {
			TR << "Performing initial minpack." << TR.Reset << std::endl;
			protocols::minimization_packing::MinPackMoverOP minpack(new protocols::minimization_packing::MinPackMover( score_fxn_, task ));
			minpack->apply(*pose_copy);
		} else {
			TR << "Performing initial pack." << TR.Reset << std::endl;
			protocols::minimization_packing::PackRotamersMoverOP pack( new protocols::minimization_packing::PackRotamersMover( score_fxn_, task, 1 ) );
			pack->apply(*pose_copy);
		}
	}

	// reset the Monte Carlo object
	mc.reset(*pose_copy);
	write_info_to_log( pose_copy );

	protocols::canonical_sampling::PDBTrajectoryRecorder trajectory;
	if ( option[ OptionKeys::coupled_moves::trajectory ] ) {
		trajectory.file_name(output_tag + "_traj.pdb" + (option[ OptionKeys::coupled_moves::trajectory_gz ] ? ".gz" : ""));
		trajectory.stride(option[ OptionKeys::coupled_moves::trajectory_stride ]);
		trajectory.reset(mc);
	}

	(*score_fxn_)(*pose_copy);
	core::Real current_score = pose_copy->energies().total_energy();

	if ( ligand_mode_ ) {
		core::Real const ligand_score_bonus = compute_ligand_score_bonus(pose_copy);
		current_score += ligand_score_bonus;
		mc.set_last_accepted_pose(*pose_copy, current_score);
	}

	//core::Size ntrials = option[ OptionKeys::coupled_moves::ntrials ];
	std::map<std::string,core::pose::Pose> unique_structures;
	unique_sequences_.clear();
	for ( core::Size i = 1; i <= ntrials_; ++i ) {
		//This assumes move_positions is not empty.  We have checked above, but let's check again.
		runtime_assert_string_msg(!move_positions_.empty(), "move_positions empty in CoupledMovesProtocol. Probably you specified no design positions or set everything to unrepackable.");
		core::Real const move_prob = numeric::random::uniform();
		setup_move_type( move_prob );
		setup_resnum( move_prob, coupled_mover, pose_copy );
		coupled_mover->apply(*pose_copy);

		// Score
		(*score_fxn_)(*pose_copy);
		current_score = pose_copy->energies().total_energy();
		core::Real ligand_score_bonus = 0.0;
		if ( ligand_mode_ ) {
			ligand_score_bonus = compute_ligand_score_bonus(pose_copy);
		}
		current_score += ligand_score_bonus;
		core::Real lowest_score = mc.lowest_score();

		// Monte Carlo
		bool accepted = mc.boltzmann(current_score, move_type_);
		if ( accepted ) {
			TR << "prev lowest_score: " << lowest_score << std::endl;
			if ( current_score < lowest_score ) {
				mc.set_lowest_score_pose(*pose_copy, current_score);
				TR << "new lowest_score: " << lowest_score << std::endl;
			}
			mc.set_last_accepted_pose(*pose_copy, current_score);
			std::string sequence = "";
			for ( core::Size index = 1; index <= design_positions_.size(); index++ ) {
				sequence += pose_copy->residue(design_positions_[index]).name1();
			}
			TR << i << " " << sequence << " " << current_score << std::endl;
			if ( option[OptionKeys::coupled_moves::save_sequences] ) {
				if ( unique_sequences_.find(sequence) == unique_sequences_.end() ) {
					TR << TR.Red << "Adding sequence " << sequence << TR.Reset << std::endl;
					unique_sequences_.insert(std::make_pair(sequence,current_score));
					unique_scores_.insert(std::make_pair(sequence,pose_copy->energies().total_energies()));
					if ( option[OptionKeys::coupled_moves::save_structures] ) {
						unique_structures.insert(std::make_pair(sequence,*pose_copy));
					}
				} else {
					if ( unique_sequences_[sequence] > current_score ) {
						unique_sequences_[sequence] = current_score;
						unique_scores_[sequence] = pose_copy->energies().total_energies();
						if ( option[OptionKeys::coupled_moves::save_structures] ) {
							unique_structures[sequence] = *pose_copy;
						}
					}
				}
			}
		} else {
			(*pose_copy) = mc.last_accepted_pose();
		}

		if ( option[ OptionKeys::coupled_moves::trajectory ] ) trajectory.update_after_boltzmann(mc);
	}
	mc.show_counters();

	// dump out the low score and last accepted poses

	//time_t current_time_time_t = time(0);
	//std::stringstream current_time_ss;
	//current_time_ss << current_time_time_t;
	//std::string current_time = current_time_ss.str();
	// in the event you wish to run multiple instances of coupled moves, for example as movers in a rosettascript, and you don't want the results to overwrite each other, you can label outputs with timestamps. just un-comment out the section above, and the dump_scored_pdb, fasta, and stats lines in the section below.

	TR << "Last Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();

	if ( option[out::pdb_gz] ) {
		pose_copy->dump_pdb(output_tag + "_last.pdb.gz");
	} else {
		pose_copy->dump_scored_pdb(output_tag + "_last.pdb", *score_fxn_);
		//pose_copy->dump_scored_pdb(output_tag + "_" + current_time + "_last.pdb", *score_fxn_);
	}

	*pose_copy = mc.lowest_score_pose();

	TR << "Low Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();

	if ( option[out::pdb_gz] ) {
		pose_copy->dump_pdb(output_tag + "_low.pdb.gz");
	} else {
		pose_copy->dump_scored_pdb(output_tag + "_low.pdb", *score_fxn_);
		//pose_copy->dump_scored_pdb(output_tag + "_" + current_time + "_low.pdb", *score_fxn_);
	}

	pose = mc.lowest_score_pose();

	if ( option[OptionKeys::coupled_moves::save_sequences] ) {
		std::ofstream out_fasta( (output_tag + ".fasta").c_str() );
		//std::ofstream out_fasta( (output_tag + "_" + current_time + ".fasta").c_str() );
		core::Size count = 1;
		for ( auto & unique_sequence : unique_sequences_ ) {
			out_fasta << ">Sequence" << count << " " << unique_sequence.second << std::endl;
			out_fasta << unique_sequence.first << std::endl;
			count++;
		}
		out_fasta.close();

		std::ofstream out_stats( (output_tag + ".stats").c_str() );
		//std::ofstream out_stats( (output_tag + "_" + current_time + ".stats").c_str() );
		count = 1;
		for ( auto & unique_sequence : unique_sequences_ ) {
			out_stats << "Sequence" << count << "\t" << unique_sequence.second << "\tsequence:\t" << unique_sequence.first << "\t" << unique_scores_[unique_sequence.first].weighted_string_of(score_fxn_->weights()) << std::endl;
			count++;
		}
		out_stats.close();

		if ( option[OptionKeys::coupled_moves::save_structures] ) {
			for ( auto & unique_structure : unique_structures ) {
				if ( option[out::pdb_gz] ) {
					unique_structure.second.dump_pdb(output_tag + "_" + unique_structure.first + "_low.pdb.gz");
				} else {
					unique_structure.second.dump_scored_pdb(output_tag + "_" + unique_structure.first + "_low.pdb", *score_fxn_);
				}
			}
		}
	}
}


// additional methods for RosettaScripts
protocols::moves::MoverOP CoupledMovesProtocolCreator::create_mover() const {
	return utility::pointer::make_shared< CoupledMovesProtocol >(); }
std::string CoupledMovesProtocolCreator::keyname() const {
	return CoupledMovesProtocol::mover_name(); }
protocols::moves::MoverOP CoupledMovesProtocol::clone() const {
	return utility::pointer::make_shared< CoupledMovesProtocol >( *this ); }
protocols::moves::MoverOP CoupledMovesProtocol::fresh_instance() const {
	return utility::pointer::make_shared< CoupledMovesProtocol >(); }
std::string CoupledMovesProtocol::get_name() const {
	return mover_name(); }
std::string CoupledMovesProtocol::mover_name() {
	return "CoupledMovesProtocol"; }
void CoupledMovesProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	CoupledMovesProtocol::provide_xml_schema( xsd ); }

void
CoupledMovesProtocol::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	main_task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data_map );
	score_fxn_ = protocols::rosetta_scripts::parse_score_function( tag, "score_fxn_", data_map );

	// Set the backbone mover
	if ( tag->getOption< bool >( "backbone_mover", false ) ) {
		std::string backbone_mover = option[ basic::options::OptionKeys::coupled_moves::backbone_mover ];

		if ( backbone_mover == "kic" ) {
			set_backbone_mover( backbone_mover );
			std::string kic_perturber = option[ basic::options::OptionKeys::coupled_moves::kic_perturber ];
			if ( kic_perturber == "fragment" ) {
				if ( ( option[ basic::options::OptionKeys::loops::frag_files ].user() ) && ( option[ basic::options::OptionKeys::loops::frag_sizes ].user() ) ) {
					utility::vector1<core::fragment::FragSetOP> frag_libs;
					protocols::loops::read_loop_fragments(frag_libs); // this function uses OptionKeys::loops::frag_sizes and OptionKeys::loops::frag_files to fill the frag_libs object, which is then used as an argument for the FragmentPerturber constructor.
					set_perturber( utility::pointer::make_shared< protocols::kinematic_closure::perturbers::FragmentPerturber >(frag_libs) );
				} else {
					std::stringstream message;
					message << "[ ERROR - fragments ] Must specify the -loops:frag_sizes and -loops:frag_files " << std::endl;
					message << "[ ERROR - fragments ] options in order to use the FragmentPerturber." << std::endl;
					throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
				}

			} else if ( kic_perturber == "walking" ) {
				set_perturber( utility::pointer::make_shared< protocols::kinematic_closure::perturbers::WalkingPerturber >( tag->getOption<core::Real>( "walking_perturber_magnitude", 2 ) ) );
			}
		} else if ( ( backbone_mover == "backrub" ) || ( backbone_mover == "" ) ) {
			set_backbone_mover( backbone_mover );
		} else {
			std::stringstream message;
			message << "[ ERROR - backbone_mover ] Specified -backbone_mover '" << backbone_mover << "' is not recognized. Try 'kic' or 'backrub' instead." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
		}
	} else {
		// If no backbone mover is specified in command line
		TR << TR.White << "[ WARNING - backbone_mover ] You did not specify -coupled_moves::backbone_mover option." << std::endl;
		TR << TR.White << "[ WARNING - backbone_mover ] Defaulting to '-coupled_moves::backbone_mover backrub'." << std::endl;
		TR << TR.White << "[ WARNING - backbone_mover ] Example usages:" << std::endl;
		TR << TR.White << "[ WARNING - backbone_mover ]     '-coupled_moves::backbone_mover kic'" << std::endl;
		TR << TR.White << "[ WARNING - backbone_mover ]     '-coupled_moves::backbone_mover backrub'" << std::endl;
		TR << TR.White << "[ WARNING - backbone_mover ] See CoupledMoves wiki or doxygen for more information on available backbone movers." << TR.Reset << std::endl;
	}

	set_loop_size( tag->getOption<core::Real>( "kic_loop_size", 9 ) );
	core::Size remainder;
	remainder = loop_size_ % 2;
	if ( remainder != 1 ) {
		std::stringstream message;
		message << "[ ERROR ] -kic_loop_size '" << loop_size_ << "' must be an odd number. See documentation." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
	}
	configure_score_fxn(); //SML this will add the extra bells and whistles to the scorefunction that the constructor used to do.  If the user actually sets those up in their scorefunction externally - this is a bad idea.  If the user has just passed REF15 as their scorefunction - this is necessary.  Anum approved at RosettaCON 2017.
}

void CoupledMovesProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute(
		"main_task_factory_", xs_string,
		"packer task");
	attlist + XMLSchemaAttribute(
		"score_fxn_", xs_string,
		"score function");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"backbone_mover", xsct_rosetta_bool,
		"make mover lisen to the coupled_moves related command line flags",
		std::to_string( false ));

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Small backbone and side chain movements",
		attlist );

}

} //coupled_moves
} //protocols
