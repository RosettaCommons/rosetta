// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file backrub.cc
/// @brief run backrub Monte Carlo
/// @author Colin A. Smith (colin.smith@ucsf.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) (Move to Protocols from apps, add Movemap support, some refactoring)
/// @details
/// Currently a work in progress. The goal is to match the features of rosetta++ -backrub_mc


// Project Headers
#include <protocols/backrub/BackrubProtocol.hh>
#include <protocols/backrub/BackrubProtocolCreator.hh>
#include <protocols/backrub/util.hh>

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/viewer/viewers.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/keys/Key3Vector.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.backrub.BackrubProtocol");

namespace protocols {
namespace backrub {


BackrubProtocol::BackrubProtocol(): Mover(),
	scorefxn_(nullptr),
	main_task_factory_(nullptr),
	backrubmover_(protocols::backrub::BackrubMoverOP( new protocols::backrub::BackrubMover() )),
	smallmover_(protocols::simple_moves::SmallMoverOP(new protocols::simple_moves::SmallMover())),
	sidechainmover_(protocols::simple_moves::sidechain_moves::SidechainMoverOP(new protocols::simple_moves::sidechain_moves::SidechainMover())),
	packrotamersmover_(protocols::simple_moves::PackRotamersMoverOP(new protocols::simple_moves::PackRotamersMover())),
	movemap_smallmover_(nullptr)
{
	read_cmd_line_options();
}

BackrubProtocol::BackrubProtocol(BackrubProtocol const & bp): Mover(bp),
	scorefxn_(bp.scorefxn_),
	main_task_factory_(bp.main_task_factory_),
	backrubmover_(bp.backrubmover_),
	smallmover_(bp.smallmover_),
	sidechainmover_(bp.sidechainmover_),
	packrotamersmover_(bp.packrotamersmover_),
	ntrials_(bp.ntrials_),
	movemap_smallmover_(bp.movemap_smallmover_),
	pivot_residues_(bp.pivot_residues_),
	pivot_atoms_(bp.pivot_atoms_),
	min_atoms_(bp.min_atoms_),
	max_atoms_(bp.max_atoms_),
	sm_prob_(bp.sm_prob_),
	sc_prob_(bp.sc_prob_),
	sc_prob_uniform_(bp.sc_prob_uniform_),
	sc_prob_withinrot_(bp.sc_prob_withinrot_),
	mc_kt_(bp.mc_kt_),
	initial_pack_(bp.initial_pack_)
{}

BackrubProtocol::~BackrubProtocol()= default;

void
BackrubProtocol::read_cmd_line_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	pivot_residues_.clear();
	pivot_atoms_.clear();

	if ( option[ OptionKeys::backrub::pivot_residues].user() ) {
		for ( core::Size i = 1; i <= option[ OptionKeys::backrub::pivot_residues ].size(); ++i ) {
			if ( option[ OptionKeys::backrub::pivot_residues ][i] >= 1 ) pivot_residues_.push_back(option[ OptionKeys::backrub::pivot_residues ][i]);
		}
	}

	pivot_atoms_ = option[ OptionKeys::backrub::pivot_atoms ];
	min_atoms_ = option[ OptionKeys::backrub::min_atoms ];
	max_atoms_ = option[ OptionKeys::backrub::max_atoms ];

	initial_pack_ = option[ OptionKeys::backrub::initial_pack ];
	sm_prob_ = option[ OptionKeys::backrub::sm_prob ];
	sc_prob_ = option[ OptionKeys::backrub::sc_prob ];

	sc_prob_uniform_ = option[ OptionKeys::backrub::sc_prob_uniform ];
	sc_prob_withinrot_ = option[ OptionKeys::backrub::sc_prob_withinrot ];
	mc_kt_ = option[ OptionKeys::backrub::mc_kt];

	ntrials_ = option[ OptionKeys::backrub::ntrials ];

}

void
BackrubProtocol::set_pivot_atoms(utility::vector1<std::string> pivot_atoms){
	pivot_atoms_ = pivot_atoms;
}

void
BackrubProtocol::set_pivot_residues(utility::vector1<core::Size> pivot_residues) {
	pivot_residues_ = pivot_residues;
}

void
BackrubProtocol::set_movemap(core::kinematics::MoveMapCOP movemap){
	pivot_residues_ = get_pivot_residues_from_movemap(movemap);
	set_movemap_smallmover(movemap);

}

void
BackrubProtocol::set_movemap_smallmover(core::kinematics::MoveMapCOP movemap){
	movemap_smallmover_ = movemap;
}

void
BackrubProtocol::set_taskfactory(core::pack::task::TaskFactoryCOP tf){
	main_task_factory_ = tf->clone();
}

void
BackrubProtocol::set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn->clone();
}

protocols::moves::MoverOP
BackrubProtocol::clone() const {
	return protocols::moves::MoverOP( new BackrubProtocol( *this ) );
}

protocols::moves::MoverOP
BackrubProtocol::fresh_instance() const {
	return protocols::moves::MoverOP( new BackrubProtocol );
}

std::string
BackrubProtocol::get_name() const {
	return "BackrubProtocol";
}

void
BackrubProtocol::set_backrub_mover(protocols::backrub::BackrubMoverOP backrub_mover){
	backrubmover_ = backrub_mover;
}

void
BackrubProtocol::write_database() {
	backrubmover_->branchopt().write_database();
}

void
BackrubProtocol::finalize_setup(core::pose::Pose & pose){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	if ( ! main_task_factory_ ) {
		main_task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory() );
		main_task_factory_->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			main_task_factory_->push_back( TaskOperationCOP( new operation::ReadResfile ) );
		} else {
			operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
			main_task_factory_->push_back( rtrop );
		}
	}

	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory_->push_back( TaskOperationCOP( new operation::PreserveCBeta ) );

	// set up the score function and add the bond angle energy term

	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}
	scorefxn_->set_weight_if_zero(core::scoring::mm_bend, option[ OptionKeys::backrub::mm_bend_weight ]);
	core::scoring::methods::EnergyMethodOptions energymethodoptions(scorefxn_->energy_method_options());
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(pivot_atoms_);
	scorefxn_->set_energy_method_options(energymethodoptions);
	if ( pose.is_centroid() ) {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	} else {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	}

	// read known and unknown optimization parameters from the database
	backrubmover_->branchopt().read_database();
	// tell the branch angle optimizer about the score function MMBondAngleResidueTypeParamSet, if any
	if ( energymethodoptions.bond_angle_residue_type_param_set() ) {
		backrubmover_->branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set());
	}

	// set up the SmallMover
	if ( sm_prob_ > 0 ) {
		if ( ! movemap_smallmover_ ) {
			core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap );
			mm->init_from_file(option[ OptionKeys::in::file::movemap ]);
			set_movemap_smallmover(mm);
		}
		smallmover_->nmoves(1);
		smallmover_->movemap(movemap_smallmover_->clone());
	}

	// set up the SidechainMover
	sidechainmover_->set_task_factory(main_task_factory_);
	sidechainmover_->set_prob_uniform(sc_prob_uniform_ );
	sidechainmover_->set_prob_withinrot(sc_prob_withinrot_ );

	// set up the PackRotamersMover
	packrotamersmover_->task_factory(main_task_factory_);
	packrotamersmover_->score_function(scorefxn_);

	if ( pivot_residues_.size() > 0 ) {
		backrubmover_->set_pivot_residues( pivot_residues_ );
	}

	backrubmover_->clear_segments();
	backrubmover_->set_pivot_atoms( pivot_atoms_ );
	backrubmover_->set_min_atoms( min_atoms_ );
	backrubmover_->set_max_atoms( max_atoms_ );
	backrubmover_->set_input_pose( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
	backrubmover_->add_mainchain_segments();

}

//void
//BackrubProtocol::parse_my_tag(
//  TagCOP tag,
//  basic::datacache::DataMap&,
//  const Filters_map&,
//  const Movers_map&,
//  const Pose&) {
//
//}

void
BackrubProtocol::apply( core::pose::Pose& pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	bool custom_fold_tree = false;
	std::string input_tag = protocols::jd2::JobDistributor::get_instance()->current_job()->inner_job()->input_tag();
	if ( pose.is_centroid() ) {
		TR.Warning << "*** This is untested with centroid mode! ***" << std::endl;
		//core::import_pose::centroid_pose_from_pdb( *input_pose, input_jobs[jobnum]->input_tag() , core::import_pose::PDB_file);
		custom_fold_tree = read_fold_tree_from_file( pose, input_tag );
		core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);
	} else {
		//core::import_pose::pose_from_file( *input_pose, input_jobs[jobnum]->input_tag() , core::import_pose::PDB_file);
		custom_fold_tree = read_fold_tree_from_file( pose, input_tag );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);
	}

	finalize_setup(pose);
	TR << "Score After PDB Load:" << std::endl;
	scorefxn_->show(TR, pose);
	TR.flush();

	backrubmover_->optimize_branch_angles(pose);
	//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postoptbranch.pdb");
	sidechainmover_->idealize_sidechains(pose);
	//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postidealizesc.pdb");
	if ( sc_prob_ && sidechainmover_->packed_residues().size() == 0 ) {
		sc_prob_ = 0;
		TR << "Warning: No side chains to move. Not using SidechainMover." << std::endl;
	}

	TR << "Score After Branch Angle Optimization/Side Chain Idealization:" << std::endl;
	scorefxn_->show(TR, pose);
	TR.flush();

	// check to see if we need to force a repack

	core::pack::task::PackerTaskOP temp_task(main_task_factory_->create_task_and_apply_taskoperations(pose));
	utility::vector1<core::Size> incompatible_positions(positions_incompatible_with_task(pose, *temp_task));
	if ( incompatible_positions.size() ) {
		TR << "Starting ResidueType not allowed in resfile at position(s):";
		for ( core::Size i = 1; i <= incompatible_positions.size(); i++ ) {
			TR << " " << incompatible_positions[i];
		}
		TR << std::endl;
		if ( !initial_pack_ ) {
			initial_pack_ = true;
			TR << "Forcing initial pack" << std::endl;
		}
	}

	protocols::moves::MonteCarlo mc(pose, *scorefxn_, mc_kt_ );

	// create viewer windows if OpenGL is enabled
	protocols::viewer::add_monte_carlo_viewer(mc, "Backrub", 600, 600);

	////////////////// Material above is scheduled for constructor //////////////////


	std::string output_tag(protocols::jd2::current_output_name());

	// start with a fresh copy of the optimized pose
	core::pose::PoseOP pose_copy( new core::pose::Pose(pose) );

	// repack/redesign at the beginning if specified/necessary
	if ( initial_pack_ ) {
		packrotamersmover_->apply(*pose_copy);
		//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postpack.pdb");

		// if a minimization movemap was specified, go through a series of minimizations
		if ( option[ OptionKeys::backrub::minimize_movemap ].user() ) {

			// setup the MoveMaps
			core::kinematics::MoveMapOP minimize_movemap( new core::kinematics::MoveMap );
			minimize_movemap->init_from_file(option[ OptionKeys::backrub::minimize_movemap ]);
			core::kinematics::MoveMapOP minimize_movemap_progressive( new core::kinematics::MoveMap );

			// setup the MinMover
			protocols::simple_moves::MinMover minmover;
			minmover.score_function(scorefxn_);
			minmover.min_type("dfpmin");

			// first minimize just the side chains
			for ( auto iter = minimize_movemap->movemap_torsion_id_begin();
					iter != minimize_movemap->movemap_torsion_id_end(); ++iter ) {
				if ( iter->first.second == core::id::CHI ) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminchi.pdb");

			// next minimize the side chains and backbone
			for ( auto iter = minimize_movemap->movemap_torsion_id_begin();
					iter != minimize_movemap->movemap_torsion_id_end(); ++iter ) {
				if ( iter->first.second == core::id::BB ) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminbb.pdb");

			// finally minimize everything
			for ( auto iter = minimize_movemap->movemap_torsion_id_begin();
					iter != minimize_movemap->movemap_torsion_id_end(); ++iter ) {
				if ( iter->first.second == core::id::JUMP ) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminjump.pdb");
		}
	}



	// reset the Monte Carlo object
	mc.reset(*pose_copy);

	protocols::canonical_sampling::PDBTrajectoryRecorder trajectory;
	if ( option[ OptionKeys::backrub::trajectory ] ) {
		trajectory.file_name(output_tag + "_traj.pdb" + (option[ OptionKeys::backrub::trajectory_gz ] ? ".gz" : ""));
		trajectory.stride(option[ OptionKeys::backrub::trajectory_stride ]);
		trajectory.reset(mc);
	}

	TR << "Running " << ntrials_ << " trials..." << std::endl;

	for ( core::Size i = 1; i <= ntrials_; ++i ) {

		std::string move_type;

		// could use random mover for this...
		core::Real move_prob = numeric::random::rg().uniform();
		if ( move_prob > sm_prob_ + sc_prob_ ) {
			backrubmover_->apply(*pose_copy);
			move_type = backrubmover_->type();
		} else if ( move_prob > sc_prob_ ) {
			smallmover_->apply(*pose_copy);
			move_type = smallmover_->type();
		} else {
			sidechainmover_->apply(*pose_copy);
			move_type = sidechainmover_->type();
		}

		mc.boltzmann(*pose_copy, move_type);

		if ( option[ OptionKeys::backrub::trajectory ] ) trajectory.update_after_boltzmann(mc);
	}

	mc.show_counters();

	// repack II: LAST
	// this could also be done IN the loop, so that we get a minimized structure right after the monte carlo run
	//packrotamersmover->apply(*pose);

	// dump out the low score and last accepted poses
	TR << "Last Score:" << std::endl;
	scorefxn_->show(TR, *pose_copy);
	TR.flush();

	*pose_copy = mc.lowest_score_pose();

	// repack II: LOW
	//packrotamersmover->apply(*pose);

	TR << "Low Score:" << std::endl;
	scorefxn_->show(TR, *pose_copy);
	TR.flush();

	pose= mc.lowest_score_pose();

	mc.last_accepted_pose().dump_pdb(output_tag + "_last.pdb");
	mc.lowest_score_pose().dump_pdb(output_tag + "_low.pdb");

	if ( custom_fold_tree ) {
		append_fold_tree_to_file(mc.lowest_score_pose().fold_tree(),  output_tag + "_low.pdb"); // this is the lowest scoring from MC trials
		append_fold_tree_to_file(mc.last_accepted_pose().fold_tree(), output_tag + "_last.pdb"); // this is the last accepted pose from MC trials
	}
}

///Creator
protocols::moves::MoverOP
BackrubProtocolCreator::create_mover() const {
	return protocols::moves::MoverOP(new BackrubProtocol);
}

std::string
BackrubProtocolCreator::keyname() const {
	return BackrubProtocolCreator::mover_name();
}

std::string
BackrubProtocolCreator::mover_name() {
	return "BackrubProtocol";
}

} //backrub
} //protocols
