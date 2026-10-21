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
/// @author Kyle Barlow
/// @details
/// Currently a work in progress. The goal is to match the features of rosetta++ -backrub_mc


// Project Headers
#include <protocols/backrub/BackrubProtocol.hh>
#include <protocols/backrub/BackrubProtocolCreator.hh>
#include <protocols/backrub/util.hh>

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jd2/util.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/viewer/viewers.hh>

// Core Headers
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Platform Headers

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR("protocols.backrub.BackrubProtocol");

namespace protocols {
namespace backrub {

BackrubProtocol::BackrubProtocol():
	Mover(),
	scorefxn_(nullptr),
	main_task_factory_(nullptr),
	backrubmover_(utility::pointer::make_shared< BackrubMover >()),
	smallmover_(utility::pointer::make_shared< protocols::simple_moves::SmallMover >()),
	sidechainmover_(utility::pointer::make_shared< protocols::simple_moves::sidechain_moves::SidechainMover >()),
	packrotamersmover_(utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >()),
	movemap_smallmover_(nullptr),
	minimize_movemap_(nullptr),
	packing_operation_(nullptr),
	trajectory_(false),
	trajectory_gz_(false),
	recover_low_(true),
	trajectory_stride_(100),
	trajectory_apply_mover_(nullptr),
	pivots_residue_selector_(nullptr)
{
	read_cmd_line_options();
}

BackrubProtocol::BackrubProtocol(BackrubProtocol const & ) = default;

BackrubProtocol::~BackrubProtocol()= default;

void
BackrubProtocol::set_options(
	utility::vector1<core::Size> pivot_residues,
	utility::vector1<std::string> pivot_atoms,
	core::kinematics::MoveMapCOP minimize_movemap,
	core::kinematics::MoveMapCOP movemap_smallmover,
	core::pack::task::operation::TaskOperationCOP packing_operation,
	core::Size min_atoms,
	core::Size max_atoms,
	bool initial_pack,
	core::Real mm_bend_weight,
	core::Real sm_prob,
	core::Real sc_prob,
	core::Real sc_prob_uniform,
	core::Real sc_prob_withinrot,
	core::Real mc_kt,
	core::Size ntrials,
	bool trajectory,
	bool trajectory_gz,
	core::Size trajectory_stride
) {
	this->set_pivot_residues( pivot_residues );
	this->set_pivot_atoms( pivot_atoms );
	this->set_movemap_smallmover(movemap_smallmover);

	minimize_movemap_ = minimize_movemap;
	packing_operation_ = packing_operation;
	min_atoms_ = min_atoms;
	max_atoms_ = max_atoms;
	initial_pack_ = initial_pack;
	mm_bend_weight_ = mm_bend_weight;
	sm_prob_ = sm_prob;
	sc_prob_ = sc_prob;
	sc_prob_uniform_ = sc_prob_uniform;
	sc_prob_withinrot_ = sc_prob_withinrot;
	mc_kt_ = mc_kt;
	ntrials_ = ntrials;
	trajectory_ = trajectory;
	trajectory_gz_ = trajectory_gz;
	trajectory_stride_ = trajectory_stride;
}

void
BackrubProtocol::read_cmd_line_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1<core::Size> pivot_residues;

	if ( option[ OptionKeys::backrub::pivot_residues].user() ) {
		for ( core::Size i = 1; i <= option[ OptionKeys::backrub::pivot_residues ].size(); ++i ) {
			if ( option[ OptionKeys::backrub::pivot_residues ][i] >= 1 ) pivot_residues.push_back(option[ OptionKeys::backrub::pivot_residues ][i]);
		}
	}

	core::kinematics::MoveMapOP minimize_movemap = nullptr;
	if ( option[ OptionKeys::backrub::minimize_movemap ].user() ) {
		minimize_movemap = utility::pointer::make_shared< core::kinematics::MoveMap >();
		minimize_movemap->init_from_file(option[ OptionKeys::backrub::minimize_movemap ]);
	}

	core::pack::task::operation::TaskOperationCOP packing_operation = nullptr;
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		packing_operation = utility::pointer::make_shared< core::pack::task::operation::ReadResfile >();
	}

	core::Real sm_prob = option[ OptionKeys::backrub::sm_prob ];
	core::kinematics::MoveMapOP movemap_smallmover = nullptr;
	if ( sm_prob > 0 ) {
		if ( ! movemap_smallmover_ ) {
			movemap_smallmover = utility::pointer::make_shared< core::kinematics::MoveMap >();
			movemap_smallmover->init_from_file( basic::options::option[ basic::options::OptionKeys::in::file::movemap ] );
		}
	}

	this->set_options(
		pivot_residues,
		option[ OptionKeys::backrub::pivot_atoms ],
		minimize_movemap,
		movemap_smallmover,
		packing_operation,
		option[ OptionKeys::backrub::min_atoms ],
		option[ OptionKeys::backrub::max_atoms ],
		option[ OptionKeys::backrub::initial_pack ],
		option[ OptionKeys::backrub::mm_bend_weight ],
		option[ OptionKeys::backrub::sm_prob ],
		option[ OptionKeys::backrub::sc_prob ],
		option[ OptionKeys::backrub::sc_prob_uniform ],
		option[ OptionKeys::backrub::sc_prob_withinrot ],
		option[ OptionKeys::backrub::mc_kt],
		option[ OptionKeys::backrub::ntrials ],
		option[ OptionKeys::backrub::trajectory ],
		option[ OptionKeys::backrub::trajectory_gz ],
		option[ OptionKeys::backrub::trajectory_stride ]
	);
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
	this->set_pivot_residues( get_pivot_residues_from_movemap(movemap) );
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

void
BackrubProtocol::set_dump_poses(bool setting) {
	dump_poses_  = setting;
}

protocols::moves::MoverOP
BackrubProtocol::clone() const {
	return utility::pointer::make_shared< BackrubProtocol >( *this );
}

protocols::moves::MoverOP
BackrubProtocol::fresh_instance() const {
	return utility::pointer::make_shared< BackrubProtocol >();
}

void BackrubProtocol::set_backrub_mover(BackrubMoverOP backrub_mover){
	backrubmover_ = backrub_mover;
}

void
BackrubProtocol::write_database() {
	backrubmover_->branchopt().write_database();
}

void
BackrubProtocol::finalize_setup(core::pose::Pose & pose){
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	if ( ! main_task_factory_ ) {
		main_task_factory_ = utility::pointer::make_shared< core::pack::task::TaskFactory >();
		main_task_factory_->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
		if ( packing_operation_ ) {
			main_task_factory_->push_back( packing_operation_ );
		} else {
			operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
			main_task_factory_->push_back( rtrop );
		}
	}

	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory_->push_back( utility::pointer::make_shared< operation::PreserveCBeta >() );

	// set up the score function and add the bond angle energy term

	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}
	scorefxn_->set_weight_if_zero(core::scoring::mm_bend, mm_bend_weight_);
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

	if ( pivots_residue_selector_ ) {
		set_pivots_from_residue_subset( pivots_residue_selector_->apply( pose ) );
	}

	if ( pivot_residues_.size() > 0 ) {
		backrubmover_->set_pivot_residues( pivot_residues_ );
	}

	backrubmover_->clear_segments();
	backrubmover_->set_pivot_atoms( pivot_atoms_ );
	backrubmover_->set_min_atoms( min_atoms_ );
	backrubmover_->set_max_atoms( max_atoms_ );
	backrubmover_->set_input_pose( utility::pointer::make_shared< core::pose::Pose >( pose ) );
	backrubmover_->add_mainchain_segments( pose );

}

void
BackrubProtocol::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap  & data
) {

	utility::vector1<std::string> pivot_atoms = utility::vector1<std::string>(1, "CA");
	if ( tag->hasOption("pivot_atoms") ) {
		std::string const pivot_atoms_string( tag->getOption<std::string>("pivot_atoms") );
		pivot_atoms = utility::string_split( pivot_atoms_string, ',' );
	}

	if ( tag->hasOption("recover_low") ) {
		recover_low_ = tag->getOption<bool>("recover_low");
	}
	dump_poses_ = tag->getOption<bool>("dump_poses", false); // For XML, default to not dumping the pose

	if ( tag->hasOption("task_operations") ) {
		main_task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
	}

	// Take any residues set to be at least packable in "pivot_task_operations" and use these as pivots
	if ( tag->hasOption("pivot_residues") ) {
		pivots_residue_selector_ = core::pose::get_resnum_selector(tag, "pivot_residues");
	}
	if ( tag->hasOption("pivot_residue_selector") ) {
		// We hold on to the residue selector until apply time (to run its apply function), as it has not been
		// initialized at this point and we would need a non-const pose to initialize
		pivots_residue_selector_ = protocols::rosetta_scripts::get_residue_selector( tag->getOption<std::string>("pivot_residue_selector"), data);
	}

	auto mc_kt = tag->getOption<core::Real>(
		"mc_kt",
		basic::options::option[ basic::options::OptionKeys::backrub::mc_kt ]()
	);

	auto ntrials = tag->getOption<core::Size>(
		"ntrials",
		basic::options::option[ basic::options::OptionKeys::backrub::ntrials ]()
	);


	if ( tag->hasOption("trajectory_apply_mover") ) {
		moves::MoverOP mover = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "trajectory_apply_mover" ), data );
		trajectory_apply_mover_ = utility::pointer::dynamic_pointer_cast < protocols::moves::Mover > (mover);
	}

	auto trajectory_stride = tag->getOption<core::Size>(
		"trajectory_stride",
		basic::options::option[ basic::options::OptionKeys::backrub::trajectory_stride ]()
	);

	bool trajectory = tag->getOption<bool>(
		"trajectory",
		basic::options::option[ basic::options::OptionKeys::backrub::trajectory ]()
	);

	bool trajectory_gz = tag->getOption<bool>(
		"trajectory_gz",
		basic::options::option[ basic::options::OptionKeys::backrub::trajectory_gz ]()
	);
	if ( trajectory_gz ) {
		trajectory = true;
	}

	utility::vector1<core::Size> pivot_residues; // Pivot residues will be set by the ResidueSelectors
	set_options(
		pivot_residues,
		pivot_atoms,
		nullptr, // minimize_movemap - not implemented in parse_my_tag
		nullptr, // movemap_smallmover - not implemented in parse_my_tag
		nullptr, // packing_operation - not needed, as main_task_factory_ set directly if task_operations passed
		basic::options::option[ basic::options::OptionKeys::backrub::min_atoms ](), // min_atoms
		basic::options::option[ basic::options::OptionKeys::backrub::max_atoms ](), // max_atoms,
		false, // initial_pack,
		basic::options::option[ basic::options::OptionKeys::backrub::mm_bend_weight ](), // mm_bend_weight,
		basic::options::option[ basic::options::OptionKeys::backrub::sm_prob ](), // sm_prob,
		basic::options::option[ basic::options::OptionKeys::backrub::sc_prob ](), // sc_prob,
		basic::options::option[ basic::options::OptionKeys::backrub::sc_prob_uniform ](), // sc_prob_uniform,
		basic::options::option[ basic::options::OptionKeys::backrub::sc_prob_withinrot ](), // sc_prob_withinrot,
		mc_kt,
		ntrials,
		trajectory,
		trajectory_gz,
		trajectory_stride
	);
}

void
BackrubProtocol::apply( core::pose::Pose& pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;



	bool custom_fold_tree = false;
	std::string input_tag = protocols::jd2::current_input_tag();
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
		TR.Warning << "No side chains to move. Not using SidechainMover." << std::endl;
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
		if ( minimize_movemap_ ) {

			// setup the MoveMaps
			core::kinematics::MoveMapOP minimize_movemap_progressive( new core::kinematics::MoveMap );

			// setup the MinMover
			protocols::minimization_packing::MinMover minmover;
			minmover.score_function(scorefxn_);
			minmover.min_type("dfpmin");

			// first minimize just the side chains
			for ( auto iter = minimize_movemap_->movemap_torsion_id_begin();
					iter != minimize_movemap_->movemap_torsion_id_end(); ++iter ) {
				if ( iter->first.second == core::id::CHI ) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminchi.pdb");

			// next minimize the side chains and backbone
			for ( auto iter = minimize_movemap_->movemap_torsion_id_begin();
					iter != minimize_movemap_->movemap_torsion_id_end(); ++iter ) {
				if ( iter->first.second == core::id::BB ) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminbb.pdb");

			// finally minimize everything
			for ( auto iter = minimize_movemap_->movemap_torsion_id_begin();
					iter != minimize_movemap_->movemap_torsion_id_end(); ++iter ) {
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
	if ( trajectory_ ) {
		TR << "Will dump a PDB file every " << trajectory_stride_ << " steps" << std::endl;
		trajectory.file_name(output_tag + "_traj.pdb" + ( trajectory_gz_ ? ".gz" : ""));
		trajectory.stride( trajectory_stride_ );
		trajectory.reset(mc);
	}

	TR << "Running " << ntrials_ << " trials..." << std::endl;

	for ( core::Size i = 1; i <= ntrials_ ; ++i ) {

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

		if ( trajectory_ ) trajectory.update_after_boltzmann(mc);

		if ( ( trajectory_apply_mover_ != nullptr ) && ( (i % trajectory_stride_) == 0 ) ) {
			TR << "Applying mover '" << trajectory_apply_mover_->get_name() << "' on step " << i << std::endl;
			core::pose::PoseOP pose_for_traj_mover( new core::pose::Pose(*pose_copy) );
			trajectory_apply_mover_->apply( *pose_for_traj_mover );
		}
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

	if ( recover_low_ ) {
		TR << "Setting pose to be lowest scoring sampled model" << std::endl;
		pose = mc.lowest_score_pose();
	} else {
		TR << "Leaving pose as last sampled model" << std::endl;
		pose = mc.last_accepted_pose();
	}

	if ( dump_poses_ ) {
		mc.last_accepted_pose().dump_pdb(output_tag + "_last.pdb");
		mc.lowest_score_pose().dump_pdb(output_tag + "_low.pdb");
	}

	if ( custom_fold_tree ) {
		append_fold_tree_to_file(mc.lowest_score_pose().fold_tree(),  output_tag + "_low.pdb"); // this is the lowest scoring from MC trials
		append_fold_tree_to_file(mc.last_accepted_pose().fold_tree(), output_tag + "_last.pdb"); // this is the last accepted pose from MC trials
	}
}

void
BackrubProtocol::set_pivots_from_residue_subset(
	core::select::residue_selector::ResidueSubset residue_subset
) {
	// A ResidueSubset is equivalent to: typedef utility::vector1< bool >
	for ( core::Size i = 1 ; i <= residue_subset.size() ; ++i ) {
		if ( residue_subset[i] ) {
			pivot_residues_.push_back( i );
		}
	}
}

utility::vector1<core::Size>
BackrubProtocol::get_pivot_residues() const
{
	return pivot_residues_;
}

std::string BackrubProtocol::get_name() const {
	return mover_name();
}

std::string BackrubProtocol::mover_name() {
	return "BackrubProtocol";
}

void BackrubProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"pivot_residues", xs_string,
		"Pivot residues to use for backbone moves. Can contain segments (comma separated). Can use PDB numbers ([resnum][chain]) or absolute Rosetta numbers (integer)");
	attlist + XMLSchemaAttribute(
		"pivot_atoms", xs_string,
		"main chain atoms usable as pivots (comma separated)");

	rosetta_scripts::attributes_for_parse_task_operations(attlist);

	attlist + XMLSchemaAttribute(
		"mc_kt", xsct_real,
		"Temperature to use for Metropolis criterion");
	attlist + XMLSchemaAttribute(
		"ntrials", xsct_real,
		"Number of trials to perform");
	attlist + XMLSchemaAttribute(
		"trajectory", xsct_rosetta_bool,
		"Set to true to dump PDBs along the trajectory (how often is controlled by trajectory_stride option)" );
	attlist + XMLSchemaAttribute(
		"trajectory_gz", xsct_rosetta_bool,
		"Set to true to dump gzipped PDBs along the trajectory (how often is controlled by trajectory_stride option)" );
	attlist + XMLSchemaAttribute(
		"recover_low", xsct_rosetta_bool,
		"Set the return Pose to be the lowest scoring structure sampled (if false, will be left as last sampled structure)" );
	attlist + XMLSchemaAttribute(
		"dump_poses", xsct_rosetta_bool,
		"Output the results of the MC trajectory to _low.pdb and _last.pdb files in the working directory" );
	attlist + XMLSchemaAttribute(
		"pivot_residue_selector", xs_string,
		"Name of residue selector to use to select pivot residues" );
	attlist + XMLSchemaAttribute(
		"trajectory_apply_mover", xs_string,
		"Name of mover to apply during trajectory. Stride (how often to apply) is affected by trajectory_stride" );
	attlist + XMLSchemaAttribute(
		"trajectory_stride", xsct_positive_integer,
		"How often (in steps) to either dump PDBs or call trajectory_apply_mover"
	);

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Performs backrub-style backbone moves",
		attlist );
}

std::string BackrubProtocolCreator::keyname() const {
	return BackrubProtocol::mover_name();
}

protocols::moves::MoverOP
BackrubProtocolCreator::create_mover() const {
	return utility::pointer::make_shared< BackrubProtocol >();
}

void BackrubProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BackrubProtocol::provide_xml_schema( xsd );
}


} //backrub
} //protocols
