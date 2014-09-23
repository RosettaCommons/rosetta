// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1
#include <protocols/loop_modeling/LoopModeler.hh>
#include <protocols/loop_modeling/LoopModelerCreator.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopBuilder.hh>
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.hh>

// Protocols headers
#include <protocols/filters/Filter.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/utilities/PeriodicMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>

// Core headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

// C++ headers
#include <iostream>
#include <cmath>

#define foreach BOOST_FOREACH

// Namespaces {{{1
using namespace std;

using core::Real;
using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::filters::FilterOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarloOP;
using core::util::switch_to_residue_type_set;
using core::chemical::CENTROID;
using core::chemical::FA_STANDARD;

using utility::tag::TagCOP;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
// }}}1

namespace protocols {
namespace loop_modeling {

MoverOP LoopModelerCreator::create_mover() const { // {{{1
	return new LoopModeler;
}

string LoopModelerCreator::keyname() const { // {{{1
	return "LoopModeler";
}
// }}}1

LoopModeler::LoopModeler() { // {{{1
	dont_manage_score_function();

	build_stage_ = utility::pointer::static_pointer_cast< LoopBuilder >( register_nested_loop_mover( LoopMoverOP(new LoopBuilder) ) );
	centroid_stage_ = utility::pointer::static_pointer_cast< LoopProtocol >( register_nested_loop_mover( LoopMoverOP(new LoopProtocol) ) );
	fullatom_stage_ = utility::pointer::static_pointer_cast< LoopProtocol >( register_nested_loop_mover( LoopMoverOP(new LoopProtocol) ) );

	is_build_stage_enabled_ = true;
	is_centroid_stage_enabled_ = true;
	is_fullatom_stage_enabled_ = true;
	repack_everything_before_fullatom_ = false;

	setup_kic_config();
}

LoopModeler::~LoopModeler() {} // {{{1

// }}}1

void LoopModeler::parse_my_tag( // {{{1
		TagCOP tag,
		DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose) {

	// Parse the 'config' option.

	string config = tag->getOption<string>("config", "");

	if (config == "")           { /* Don't change anything. */ }
	else if (config == "empty") { setup_empty_config(); }
	else if (config == "kic")   { setup_kic_config(); }
	else if (config == "ngk")   { setup_next_gen_kic_config(); }
	else { 
		stringstream message;
		message << "Unknown <LoopModeler config option: '" << config << "'";
		throw utility::excn::EXCN_Msg_Exception(message.str());
	}

	// Parse the 'auto_refine' option.

	if (! tag->getOption<bool>("auto_refine", true)) {
		clear_refiners();
	}

	// Parse the 'loops_file' option.

	if (tag->hasOption("loops_file")) {
		string loops_file = tag->getOption<string>("loops_file");
		set_loops(Loops(loops_file));
	}

	// Parse subtags.

	foreach (TagCOP subtag, tag->getTags()) {

		// Parse <Build> subtags.

		if (subtag->getName() == "Build") {
			build_stage_->parse_my_tag(subtag, data, filters, movers, pose);
			is_build_stage_enabled_ = ! subtag->getOption<bool>("skip", ! is_build_stage_enabled_);
		}

		// Parse <Centroid> subtags.

		else if (subtag->getName() == "Centroid") {
			centroid_stage_->parse_my_tag(subtag, data, filters, movers, pose);
			is_centroid_stage_enabled_ = ! subtag->getOption<bool>("skip", ! is_centroid_stage_enabled_);
		}

		// Parse <Fullatom> subtags.

		else if (subtag->getName() == "Fullatom") {
			fullatom_stage_->parse_my_tag(subtag, data, filters, movers, pose);
			is_fullatom_stage_enabled_ = ! subtag->getOption<bool>("skip", ! is_fullatom_stage_enabled_);
		}

		// Parse LoopMover subtags.

		else {
			LoopMoverOP loop_mover = utilities::loop_mover_from_tag(subtag, data, filters, movers, pose);
			add_shared_mover(loop_mover);
		}
	}
}
// }}}1

bool LoopModeler::do_apply(Pose & pose) { // {{{1
	Pose original_pose(pose);

	// Converting the pose to centroid mode can sometimes cause headaches with 
	// non-canonical reside type sets, so this conditional makes sure that this 
	// conversion only happens if necessary.
	
	if (is_build_stage_enabled_ || is_centroid_stage_enabled_) {
		convert_to_centroid(pose);
		if (is_build_stage_enabled_) { build_stage_->apply(pose); }
		if (is_centroid_stage_enabled_) { centroid_stage_->apply(pose); }
	}

	convert_to_fullatom(pose, original_pose);
	if (is_fullatom_stage_enabled_) { fullatom_stage_->apply(pose); }

	return true;
}

void LoopModeler::convert_to_centroid(Pose & pose) { // {{{1
	if (! pose.is_centroid()) {
		switch_to_residue_type_set(pose, CENTROID);
	}
}

void LoopModeler::convert_to_fullatom(Pose & pose, Pose & original_pose) { // {{{1
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using core::chemical::DISULFIDE;
	using core::kinematics::MoveMap;
	using core::kinematics::MoveMapOP;
	using core::optimization::AtomTreeMinimizer;
	using core::optimization::AtomTreeMinimizerOP;
	using core::optimization::MinimizerOptions;
	using core::optimization::symmetry::SymAtomTreeMinimizer;
	using core::pose::symmetry::is_symmetric;
	using core::pose::symmetry::make_residue_mask_symmetric;
	using core::pose::symmetry::make_symmetric_movemap;
	using protocols::simple_moves::PackRotamersMover;
	using protocols::simple_moves::PackRotamersMoverOP;
	using protocols::simple_moves::symmetry::SymPackRotamersMover;

	// If the pose is already fullatom, this means that the centroid stages were 
	// skipped and that the pose has been fullatom the whole time.  So we don't 
	// need to do anything unless a full repack has been requested.

	if (pose.is_fullatom() && ! repack_everything_before_fullatom_) { return; }

	// Convert the pose to fullatom mode.

	switch_to_residue_type_set(pose, FA_STANDARD);

	// Decide which sidechains (if any) to copy directly from the input pose.  No 
	// sidechains will be copied if the input pose is in centroid mode or if a 
	// full repack was requested.  Otherwise, all sidechains that are from 
	// outside the loop will be copied.
	
	vector1<bool> residues_to_repack(pose.total_residue(), true);

	if (original_pose.is_fullatom() && ! repack_everything_before_fullatom_) {
		for (Size i = 1; i < pose.total_residue(); i++) {
			if (! get_loops().is_loop_residue(i)) {
				pose.replace_residue(i, original_pose.residue(i), true);
				residues_to_repack[i] = false;
			}
		}
	}

	// Setup a packer task that will pack any residues that weren't copied from 
	// the input pose.  Disulfides are explicitly not packed.

	pose.conformation().detect_disulfides();

	TaskFactoryOP task_factory = new TaskFactory;
	task_factory->push_back(new InitializeFromCommandline);
	task_factory->push_back(new RestrictToRepacking);
	task_factory->push_back(new NoRepackDisulfides);

	PackerTaskOP task = task_factory->create_task_and_apply_taskoperations(pose);
	make_residue_mask_symmetric(pose, residues_to_repack);
	task->restrict_to_residues(residues_to_repack);

	// Setup a move map for the minimizer that will allow only sidechain DOFs to 
	// move.  Again, disulfides are explicitly left in place.

	MoveMap move_map;
	move_map.set_bb(false);
	move_map.set_chi(true);

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if (pose.residue(i).has_variant_type(DISULFIDE)) {
			move_map.set_chi(i, false);
		}
	}

	// Apply the packer and the minimizer, accounting for symmetry.

	PackRotamersMoverOP packer;
	AtomTreeMinimizerOP minimizer;
	MinimizerOptions min_options("dfpmin_armijo_nonmonotone", 1e-5, true, false);
	ScoreFunctionCOP fa_score_function = fullatom_stage_->get_score_function();

	if (is_symmetric(pose)) {
		packer = new SymPackRotamersMover(fa_score_function, task);
		minimizer = new SymAtomTreeMinimizer;
		make_symmetric_movemap(pose, move_map);
	} else {
		packer = new PackRotamersMover(fa_score_function, task);
		minimizer = new AtomTreeMinimizer;
	}

	packer->apply(pose);
	minimizer->run(pose, move_map, *fa_score_function, min_options);
}

// }}}1

void LoopModeler::setup_empty_config() { // {{{1
	build_stage_->set_score_function(loops::get_cen_scorefxn());

	centroid_stage_->set_sfxn_cycles(3);
	centroid_stage_->set_temp_cycles(20, true);
	centroid_stage_->set_mover_cycles(1);
	centroid_stage_->set_score_function(loops::get_cen_scorefxn());
	centroid_stage_->clear_movers_and_refiners();

	fullatom_stage_->set_sfxn_cycles(3);
	fullatom_stage_->set_temp_cycles(10, true);
	fullatom_stage_->set_mover_cycles(2);
	fullatom_stage_->set_score_function(loops::get_fa_scorefxn());
	fullatom_stage_->clear_movers_and_refiners();
}

void LoopModeler::setup_kic_config() { // {{{1
	using namespace protocols::kinematic_closure;
	using namespace protocols::loop_modeling::refiners;
	using namespace protocols::loop_modeling::utilities;

	setup_empty_config();

	add_shared_mover(new KicMover);

	centroid_stage_->add_refiner(new MinimizationRefiner);
	fullatom_stage_->add_refiner(new PeriodicMover(LoopMoverOP( new RepackingRefiner ), 20));
	fullatom_stage_->add_refiner(new RotamerTrialsRefiner);
	fullatom_stage_->add_refiner(new MinimizationRefiner);

	mark_as_default();
}

void LoopModeler::setup_next_gen_kic_config() { // {{{1
	using namespace protocols::kinematic_closure;

	utility_exit_with_message("Next-generation KIC not yet supported by LoopModeler.");

	/*
	setup_kic_config();

	KicMoverOP kic_mover = new KicMover;
	kic_mover->add_perturber(new BondAnglePerturber);
	kic_mover->add_perturber(new TabooPerturber);
	kic_mover->add_perturber(new NeighborDependentRamaPerturber);

	add_shared_mover(kic_mover);

	fullatom_stage_->set_ramp_rama(true);
	fullatom_stage_->set_ramp_repulsive(true);

	mark_as_default();
	*/
}
// }}}1

void LoopModeler::enable_build_stage() { // {{{1
	is_build_stage_enabled_ = true;
}

void LoopModeler::disable_build_stage() { // {{{1
	is_build_stage_enabled_ = false;
}

void LoopModeler::enable_centroid_stage() { // {{{1
	is_centroid_stage_enabled_ = true;
}

void LoopModeler::disable_centroid_stage() { // {{{1
	is_centroid_stage_enabled_ = false;
}

void LoopModeler::enable_fullatom_stage() { // {{{1
	is_fullatom_stage_enabled_ = true;
}

void LoopModeler::disable_fullatom_stage() { // {{{1
	is_fullatom_stage_enabled_ = false;
}

LoopBuilderOP LoopModeler::build_stage() { // {{{1
	return build_stage_;
}

LoopProtocolOP LoopModeler::centroid_stage() { // {{{1
	return centroid_stage_;
}

LoopProtocolOP LoopModeler::fullatom_stage() { // {{{1
	return fullatom_stage_;
}
// }}}1

void LoopModeler::repack_everything_before_fullatom(bool setting) { // {{{1
	repack_everything_before_fullatom_ = setting;
}
// }}}1

void LoopModeler::add_shared_mover(LoopMoverOP mover) { // {{{1
	centroid_stage_->add_mover(mover);
	fullatom_stage_->add_mover(mover);
}

void LoopModeler::add_shared_refiner(LoopMoverOP refiner) { // {{{1
	centroid_stage_->add_refiner(refiner);
	fullatom_stage_->add_refiner(refiner);
}

void LoopModeler::add_shared_filter(FilterOP filter) { // {{{1
	centroid_stage_->add_filter(filter);
	fullatom_stage_->add_filter(filter);
}

void LoopModeler::clear_movers() { // {{{1
	centroid_stage_->clear_movers();
	fullatom_stage_->clear_movers();
}

void LoopModeler::clear_refiners() { // {{{1
	centroid_stage_->clear_refiners();
	fullatom_stage_->clear_refiners();
}

void LoopModeler::clear_movers_and_refiners() { // {{{1
	centroid_stage_->clear_movers_and_refiners();
	fullatom_stage_->clear_movers_and_refiners();
}

void LoopModeler::mark_as_default() { // {{{1
	centroid_stage_->mark_as_default();
	fullatom_stage_->mark_as_default();
}
// }}}1

}
}

