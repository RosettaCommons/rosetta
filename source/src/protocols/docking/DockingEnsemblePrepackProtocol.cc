// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingEnsemblePrepackProtocol.cc
/// @brief Prepacking of the bound structure before docking with ensembles
/// @author Monica Berrondo

// Unit Headers
#include <protocols/docking/DockingEnsemblePrepackProtocol.hh>

// Package headers
#include <protocols/docking/util.hh>
#include <protocols/docking/DockingEnsemble.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <protocols/docking/SidechainMinMover.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>


#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>


using basic::T;
using namespace protocols::moves;
using namespace core;
using namespace pack::task;
using protocols::jd2::JobDistributor;

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.DockingEnsemblePrepackProtocol" );

namespace protocols {
namespace docking {

DockingEnsemblePrepackProtocol::DockingEnsemblePrepackProtocol(): DockingHighRes()
{
	Mover::type( "DockingEnsemblePrepackProtocol" );
	setup_defaults();
	register_options();
	init_from_options();
}


void DockingEnsemblePrepackProtocol::setup_defaults()
{
	trans_magnitude_ = 1000.0;
	pack_operations_ = SequenceMoverOP( new SequenceMover() );

	ensemble1_ = NULL;
	ensemble2_ = NULL;
	ensemble1_filename_ = "";
	ensemble2_filename_ = "";
}

DockingEnsemblePrepackProtocol::~DockingEnsemblePrepackProtocol(){
}

void DockingEnsemblePrepackProtocol::init_from_options()
{
	using namespace basic::options;
	if ( option[ OptionKeys::docking::dock_rtmin ].user() ) {
		set_rt_min(option[ OptionKeys::docking::dock_rtmin ]());
	}

	if ( option[ OptionKeys::docking::sc_min ].user() ) {
		set_sc_min(option[ OptionKeys::docking::sc_min ]());
	}

	if ( option[ OptionKeys::docking::partners ].user() ) {
		set_partners(option[ OptionKeys::docking::partners ]());
	}

	if ( option[ OptionKeys::docking::ensemble1 ].user() ) {
		set_ensemble1(option[ OptionKeys::docking::ensemble1 ]());
	}

	if ( option[ OptionKeys::docking::ensemble2 ].user() ) {
		set_ensemble2(option[ OptionKeys::docking::ensemble2 ]());
	}
}

void DockingEnsemblePrepackProtocol::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::dock_rtmin );
	option.add_relevant( OptionKeys::docking::sc_min );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::docking::ensemble1 );
	option.add_relevant( OptionKeys::docking::ensemble2 );
}

void DockingEnsemblePrepackProtocol::setup_pack_operation_movers()
{
	prepack_full_repack_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover() );
	prepack_full_repack_->score_function( scorefxn_pack() );
	prepack_full_repack_->task_factory( task_factory() );
	pack_operations_->add_mover(prepack_full_repack_);

	if ( rt_min() ) {
		rtmin_mover_ = protocols::simple_moves::RotamerTrialsMinMoverOP( new protocols::simple_moves::RotamerTrialsMinMover( ) );
		rtmin_mover_->score_function( scorefxn_pack() );
		rtmin_mover_->task_factory( task_factory() );
		pack_operations_->add_mover( rtmin_mover_ );
	}
	if ( sc_min() ) {
		scmin_mover_ = SidechainMinMoverOP( new SidechainMinMover() );
		scmin_mover_->set_scorefxn( scorefxn_pack() );
		scmin_mover_->set_task_factory( task_factory() );
		pack_operations_->add_mover( scmin_mover_ );
	}

}

void DockingEnsemblePrepackProtocol::finalize_setup( pose::Pose & pose ) {
	setup_foldtree( pose, partners(), movable_jumps() );
	tf2()->set_prepack_only(true);
	tf2()->create_and_attach_task_factory( this, pose );
	setup_pack_operation_movers();

	core::scoring::ScoreFunctionOP scorefxn_low = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );

	if ( ensemble1_filename_ == "" || ensemble2_filename_ == "" ) utility_exit_with_message( "Must define ensemble file for both partners");
	runtime_assert( movable_jumps().size() == 1 ); // ensemble mode is only allowed for traditional docking
	core::Size const rb_jump = movable_jumps()[1];
	Size start_res(1), end_res(1), cutpoint(pose.fold_tree().cutpoint_by_jump( rb_jump ));

	TR << "Ensemble 1: " << ensemble1_filename_ << std::endl;
	start_res = 1;
	end_res = cutpoint;

	ensemble1_ = DockingEnsembleOP( new DockingEnsemble( start_res, end_res, rb_jump, ensemble1_filename_, "dock_ens_conf1", scorefxn_low, scorefxn() ) );

	TR << "Ensemble 2: " << ensemble2_filename_ << std::endl;
	start_res = cutpoint + 1;
	end_res = pose.total_residue();

	ensemble2_ = DockingEnsembleOP( new DockingEnsemble( start_res, end_res, rb_jump, ensemble2_filename_, "dock_ens_conf2", scorefxn_low, scorefxn() ) );
}

void DockingEnsemblePrepackProtocol::apply( core::pose::Pose & pose )
{
	finalize_setup(pose);
	protocols::docking::ConformerSwitchMoverOP switch_mover;
	core::pose::Pose starting_pose;

	starting_pose = pose;

	switch_mover = protocols::docking::ConformerSwitchMoverOP( new protocols::docking::ConformerSwitchMover( ensemble1_ ) );
	for ( Size i=1; i<=ensemble1_->size(); ++i ) {
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( pose );
		switch_mover->switch_conformer( pose, i );

		//Move each partners away from the others
		for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		}
		ensemble1_->calculate_lowres_ref_energy( pose );
		//bringing the packed structures together
		for (  DockJumps::const_iterator jump= movable_jumps().begin() ; jump != movable_jumps().end(); ++jump ) {
			rigid::RigidBodyTransMoverOP translate_back( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_back->step_size( trans_magnitude_ );
			translate_back->trans_axis().negate();
			translate_back->apply(pose);
		}

		ensemble1_->set_packer( pack_operations_ );

		ensemble1_->calculate_highres_ref_energy( i ); // also does the dump_pdb
	}
	ensemble1_->update_pdblist_file();

	// reset to starting pose
	pose = starting_pose;
	switch_mover = protocols::docking::ConformerSwitchMoverOP( new protocols::docking::ConformerSwitchMover( ensemble2_ ) );
	for ( Size i=1; i<=ensemble2_->size(); ++i ) {
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( pose );

		switch_mover->switch_conformer( pose, i );

		//Move each partners away from the others
		for ( DockJumps::const_iterator jump = movable_jumps().begin() ; jump != movable_jumps().end() ; ++jump ) {
			rigid::RigidBodyTransMoverOP translate_away( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_away->step_size( trans_magnitude_ );
			translate_away->apply(pose);
		}
		ensemble2_->calculate_lowres_ref_energy( pose );
		//bringing the packed structures together
		for (  DockJumps::const_iterator jump= movable_jumps().begin() ; jump != movable_jumps().end(); ++jump ) {
			rigid::RigidBodyTransMoverOP translate_back( new rigid::RigidBodyTransMover(pose, *jump) );
			translate_back->step_size( trans_magnitude_ );
			translate_back->trans_axis().negate();
			translate_back->apply(pose);
		}

		ensemble2_->set_packer( pack_operations_ );

		ensemble2_->calculate_highres_ref_energy( i ); // also does the dump_pdb
	}
	ensemble2_->update_pdblist_file();
}

std::string DockingEnsemblePrepackProtocol::get_name() const {
	return "DockingEnsemblePrepackProtocol";
}

}
}
