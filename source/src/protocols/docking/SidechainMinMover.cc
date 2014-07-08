// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Sid Chaudhury

// Unit headers
#include <protocols/docking/SidechainMinMover.hh>

// Package Headers
#include <protocols/docking/DockTaskFactory.hh>

// Project headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/scoring/Interface.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/MinMover.hh>

// Utility Headers
#include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.docking.SidechainMinMover");

namespace protocols {
namespace docking {

//default constructor
SidechainMinMover::SidechainMinMover() : DockingHighRes()
{
	DockingHighRes::type( "SidechainMin" );
	update_movemap_ = true;
	core::scoring::ScoreFunctionOP temp_sf = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	set_scorefxn( temp_sf );
	set_default();
}

//construtor with arguments
SidechainMinMover::SidechainMinMover( core::scoring::ScoreFunctionOP scorefxn) : DockingHighRes( 1, scorefxn )
{
	DockingHighRes::type( "SidechainMin" );
	//scorefxn_ = new core::scoring::ScoreFunction( *scorefxn() );
	set_default();
}

//construtor with arguments
SidechainMinMover::SidechainMinMover(
	core::scoring::ScoreFunctionOP scorefxn,
	core::kinematics::MoveMapOP movemap
) : DockingHighRes( 1, scorefxn ), movemap_(movemap), update_movemap_(false)
{
	DockingHighRes::type( "SidechainMin" );
	//scorefxn_ = new core::scoring::ScoreFunction( *scorefxn );
	set_default();
}

//constructor with arguments
SidechainMinMover::SidechainMinMover(
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::PackerTaskOP task) : DockingHighRes( 1, scorefxn ), task_(task), update_movemap_(true)
{
	DockingHighRes::type( "SidechainMin" );
	//scorefxn_ = new core::scoring::ScoreFunction( *scorefxn() );
	set_default();
}

//constructor with arguments
SidechainMinMover::SidechainMinMover(
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryCOP tf) : DockingHighRes( 1, scorefxn), update_movemap_(true)
{
	DockingHighRes::type( "SidechainMin" );
	//scorefxn_ = new core::scoring::ScoreFunction( *scorefxn() );
	//tf_ = new core::pack::task::TaskFactory( *tf );
	set_task_factory( tf );
	set_default();
}

SidechainMinMover::SidechainMinMover(
	core::Size rb_jump, core::scoring::ScoreFunctionOP scorefxn ): DockingHighRes( rb_jump, scorefxn)
{
	DockingHighRes::type( "SidechainMin" );
	set_default();
}
//destructor
SidechainMinMover::~SidechainMinMover() {}

void SidechainMinMover::set_minmover( protocols::simple_moves::MinMoverOP minmover ){ minmover_ = minmover; }

//default options setup for SidechainMinMover
void SidechainMinMover::set_default()
{
	if(update_movemap_){
		movemap_ = new core::kinematics::MoveMap();
		movemap_->set_chi( true );
	}
	minmover_ = new protocols::simple_moves::MinMover(movemap_, scorefxn(), "dfpmin_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/  );
}

std::string
SidechainMinMover::get_name() const {
	return "SidechainMin";
}

//SidechainMinMover updates the movemap based on a packer task
void SidechainMinMover::update_movemap( core::pose::Pose & pose)
{
	movemap_->set_chi( true );
	if ( task_factory() ) task_ = task_factory()->create_task_and_apply_taskoperations( pose );
	if ( task_ ){
		for(Size i = 1; i <= pose.total_residue(); i++){
			if (!task_->nonconst_residue_task(i).being_packed()) movemap_->set_chi(i, false);
		}
	}
}

//SidechainMinMover apply function
void SidechainMinMover::apply( core::pose::Pose & pose )
{
	//runtime_assert(pose.is_fullatom());
	if(update_movemap_) update_movemap( pose );
	minmover_->apply( pose );
}

/////////////////////////// InterfaceSidechainMinMover ////////////////////////////////

//constructor
InterfaceSidechainMinMover::InterfaceSidechainMinMover() : SidechainMinMover(), interface_dist_(8.0)
{
	Mover::type( "InterfaceSidechainMinMover" );
	set_default();
}

InterfaceSidechainMinMover::InterfaceSidechainMinMover(
	core::Size rb_jump,
	core::scoring::ScoreFunctionOP scorefxn,
	core::Real interface_dist
) : SidechainMinMover(rb_jump, scorefxn), interface_dist_(interface_dist)
{
	Mover::type( "InterfaceSidechainMin" );
	set_default();
}

//destructor
InterfaceSidechainMinMover::~InterfaceSidechainMinMover() {}

void InterfaceSidechainMinMover::set_interface_dist( core::Real interface_dist)
{
	interface_dist_ = interface_dist;
	interface_->distance(interface_dist_);
}

//default setup for InterfaceSidechainMinMover
void
InterfaceSidechainMinMover::set_default()
{
	interface_ = new protocols::scoring::Interface( movable_jumps()[1] );
}

//apply function for InterfaceSidechainMinMover
void
InterfaceSidechainMinMover::apply( core::pose::Pose & pose )
{
	//reset movemap
	movemap_->set_bb(false);
	movemap_->set_jump(false);
	movemap_->set_chi(false);

	//calculate interface
	interface_->calculate( pose );

	Size cutpoint ( pose.fold_tree().cutpoint_by_jump( movable_jumps()[1] ) );

	//set move map for interface residues only
	if ( !( tf2() )->get_norepack1() ){
		for(core::Size i = 1; i <= cutpoint; i++){
			if (interface_->is_interface(i)) movemap_->set_chi(i, true);
		}
	}

	if ( !( tf2() )->get_norepack2() ){
		for(core::Size i = cutpoint+1; i <= pose.total_residue(); i++){
			if (interface_->is_interface(i)) movemap_->set_chi(i, true);
		}
	}

	minmover_->apply(pose);

}

std::string
InterfaceSidechainMinMover::get_name() const {
	return "InterfaceSidechainMin";
}

} // docking
} // protocols
