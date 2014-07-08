// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingHighRes
/// @brief protocols that are specific to high resolution docking
/// @detailed
///		This contains the functions that create initial positions for docking
///		You can either randomize partner 1 or partner 2, spin partner 2, or
///		perform a simple perturbation.
/// 	Also contains docking mcm protocol
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Sid Chaudhury
/// @author Modified by Jacob Corn


// Package Headers
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/docking/DockTaskFactory.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// For some reason, the full header has to be included in the header. since it's there, it doesn't need to be included here too.
// Ideally the forward declaration would be included in the header and the full header would be here.
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.hh>

// ObjexxFCL Headers

//Utility Headers
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


// C++ Headers

using basic::T;

using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.docking.DockingHighRes");

//     originally from dock_structure.cc Jeff Gray April 2001

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockingHighRes::DockingHighRes() : Mover()
{
	init( utility::tools::make_vector1<core::SSize>(1) ); // operate on the first jump
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	scorefxn_pack_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
}
// default constructor
DockingHighRes::DockingHighRes(
	core::Size const rb_jump
) : Mover()
{
	init( utility::tools::make_vector1<core::SSize>(rb_jump) ); // operate on the first jump
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	scorefxn_pack_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
}

// constructor with arguments
// only one movable jump
DockingHighRes::DockingHighRes(
	core::Size const rb_jump,
	core::scoring::ScoreFunctionOP scorefxn
) : Mover(), scorefxn_(scorefxn), scorefxn_pack_(scorefxn)
{
	init( utility::tools::make_vector1<core::SSize>(rb_jump) );
}

// constructor with arguments
// only one movable jump, scoring and packing defined
DockingHighRes::DockingHighRes(
	core::Size const rb_jump,
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
) : Mover(), scorefxn_(scorefxn), scorefxn_pack_(scorefxn_pack)
{
	init( utility::tools::make_vector1<core::SSize>(rb_jump) );
}

// constructor with arguments
// vector of  movable jumps, scoring and packing defined
DockingHighRes::DockingHighRes(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionOP scorefxn,
	core::scoring::ScoreFunctionOP scorefxn_pack
) : Mover(), scorefxn_(scorefxn), scorefxn_pack_(scorefxn_pack)
{
	init( movable_jumps );
}

DockingHighRes::DockingHighRes( DockingHighRes const & old_instance ) :
	//utility::pointer::ReferenceCount(),
	Mover( old_instance )
{
	sc_min_ = old_instance.sc_min_;
	rt_min_ = old_instance.rt_min_;
	partners_ = old_instance.partners_;

	scorefxn_ = old_instance.scorefxn_->clone();
	scorefxn_pack_ = old_instance.scorefxn_pack_->clone();

	movable_jumps_ = old_instance.movable_jumps_;
	if( old_instance.tf_ ){
		tf_ = new core::pack::task::TaskFactory( *old_instance.tf_ );
	}
	tf2_ = new DockTaskFactory( *old_instance.tf2_ );
}

void DockingHighRes::init( DockJumps const movable_jumps )
{
	moves::Mover::type( "DockingHighRes" );
	partners_ = "_";
	sc_min_ = false;
	rt_min_ = false;
	ignore_default_task_ = false; //needs to be false by default

	movable_jumps_ = movable_jumps;
	tf_ = NULL; //might need this
	tf2_ = new DockTaskFactory();
}

//destructor
DockingHighRes::~DockingHighRes() {}

void
DockingHighRes::set_task_factory( core::pack::task::TaskFactoryCOP tf )
{
	tf_ = new core::pack::task::TaskFactory( *tf );
}

void  DockingHighRes::set_ignore_default_task( bool ignore_default_task )
{
	ignore_default_task_ = ignore_default_task;
}

bool  DockingHighRes::ignore_default_task()
{
	return ignore_default_task_;
}

protocols::docking::DockTaskFactoryOP DockingHighRes::tf2()    //JQX: change COP to OP
{
	return tf2_;
}

void DockingHighRes::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
}

void DockingHighRes::set_scorefxn_pack( core::scoring::ScoreFunctionOP scorefxn_pack )
{
	scorefxn_pack_ = scorefxn_pack;
}

void DockingHighRes::set_interface_definition_task_operation( protocols::toolbox::task_operations::InterfaceTaskOperationOP interface_definition )
{
    tf2_->set_interface_definition_task_operation( interface_definition );
}

void DockingHighRes::set_additional_task_operarations( utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations )
{
    tf2_->set_additional_task_operarations( additional_task_operations );
}

void DockingHighRes::add_additional_task_operaration( core::pack::task::operation::TaskOperationOP task_operation )
{
    tf2_->add_additional_task_operaration( task_operation );
}

utility::vector1< core::pack::task::operation::TaskOperationOP > DockingHighRes::get_additional_task_operarations()
{
    return tf2_->get_additional_task_operarations();
}

core::scoring::ScoreFunctionOP DockingHighRes::scorefxn() const
{
	return scorefxn_;
}

core::scoring::ScoreFunctionOP DockingHighRes::scorefxn_pack() const
{
	return scorefxn_pack_;
}

core::pack::task::TaskFactoryCOP DockingHighRes::task_factory()
{
	return tf_;
}


} // namespace docking
} // namespace protocols

