// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Sid Chaudhury (Carbon copy by Ingemar André)

// Unit headers
#include <protocols/symmetric_docking/SymSidechainMinMover.hh>

// Project headers
#include <protocols/scoring/Interface.hh>

#include <core/types.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.SymSidechainMinMover" );

namespace protocols {
namespace symmetric_docking {

//default constructor
SymSidechainMinMover::SymSidechainMinMover() : moves::Mover()
{
	Mover::type( "SymSidechainMin" );
	update_movemap_ = true;
	set_default_options();
}

//construtor with arguments
SymSidechainMinMover::SymSidechainMinMover( core::scoring::ScoreFunctionCOP scorefxn_in) : moves::Mover(), scorefxn_(std::move(scorefxn_in))
{
	Mover::type( "SymSidechainMin" );
	set_default_options();
}

//construtor with arguments
SymSidechainMinMover::SymSidechainMinMover(
	core::scoring::ScoreFunctionCOP scorefxn_in,
	core::kinematics::MoveMapOP movemap_in
) : moves::Mover(), scorefxn_(std::move(scorefxn_in)), movemap_(std::move(movemap_in)), update_movemap_(false)
{
	Mover::type( "SymSidechainMin" );
	set_default_options();
}

//constructor with arguments
SymSidechainMinMover::SymSidechainMinMover(
	core::scoring::ScoreFunctionCOP scorefxn_in,
	core::pack::task::PackerTaskOP task_in) : moves::Mover(), scorefxn_(std::move(scorefxn_in)), task_(std::move(task_in)), update_movemap_(true)
{
	Mover::type( "SymSidechainMin" );
	set_default_options();
}

//constructor with arguments
SymSidechainMinMover::SymSidechainMinMover(
	core::scoring::ScoreFunctionCOP scorefxn_in,
	core::pack::task::TaskFactoryOP tf_in) : moves::Mover(), scorefxn_(std::move(scorefxn_in)), tf_(std::move(tf_in)), update_movemap_(true)
{
	Mover::type( "SymSidechainMin" );
	set_default_options();
}

//destructor
SymSidechainMinMover::~SymSidechainMinMover() = default;

void SymSidechainMinMover::set_minmover( protocols::simple_moves::MinMoverOP minmover_in ){ minmover_ = minmover_in; }

//default options setup for SymSidechainMinMover
void SymSidechainMinMover::set_default_options()
{
	if ( update_movemap_ ) {
		movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );
		movemap_->set_chi( true );
	}
	scorefxn_ = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	minmover_ = protocols::simple_moves::MinMoverOP( new simple_moves::symmetry::SymMinMover(movemap_, scorefxn_, "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/  ) );
}

std::string
SymSidechainMinMover::get_name() const {
	return "SidechainMin";
}

//SymSidechainMinMover updates the movemap based on a packer task
void SymSidechainMinMover::update_movemap( core::pose::Pose & pose)
{
	movemap_->set_chi( true );
	if ( tf_ ) task_ = tf_->create_task_and_apply_taskoperations( pose );
	if ( task_ ) {
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			if ( !task_->nonconst_residue_task(i).being_packed() ) movemap_->set_chi(i, false);
		}
	}
}

//SymSidechainMinMover apply function
void SymSidechainMinMover::apply( core::pose::Pose & pose )
{
	//runtime_assert(pose.is_fullatom());
	if ( update_movemap_ ) update_movemap( pose );
	core::pose::symmetry::make_symmetric_movemap(pose, *movemap_);
	minmover_->apply( pose );
}

//constructor
SymInterfaceSidechainMinMover::SymInterfaceSidechainMinMover() : SymSidechainMinMover(), interface_dist_(8.0)
{
	Mover::type( "SymInterfaceSidechainMinMover" );
	set_default_options();
}

SymInterfaceSidechainMinMover::SymInterfaceSidechainMinMover(
	core::scoring::ScoreFunctionCOP scorefxn_in,
	core::Real interface_dist_in
) : SymSidechainMinMover(scorefxn_in), interface_dist_(interface_dist_in)
{
	Mover::type( "SymInterfaceSidechainMin" );
	set_default_options();
}

//destructor
SymInterfaceSidechainMinMover::~SymInterfaceSidechainMinMover() = default;

void SymInterfaceSidechainMinMover::set_interface_dist( core::Real interface_dist_in)
{
	interface_dist_ = interface_dist_in;
	interface_->distance(interface_dist_);
}

//default options setup for SymInterfaceSidechainMinMover
void
SymInterfaceSidechainMinMover::set_default_options()
{
	interface_ = protocols::scoring::InterfaceOP( new protocols::scoring::Interface( 1 ) );

}

//apply function for SymInterfaceSidechainMinMover
void
SymInterfaceSidechainMinMover::apply( core::pose::Pose & pose )
{
	//reset movemap
	movemap_->set_bb(false);
	movemap_->set_jump(false);
	movemap_->set_chi(false);

	core::pose::symmetry::make_symmetric_movemap(pose, *movemap_);

	//calculate interface
	interface_->calculate( pose );

	minmover_->apply(pose);

}

std::string
SymInterfaceSidechainMinMover::get_name() const {
	return "SymInterfaceSidechainMin";
}

} // docking
} // protocols
