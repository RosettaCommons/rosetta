// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SymRestrictTaskForDocking.cc
/// @brief When passed to a PackerTask, pack/design is limited to the protein interface
/// @author Ingemar Andre

#include <protocols/symmetric_docking/SymRestrictTaskForDocking.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
// Symmetry


#include <protocols/scoring/Interface.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

using namespace core;
using namespace scoring;
using namespace pack;

SymRestrictTaskForDocking::SymRestrictTaskForDocking()
: scorefxn_( /* 0 */ ),
	include_current_( true ),
	distance_( 0 )
{}

SymRestrictTaskForDocking::SymRestrictTaskForDocking(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool include_current,
	core::Real distance
) : scorefxn_( scorefxn ),
	include_current_( include_current ),
	distance_( distance )
{}

SymRestrictTaskForDocking::~SymRestrictTaskForDocking(){}


task::operation::TaskOperationOP SymRestrictTaskForDocking::clone() const
{
	return task::operation::TaskOperationOP( new SymRestrictTaskForDocking( *this ) );
}

void
SymRestrictTaskForDocking::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const
{
	task.initialize_from_command_line().restrict_to_repacking().or_include_current( include_current_ );

	assert( scorefxn_ != 0 );
	// (existing comment) /// why is this still necessary???
	// (*scorefxn_)(pose);
	// scorefxn_->accumulate_residue_total_energies( pose );

	runtime_assert( scorefxn_ != 0 );
	runtime_assert( distance_ );

	protocols::scoring::Interface interface( 1 );
	interface.distance( distance_ );
	interface.calculate( pose );

	core::pack::task::PackerTaskOP task_copy = task.get_self_ptr();
	interface.set_symmetric_pack( pose, task_copy );
}

} // namespace symmetric_docking
} // namespace protocols

