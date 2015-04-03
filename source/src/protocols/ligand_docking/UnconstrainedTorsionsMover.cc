// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/UnconstrainedTorsionsMover.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <protocols/ligand_docking/MinimizeLigand.hh>


// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {


UnconstrainedTorsionsMover::UnconstrainedTorsionsMover(
	protocols::moves::MoverOP child_mover,
	Restraints restraints
):
	Mover(),
	child_mover_(child_mover),
	restraints_(restraints)
{
}

UnconstrainedTorsionsMover::UnconstrainedTorsionsMover(
	protocols::moves::MoverOP child_mover,
	std::set<ResidueTorsionRestraintsOP> restraints
):
	Mover(),
	child_mover_(child_mover)
{
	BOOST_FOREACH(ResidueTorsionRestraintsOP restraint, restraints){
		restraints_.push_back(restraint);
	}
}

UnconstrainedTorsionsMover::UnconstrainedTorsionsMover(
	protocols::moves::MoverOP child_mover,
	MinimizeLigandOPs minimize_ligands
):
	Mover(),
	child_mover_(child_mover)
{
	BOOST_FOREACH(MinimizeLigandOP minimize_ligand, minimize_ligands){
		restraints_.insert( restraints_.end(), minimize_ligand->begin(), minimize_ligand->end() );
	}
}

void UnconstrainedTorsionsMover::apply( core::pose::Pose & pose )
{
	BOOST_FOREACH(ResidueTorsionRestraintsOP restraint, restraints_){
		restraint->disable( pose );
	}
	child_mover_->apply(pose);

	BOOST_FOREACH(ResidueTorsionRestraintsOP restraint, restraints_){
		restraint->enable( pose );
	}
}

std::string
UnconstrainedTorsionsMover::get_name() const {
	return "UnconstrainedTorsionsMover";
}


} // namespace ligand_docking
} // namespace protocols
