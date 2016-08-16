// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/UnconstrainedTorsionsMover.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_UnconstrainedTorsionsMover_hh
#define INCLUDED_protocols_ligand_docking_UnconstrainedTorsionsMover_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh> ///TODO make this a .fwd
#include <protocols/ligand_docking/UnconstrainedTorsionsMover.fwd.hh>
#include <protocols/ligand_docking/MinimizeLigand.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

/// @brief Juggles torsional constraints with packing or rotamer trials.
///
/// @details Adds torsional constraints to the specified residue in the pose,
/// but then removes them before running the supplied mover.
/// They are then either restored (if no conformational change) or re-initialized.
/// The supplied child_mover is expected to be either a full repack or rotamer trials.
///
class UnconstrainedTorsionsMover : public protocols::moves::Mover
{
public:
	typedef utility::vector1< ResidueTorsionRestraintsOP > Restraints;

	UnconstrainedTorsionsMover(
		protocols::moves::MoverOP child_mover,
		Restraints restraints
	);

	UnconstrainedTorsionsMover(
		protocols::moves::MoverOP child_mover,
		std::set<ResidueTorsionRestraintsOP> restraints
	);

	UnconstrainedTorsionsMover(
		protocols::moves::MoverOP child_mover,
		MinimizeLigandOPs
	);


	virtual ~UnconstrainedTorsionsMover() {}

	/// @brief Removes its constraints, runs mover, restores constraints.
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	protocols::moves::MoverOP child_mover_;
	Restraints restraints_;

}; // UnconstrainedTorsionsMover


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_UnconstrainedTorsionsMover_hh
