// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ResidueTorsionRestraints.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_ResidueTorsionRestraints_hh
#define INCLUDED_protocols_ligand_docking_ResidueTorsionRestraints_hh

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.hh> ///TODO .fwd fails (why?)
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <set>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

///@brief Manages harmonic restraints on torsions, so they can be turned off for packing.
///
///@details Restraints are created when object is created, so they start off enabled.
/// I had to change from a PoseOP in the constructor to Pose references in enable/disable
/// because Movers only get Pose references, not PoseOPs.
/// Do not try to use one of these with multiple different poses, or surely the C++ gods will smite thee.
///
class ResidueTorsionRestraints : public utility::pointer::ReferenceCount
{
public:

	///@brief Establishes initial constraints set -- constraints start off enabled.
	ResidueTorsionRestraints(
		core::pose::Pose & pose,
		core::Size resid,
		core::Real stddev_degrees
	);
	virtual ~ResidueTorsionRestraints() {}

	///@brief Constrain residue torsions for specified pose.
	virtual void enable( core::pose::Pose & pose );

	///@brief Remove residue torsions constraints added by this object (if any).
	virtual void disable( core::pose::Pose & pose );

	bool operator==(const ResidueTorsionRestraints &other);

private:

	///@brief Shared logic for creating torsional constraints
	virtual void setup_constraints(core::pose::Pose & pose);

	///@brief Shared logic; returns old_constraints without my_constraints_
	virtual
	core::scoring::constraints::ConstraintSetOP
	without_my_constraints(
		core::scoring::constraints::ConstraintSetCOP old_constraints,
		std::set< core::scoring::constraints::ConstraintCOP > & removed_constraints
	);

	core::Size resid_;
	core::Real stddev_degrees_;
	/// Constraints that were created and added by this object.
	/// We compare by to_string() names rather than object identity
	/// because copying a Pose creates a deep copy of its constraints.
	std::set< std::string > my_constraints_;
	// Values from the last time a pose was disabled():
	utility::vector1< core::Real > old_chi_;
	std::set< core::scoring::constraints::ConstraintCOP > old_constraints_;

}; // ResidueTorsionRestraints


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_ResidueTorsionRestraints_hh
