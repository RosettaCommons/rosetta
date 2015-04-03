// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ligand_functions.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_ligand_functions_hh
#define INCLUDED_protocols_ligand_docking_ligand_functions_hh

//#include <protocols/ligand_docking/ligand_functions.fwd.hh>
//#include <utility/pointer/ReferenceCount.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {


/// @brief  Produce an ambiguous dihedral restraint for the specified chi angle,
/// assuming that the provided conformations represent (all) energetic minima.
core::scoring::constraints::ConstraintOP
torsion_constraints_from_rotamers(
	core::Size rsd_no,
	core::Size chino,
	utility::vector1< core::conformation::ResidueCOP > const & rsds,
	core::Real stddev_degrees
);

/// @brief  Produce an ambiguous dihedral restraint for the specified chi angle,
/// assuming that ResidueType.chi_rotamers() lists (all) energetic minima.
core::scoring::constraints::ConstraintOP
torsion_constraints_from_chi_rotamers(
	core::Size rsd_no,
	core::Size chino,
	core::chemical::ResidueType const & rsdtype
);

/// @brief Produce dihedral restraints for all chi angles in the specified
/// residue, from chi_rotamers() if available, and from the rotamer library otherwise.
void
get_ligand_torsion_constraints(
	core::pose::Pose & pose,
	core::Size rsd_no,
	core::Real stddev_degrees,
	utility::vector1< core::scoring::constraints::ConstraintOP > & csts_out,
	bool const constrain_all_torsions_equally
);

/// @brief Call get_ligand_torsion_constraints() for all non-polymer residues
/// and add the resulting constraints to the Pose.
void
constrain_ligand_torsions(
	core::pose::Pose & pose,
	core::Real stddev_degrees,
	bool constrain_all_torsions_equally = true
);

/// @brief simple function to scan the pose for all ligand residues
utility::vector1< core::Size >
get_ligand_seqpos(
	core::pose::Pose const & pose
);


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_ligand_functions_HH
