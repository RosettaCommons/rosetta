// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/util.hh
/// @brief definition of classes for iterations of docking/design.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_util_hh
#define INCLUDED_protocols_protein_interface_design_util_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility Headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {

core::kinematics::FoldTree
star_fold_tree( core::pose::Pose & pose );

/// @brief removes ALL coordinate constraints from a pose. returns the constraints that were removed
core::scoring::constraints::ConstraintCOPs remove_coordinate_constraints_from_pose( core::pose::Pose & pose );

// @brief returns ConstraintCOPs matching only the backbone_stub_constraints currently in the pose
core::scoring::constraints::ConstraintCOPs get_bbcsts( core::pose::Pose const & pose );

/// @brief evaluate backbone_stub_constraints for each residue in a chain and return a vector with the top n_return residue numbers by cst score
/// note that this function is NOT guaranteed to return n_return residues! It will return the best n<=n_return
utility::vector1< core::Size >
best_bbcst_residues( core::pose::Pose const & pose, core::Size const chain, core::Size const n_return );

void find_lowest_constraint_energy_residue( core::pose::Pose const & pose, core::Size const chain, core::Size & resi, core::Real & lowest_energy );

/// @brief utility function for stub_based_atom_tree. tries to find an optimal cutpoint in a pose given two different boundaries.
core::Size best_cutpoint( core::pose::Pose & pose, core::Size const prev_u, core::Size const prev_d, core::Size const u, core::Size const d );

/// @brief find nearest residue on target_chain to res
core::Size
find_nearest_residue( core::pose::Pose const & pose, core::Size const target_chain, core::Size const res, std::string const & atom="CA" );

/// @brief what is the optimal connection point for an atom tree, given a residue type (puts the connection point
/// at the beginning of a functional group
std::string optimal_connection_point( std::string const & residue_type );

core::kinematics::FoldTree
make_hotspot_foldtree( core::pose::Pose const & pose );

} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_util_HH*/

