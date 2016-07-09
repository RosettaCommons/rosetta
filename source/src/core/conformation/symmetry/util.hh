// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/pose/symmetry/util.hh
/// @brief utility functions for handling with symmetric conformations
/// @author Ingemar Andre

#ifndef INCLUDED_core_conformation_symmetry_util_hh
#define INCLUDED_core_conformation_symmetry_util_hh


// Unit headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>

namespace core {
namespace conformation {
namespace symmetry {

bool
is_symmetric( conformation::Conformation const & conf );

bool
is_mirror_symmetric( conformation::Conformation const & conf );

bool
is_symmetric( conformation::symmetry::SymmetryInfo const & symminfo );

conformation::symmetry::SymmetricConformationOP
setup_symmetric_conformation(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata,
	std::map< int, char > conf2pdb_chain = std::map< int, char >()
);

kinematics::FoldTree
set_fold_tree_from_symm_data(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata,
	std::map< int, char > conf2pdb_chain = std::map< int, char >()
);

kinematics::FoldTree
replaced_symmetric_foldtree_with_new_monomer(
	kinematics::FoldTree symm_f,
	conformation::symmetry::SymmetryInfo symmetry_info,
	kinematics::FoldTree monomer_f
);

void
recenter(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata
);

// @details shift jump numbers in dof
void
shift_jump_numbers_in_dofs(
	conformation::Conformation & src_conformation,
	Size shift
);

kinematics::FoldTree
get_asymm_unit_fold_tree( core::conformation::Conformation const &conf );

/// @brief    Converts an asymmetric foldtree (f) with virtual root into a
///           symmetric foldtree compatible with symmetric conformation (conf)
/// @param    conf - A symmetric conformation
/// @param    f - An asymmetric foldtree. This foldtree MUST have a virtual root
/// @details  This function does not require the symm data
void
symmetrize_fold_tree( core::conformation::Conformation const &conf, kinematics::FoldTree &f );

void
set_asymm_unit_fold_tree( core::conformation::Conformation &p, kinematics::FoldTree const &f);

int
residue_center_of_mass(
	conformation::Conformation const & conformation,
	int const start,
	int const stop
);

int
return_nearest_residue(
	conformation::Conformation const & conformation,
	int const begin,
	int const end,
	core::Vector center
);

std::string
show_foldtree(
	core::conformation::symmetry::SymmetricConformation const & symm_conf,
	SymmData const & symdata,
	std::map<char,std::pair<Size,Size> > const & chain2range
);

/// @breif A very specific helper function that counts the number of mirror ops from root->each subunit
void
calculate_inverting_virtuals(
	core::kinematics::FoldTree const & ft,
	Conformation const & conf,
	SymmetryInfo const &symm_info,
	utility::vector1<bool> &mirrored_subs,
	utility::vector1< std::pair<bool,bool> > &mirrored_jumps );


} // symmetry
} // conformation
} // core


#endif
