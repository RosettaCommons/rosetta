// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/util.hh
/// @brief Utilities for modifying and utilizing Residues and other core::chemical classes.

#ifndef INCLUDED_core_chemical_util_HH
#define INCLUDED_core_chemical_util_HH

// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ header
#include <string>
#include <map>


namespace core {
namespace chemical {

/// @brief  Return a constant access pointer to the ResidueTypeSet specified by the command-line options.
core::chemical::ResidueTypeSetCAP rsd_set_from_cmd_line();

/// @brief  Add additional parameter files not present in <atom-set-name>/extras.txt.
void add_atom_type_set_parameters_from_command_line(
	std::string const & atom_type_set_tag,
	AtomTypeSet & atom_type_set );

/// @brief  Modify atom_type properties from the command line.
void modify_atom_properties_from_command_line(
	std::string const & atom_type_set_tag,
	AtomTypeSet & atom_type_set );


/// @brief Return a string representing the internal coordinates tree of this ResidueType.
std::string formatted_icoord_tree( core::chemical::ResidueType const & restype );

/// @brief Utility to examine chi output.
void print_chis( std::ostream & out, ResidueType const & res );


/// @brief Replaces the deprecated "_p:" linker connecting ResidueType base names with their patch names with ":".
std::string fixup_patches( std::string string_in );


/// @brief  Are these two residues patched in exactly the same way, ignoring any VariantTypes in the list of exceptions?
bool variants_match_with_exceptions(
	ResidueType const & res1,
	ResidueType const & res2,
	utility::vector1< VariantType > list_of_variants_to_ignore );


utility::vector1< VariantType > pH_mode_exceptions();

/// @brief  Are these two residues patched in exactly the same way?
bool variants_match( ResidueType const & res1, ResidueType const & res2 );

/// @brief  Similar to variants_match(), but allows different adduct-modified states.
bool nonadduct_variants_match( ResidueType const & res1, ResidueType const & res2 );

/// @brief look for best match to atom_names
ResidueTypeCOP find_best_match( ResidueTypeCOPs const & rsd_type_list,
	utility::vector1< std::string > const & atom_names,
	bool const ignore_atom_named_H = false );

/// @brief Fang-Chieh Chou 8/10/2012. Use larger LJ_WDEPTH for protons to avoid clashes in RNA
void enlarge_h_lj_wdepth( utility::vector1< Real > & lj_wdepth, AtomTypeSet const & atom_type_set );

/// @brief Fang-Chieh Chou 8/10/2012. Use larger LJ_WDEPTH for protons to avoid clashes in RNA
void enlarge_h_lj_wdepth( AtomTypeSet & atom_type_set );

/// @brief Rhiju. O3', O4', O5' in nucleic acids are ethers -- turn them off as acceptors for H-bonds.
void unset_acceptor_ether_oxygens( AtomTypeSet & atom_type_set );

void detect_ld_chirality_from_polymer_residue(
	std::map< std::string, Vector > const & xyz,
	std::string const & name3,
	bool & is_d_aa,
	bool & is_l_aa );

/// @brief Return true if the two residues have the same number and name of heavy atoms.
bool heavy_atom_names_match( ResidueType const & first, ResidueType const & second );


/// @brief  Are these main-chain torsions also ring torsions?
bool is_mainchain_torsion_also_ring_torsion( ResidueType const & res_type, uint torsion_index );

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_util_HH
