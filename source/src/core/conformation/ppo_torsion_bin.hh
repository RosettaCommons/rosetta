// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/ppo_torsion_bin.hh
/// @brief  An enumeration to represent a binning of the phi/psi/omega torsions
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_conformation_ppo_torsion_bin_HH
#define INCLUDED_core_conformation_ppo_torsion_bin_HH

// Unit headers
#include <core/conformation/ppo_torsion_bin.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector0.hh>

// C++ headers
#include <string>

namespace core {
namespace conformation {

/// @brief determine the torsion bin for a given phi/psi/omega combination, assuming that omega is 180 if not specified
/// @author Amelie Stein (amelie.stein@ucsf.edu)
/// @date Wed May  2 11:18:29 PDT 2012
ppo_torsion_bin
get_torsion_bin (core::Real phi, core::Real psi, core::Real omega = 180);

ppo_torsion_bin
remap_cis_omega_torsion_bins_to_trans( ppo_torsion_bin torbin );

bool
cis_omega_torsion_bin( ppo_torsion_bin torbin );

/// @brief convert a string of characters into a vector of the internally recognized ppo_torsion_bin enumeration
/// @throws utility::excn::EXCN_MsgException if any of the input characters in this string are invalid
torsion_bin_string
map_string_to_torsion_bin_string( std::string const & torstring );

/// @brief returns true if the input character represents a valid torsion bin
bool
char_valid_as_torsion_bin( char torbin );

/// @brief returns the torsion bin that the input character represents
ppo_torsion_bin
map_char_to_torsion_bin( char torbin );

/// @brief convert a torsion bin to a character s.t. that character can be converted back to a torsion bin
char
map_torsion_bin_to_char( ppo_torsion_bin torbin );

} // conformation
} // core

#endif
