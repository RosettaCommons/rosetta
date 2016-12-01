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
/// @author

#ifndef INCLUDED_basic_pymol_chains_hh
#define INCLUDED_basic_pymol_chains_hh

#include <string>
#include <numeric/types.hh>

namespace basic {

//fd removing '{' and '}' since they mess with pymol selection
static std::string const pymol_chains(
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz!@#$&.<>?]|-_\\~=%" );


numeric::Size get_pymol_num_unique_ids();

char get_pymol_chain( numeric::Size i );

numeric::Size get_pymol_chain_index_1( char c );

}

#endif  // INCLUDED_basic_pymol_chains_hh
