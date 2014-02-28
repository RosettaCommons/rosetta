// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/pymol_chains.cc
/// @brief  method definitions for a couple PyMOL helper functions
/// @author Labonte
/// @note   I am simply moving the definitions out of the header file where they used to reside.

// Unit header
#include <basic/pymol_chains.hh>

// Numeric header
#include <numeric/types.hh>

// C++ header
#include <string>

namespace basic {

numeric::Size
get_pymol_num_unique_ids()
{
	return pymol_chains.length();
}

char
get_pymol_chain( numeric::Size i )
{
	return pymol_chains[ ( i - 1 ) % pymol_chains.size() ];
}

numeric::Size
get_pymol_chain_index_1( char c )
{
	numeric::Size i = pymol_chains.find( c );
	if ( i == pymol_chains.size() ) return 0;
	return i + 1;
}

}  // namespace basic
