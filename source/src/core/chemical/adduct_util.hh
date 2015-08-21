// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Jim Havranek


#ifndef INCLUDED_core_chemical_adduct_util_hh
#define INCLUDED_core_chemical_adduct_util_hh

// ObjexxFCL Headers

// Unit headers
#include <core/chemical/ResidueType.fwd.hh>

// Project headers

// Utility headers

#include <utility/options/StringVectorOption.fwd.hh>
#include <map>


// C++ headers

namespace core {
namespace chemical {

typedef std::map< std::string, int > AdductMap;

/// @brief Convert input string to map of adducts->max usage
AdductMap
parse_adduct_string(
	utility::options::StringVectorOption & add_vec
);

/// @brief Make sure requested adducts exist in some residue
void
error_check_requested_adducts( AdductMap const & add_map,
	ResidueTypeCOPs const & rsd_types );

/// @brief Apply adducts to residue using a boolean mask
ResidueTypeOP apply_adducts_to_residue(
	ResidueType const & rsd,
	utility::vector1< bool > & add_mask
);

} // chemical
} // core


#endif
