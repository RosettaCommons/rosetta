// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/legacy_types.hh
/// @brief  types that were in use for patching by mm_params in 2015.
/// @author Rhiju Das
#ifndef INCLUDED_core_chemical_legacy_types_hh
#define INCLUDED_core_chemical_legacy_types_hh

#include <utility/vector1.hh>
#include <core/chemical/VariantType.hh>

namespace core {
namespace chemical {

	/// @details only in use by name3_map_DO_NOT_USE and aa_map_DO_NOT_USE.
	/// list of variant_types that were "standard" in Rosetta in late 2015, and
	/// only these will be supported by those DO_NOT_USE legacy functions.
	/// Newer variant_ types should *not* be placed into this list, or apps that
	///  ask for *all* variants with an aa or name3 will get super-bloated.
	/// If those apps need new variant types, please make use of ResidueTypeFinder
	///  to request specific subsets of ResidueTypes as needed.
	utility::vector1< VariantType > const &
	variant_types_list_LEGACY();

} // namespace chemical
} // namespace core


#endif // INCLUDED_core_legacy_types_HH
