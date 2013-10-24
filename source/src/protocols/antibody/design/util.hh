// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/util.hh
/// @brief Utility functions for antibody design namespace.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_UTIL_HH
#define INCLUDED_protocols_antibody_design_UTIL_HH

#include <protocols/antibody/design/AntibodyDesignEnum.hh>

#include <string>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace antibody{
namespace design{
using namespace protocols::antibody;
using namespace utility;

///@brief Returns (?,?,?) With question marks of length n to help create database query using IN operator
std::string
get_string_for_IN(core::Size const n);

///@brief Gets all possible graft permutations.
///@details all_permutations is a list of vectors corresponding to cdrs_to_design vector.  Essentially, each inner index describes a position in the cdr_set.
/// Indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
/// Example: <1, 0, 1, 1, 1, 1>.  This is a possible combination to try graft, the second CDR, H2 is not part of the combination.
void
get_all_graft_permutations(
	vector1<core::Size > & total_cdr_set,
	vector1<vector1< core::Size > > & all_permutations,
	vector1< core::Size >current_index,
	core::Size const cdr_num);

DesignTypeEnum
design_type_from_string(std::string const design_type);

// Undefined, commenting out to fix PyRosetta build  std::string design_type_from_enum(DesignTypeEnum const design_type);

} //design
} //antibody
} //protocols


#endif	//#ifndef INCLUDED_protocols/antibody/design_UTIL_HH
