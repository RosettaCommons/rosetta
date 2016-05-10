// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/OptionCollection.fwd.hh
/// @brief  utility::options::OptionCollection forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_OptionCollection_fwd_hh
#define INCLUDED_utility_options_OptionCollection_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace options {


// Forward
class OptionCollection;

typedef utility::pointer::shared_ptr< OptionCollection > OptionCollectionOP;
typedef utility::pointer::shared_ptr< OptionCollection const > OptionCollectionCOP;

/// @brief Option types enumeration
enum OptionTypes {
	UNKNOWN_OPTION,
	BOOLEAN_OPTION,
	INTEGER_OPTION,
	REAL_OPTION,
	STRING_OPTION,
	FILE_OPTION,
	PATH_OPTION,
	ANY_OPTION,
	BOOLEAN_VECTOR_OPTION,
	INTEGER_VECTOR_OPTION,
	REAL_VECTOR_OPTION,
	RESIDUE_CHAIN_VECTOR_OPTION,
	STRING_VECTOR_OPTION,
	FILE_VECTOR_OPTION,
	PATH_VECTOR_OPTION,
	ANY_VECTOR_OPTION
};


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_OptionCollection_FWD_HH
