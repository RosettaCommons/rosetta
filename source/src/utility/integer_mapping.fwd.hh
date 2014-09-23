// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/integer_mapping.hh
/// @brief  A set of useful classes to map between two enumerations.  So far, only a subset mapping is implemented.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_integer_mapping_fwd_hh
#define INCLUDED_utility_integer_mapping_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility {

class subset_mapping;
typedef pointer::shared_ptr< subset_mapping      > subset_mappingOP;
typedef pointer::shared_ptr< subset_mapping const> subset_mappingCOP;

}

#endif
