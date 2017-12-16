// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/pointer/std/memory.hh
/// @brief  header for std make_shared
/// @author Sergey Lyskov


#ifndef INCLUDED_utility_pointer_std_memory_hh
#define INCLUDED_utility_pointer_std_memory_hh

#ifdef PTR_STD

#include <memory>

namespace utility {
namespace pointer {

using std::make_shared;

} // namespace pointer
} // namespace utility

#endif // PTR_STD

#endif // INCLUDED_utility_pointer_std_memory_hh
