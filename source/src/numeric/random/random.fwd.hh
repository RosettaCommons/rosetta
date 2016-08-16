// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/random.fwd.hh
/// @brief  Random number generator system
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
///
/// @remarks
///  @li -


#ifndef INCLUDED_numeric_random_random_fwd_hh
#define INCLUDED_numeric_random_random_fwd_hh

#include <platform/types.hh>

// C++ headers
#include <cstddef>

namespace numeric {
namespace random {


typedef platform::Size Size;

double uniform(void);
double gaussian(void);

// Forward
class RandomGenerator;
class uniform_RG;


} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_FWD_HH
