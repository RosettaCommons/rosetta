// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/unti/tracer.fwd.hh
/// @brief  Uniform random number generator
/// @author Sergey Lyskov


#ifndef INCLUDED_numeric_random_uniform_fwd_hh
#define INCLUDED_numeric_random_uniform_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace numeric {
namespace random {

class uniform_RG;

typedef utility::pointer::shared_ptr< uniform_RG > uniform_RG_OP;
typedef utility::pointer::shared_ptr< uniform_RG const > uniform_RG_COP;

} // namespace random
} // namespace numeric

#endif // INCLUDED_numeric_random_uniform_fwd_hh
