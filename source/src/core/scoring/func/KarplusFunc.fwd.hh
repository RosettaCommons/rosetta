// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/func/CircularHarmonicFunc.fwd.hh
/// @brief forward declaration for CircularHarmonicFunc class
/// @author Nikolas Sgourakis


#ifndef INCLUDED_core_scoring_func_KarplusFunc_fwd_hh
#define INCLUDED_core_scoring_func_KarplusFunc_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace func {

class KarplusFunc;
typedef utility::pointer::shared_ptr< KarplusFunc > KarplusFuncOP;
typedef utility::pointer::shared_ptr< KarplusFunc const > KarplusFuncCOP;

}
}
}

#endif
