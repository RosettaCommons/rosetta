// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/SymMinimizerMap.fwd.hh
/// @brief  MinimizerMap for symmetric minimization forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_optimization_symmetry_SymMinimizerMap_fwd_hh
#define INCLUDED_core_optimization_symmetry_SymMinimizerMap_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace optimization {
namespace symmetry {

class SymMinimizerMap;

typedef utility::pointer::shared_ptr< SymMinimizerMap > SymMinimizerMapOP;
typedef utility::pointer::shared_ptr< SymMinimizerMap const > SymMinimizerMapCOP;

} // symmetry
} // namespace optimization
} // namespace core

#endif
