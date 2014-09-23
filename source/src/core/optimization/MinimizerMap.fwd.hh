// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/MinimizerMap.fwd.hh
/// @brief  forward declaration of MinimizerMap class
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_MinimizerMap_fwd_hh
#define INCLUDED_core_optimization_MinimizerMap_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace optimization {

class MinimizerMap;

typedef utility::pointer::shared_ptr< MinimizerMap > MinimizerMapOP;
typedef utility::pointer::shared_ptr< MinimizerMap const > MinimizerMapCOP;

} // namespace
} // namespace core


#endif
