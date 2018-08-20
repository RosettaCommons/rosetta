// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/atomic_depth/AtomicDepth.fwd.hh
/// @brief  AtomicDepth forward declaration
/// @author Brian Coventry

#ifndef INCLUDED_core_scoring_atomic_depth_AtomicDepth_fwd_hh
#define INCLUDED_core_scoring_atomic_depth_AtomicDepth_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace atomic_depth {

// Forward
class AtomicDepth;

typedef utility::pointer::shared_ptr< AtomicDepth > AtomicDepthOP;
typedef utility::pointer::shared_ptr< AtomicDepth const > AtomicDepthCOP;

} // namespace atomic_depth
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_atomic_depth_AtomicDepth_fwd_hh
