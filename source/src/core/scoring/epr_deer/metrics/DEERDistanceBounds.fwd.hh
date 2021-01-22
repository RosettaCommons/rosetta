// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/epr_deer/metrics/DEERDistanceBounds.hh
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  This class contains upper and lower bounds for experimental
///    DEER data using a harmonic function.
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERDistanceBounds_fwd_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERDistanceBounds_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

class DEERDistanceBounds;

typedef utility::pointer::shared_ptr< DEERDistanceBounds > DEERDistanceBoundsOP;
typedef utility::pointer::shared_ptr< DEERDistanceBounds const > DEERDistanceBoundsCOP;

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
