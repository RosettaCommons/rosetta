// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/epr_deer/DEERData.hh
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  These classes contain the base and derived types for various DEER data containers.
///      The DEERData parent class stores generic information. The DEERDistanceBounds type
///      stores a distance value of interest and evaluates as a harmonic function. The
///      DEERDistanceDistribution type store the data as a probability distribution and
///      tries to maximize overlap. The DEERDecayData type stores the raw data and tries
///      to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_DEERData_fwd_hh
#define INCLUDED_core_scoring_epr_deer_DEERData_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace epr_deer {

class DEERData;

typedef utility::pointer::shared_ptr< DEERData > DEERDataOP;
typedef utility::pointer::shared_ptr< DEERData const > DEERDataCOP;

class DEERDistanceBounds;

typedef utility::pointer::shared_ptr< DEERDistanceBounds > DEERDistanceBoundsOP;
typedef utility::pointer::shared_ptr< DEERDistanceBounds const > DEERDistanceBoundsCOP;

class DEERDistanceDistribution;

typedef utility::pointer::shared_ptr< DEERDistanceDistribution > DEERDistanceDistributionOP;
typedef utility::pointer::shared_ptr< DEERDistanceDistribution const > DEERDistanceDistributionCOP;

class DEERDecayData;

typedef utility::pointer::shared_ptr< DEERDecayData > DEERDecayDataOP;
typedef utility::pointer::shared_ptr< DEERDecayData const > DEERDecayDataCOP;

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
