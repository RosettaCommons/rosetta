// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/Simulated4PDEERTrace.cc
/// @brief  Class containing simulated DEER dipolar coupling decay data
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <map>
#include <tuple>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.Simulated4PDEERTrace" );

/// @brief Default copy constructor
/// @param rhs: To copy
Simulated4PDEERTrace::Simulated4PDEERTrace(
	Simulated4PDEERTrace const & rhs
) : deer_trace_( rhs.deer_trace_ ),
	deer_trace_intra_( rhs.deer_trace_intra_ ),
	time_pts_( rhs.time_pts_ ),
	depth_( rhs.depth_ ),
	slope_( rhs.slope_ ),
	dim_( rhs.dim_ )
{
	if ( deer_trace_.size() != deer_trace_intra_.size()
			|| deer_trace_.size() != time_pts_.size()
			|| deer_trace_intra_.size() != time_pts_.size()
			) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Sizes of DEER data and time points vectors do not match!" );
	}
}

/// @brief Constructor that instantiates everything
/// @param deer_trace: The final DEER trace
/// @param deer_trace_intra: The intramolecular DEER trace
/// @param time_pts: Time points for the DEER trace
/// @param depth: Modulation depth
/// @param slope: Background slope
/// @param dim: Background dimensionality
Simulated4PDEERTrace::Simulated4PDEERTrace(
	utility::vector1< Real > const & deer_trace,
	utility::vector1< Real > const & deer_trace_intra,
	utility::vector1< Real > const & time_pts,
	Real const & depth,
	Real const & slope,
	Real const & dim
) :
	// Assign variables
	deer_trace_( deer_trace ),
	deer_trace_intra_( deer_trace_intra ),
	time_pts_( time_pts ),
	depth_( depth ),
	slope_( slope ),
	dim_( dim )
{

	// Check that these line up
	if ( deer_trace_.size() != deer_trace_intra_.size()
			|| deer_trace_.size() != time_pts_.size()
			|| deer_trace_intra_.size() != time_pts_.size()
			) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Sizes of DEER data and time points vectors do not match!" );
	}
}

/// @brief Copy assignment operator
/// @param rhs: Trace object to copy
/// @return Copy
Simulated4PDEERTrace &
Simulated4PDEERTrace::operator=(
	Simulated4PDEERTrace const & rhs
) {
	deer_trace( rhs.deer_trace() );
	deer_trace_intra( rhs.deer_trace_intra() );
	time_pts( rhs.time_pts() );
	depth( rhs.depth() );
	slope( rhs.slope() );
	dim( rhs.dim() );
	return *this;
}

std::tuple< Real, Real, Real >
Simulated4PDEERTrace::operator[](
	Size const & i
) const {

	// Error in case the size is not large enough
	if ( i > time_pts_.size() ) {
		throw CREATE_EXCEPTION( utility::excn::KeyError,
			"Desired time point does not exist!" );
	} else {
		return std::make_tuple(
			time_pts_[ i ], deer_trace_intra_[ i ], deer_trace_[ i ] );
	}
}

/// @brief DEER trace getter
/// @return The final DEER trace
utility::vector1< Real >
Simulated4PDEERTrace::deer_trace() const {
	return deer_trace_;
}

/// @brief Intramolecular DEER trace getter
/// @return The intramolecular DEER trace
utility::vector1< Real >
Simulated4PDEERTrace::deer_trace_intra() const {
	return deer_trace_intra_;
}

/// @brief Time point getter
/// @return Time points for the DEER trace
utility::vector1< Real >
Simulated4PDEERTrace::time_pts() const {
	return time_pts_;
}

/// @brief Modulation depth getter
/// @return Modulation depth
Real
Simulated4PDEERTrace::depth() const {
	return depth_;
}

/// @brief Background slope getter
/// @return Background slope
Real
Simulated4PDEERTrace::slope() const {
	return slope_;
}

/// @brief Background dimensionality getter
/// @return Background dimensionality
Real
Simulated4PDEERTrace::dim() const {
	return dim_;
}

/// @brief  Return size of DEER trace
/// @return Size of DEER trace
Size
Simulated4PDEERTrace::size() const {
	return deer_trace_.size();
}

/// @brief Final DEER trace setter
/// @param deer_trace: The final DEER trace
void
Simulated4PDEERTrace::deer_trace(
	utility::vector1< Real > const & val
) {
	deer_trace_ = val;
}

/// @brief Intramolecular DEER trace setter
/// @param deer_trace_intra: The intramolecular DEER trace
void
Simulated4PDEERTrace::deer_trace_intra(
	utility::vector1< Real > const & val
) {
	deer_trace_intra_ = val;
}

/// @brief Time point vector setter
/// @param time_pts: Time points for the DEER trace
void
Simulated4PDEERTrace::time_pts(
	utility::vector1< Real > const & val
) {
	time_pts_ = val;
}

/// @brief Modulation depth setter
/// @param depth: Modulation depth
void
Simulated4PDEERTrace::depth(
	Real const & val
) {
	depth_ = val;
}

/// @brief Background slope setter
/// @param slope: Background slope
void
Simulated4PDEERTrace::slope(
	Real const & val
) {
	slope_ = val;
}

/// @brief Dimensionality setter
/// @param dim: Background dimensionality
void
Simulated4PDEERTrace::dim(
	Real const & val
) {
	dim_ = val;
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
