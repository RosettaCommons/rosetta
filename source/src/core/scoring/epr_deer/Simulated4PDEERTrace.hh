// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/Simulated4PDEERTrace.cc
/// @brief  Class that simulates the observable signal resulting from electron dipolar coupling
/// @details The dipolar coupling signal between two electrons can be simulated from a distance
///       distribution. This class does that simulation by pulling values from the datacache
///       container objects (DEERData) and effectively storing a vector with values of interest.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_Simulated4PDEERTrace_hh
#define INCLUDED_core_scoring_epr_deer_Simulated4PDEERTrace_hh

#include <core/scoring/epr_deer/Simulated4PDEERTrace.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>

#include <map>
#include <tuple>

namespace core {
namespace scoring {
namespace epr_deer {

class Simulated4PDEERTrace {

public:

	/// @brief Default constructor
	Simulated4PDEERTrace() = default;

	/// @brief Destructor
	~Simulated4PDEERTrace() = default;

	/// @brief Default copy constructor
	/// @param rhs: To copy
	Simulated4PDEERTrace(
		Simulated4PDEERTrace const & rhs
	);

	/// @brief Constructor that instantiates everything
	/// @param deer_trace: The final DEER trace
	/// @param deer_trace_intra: The intramolecular DEER trace
	/// @param time_pts: Time points for the DEER trace
	/// @param depth: Modulation depth
	/// @param slope: Background slope
	/// @param dim: Background dimensionality
	Simulated4PDEERTrace(
		utility::vector1< Real > const & deer_trace,
		utility::vector1< Real > const & deer_trace_intra,
		utility::vector1< Real > const & time_pts,
		Real const & depth,
		Real const & slope,
		Real const & dim
	);

	/// @brief Copy assignment operator
	/// @param rhs: Trace object to copy
	/// @return Copy
	Simulated4PDEERTrace &
	operator=(
		Simulated4PDEERTrace const & rhs
	);

	std::tuple< Real, Real, Real >
	operator[](
		Size const & i
	) const;

	/// @brief DEER trace getter
	/// @return The final DEER trace
	utility::vector1< Real >
	deer_trace() const;

	/// @brief Intramolecular DEER trace getter
	/// @return The intramolecular DEER trace
	utility::vector1< Real >
	deer_trace_intra() const;

	/// @brief Time point getter
	/// @return Time points for the DEER trace
	utility::vector1< Real >
	time_pts() const;

	/// @brief Modulation depth getter
	/// @return Modulation depth
	Real
	depth() const;

	/// @brief Background slope getter
	/// @return Background slope
	Real
	slope() const;

	/// @brief Background dimensionality getter
	/// @return Background dimensionality
	Real
	dim() const;

	/// @brief  Return size of DEER trace
	/// @return Size of DEER trace
	Size
	size() const;

	/// @brief Final DEER trace setter
	/// @param deer_trace: The final DEER trace
	void
	deer_trace(
		utility::vector1< Real > const & val
	);

	/// @brief Intramolecular DEER trace setter
	/// @param deer_trace_intra: The intramolecular DEER trace
	void
	deer_trace_intra(
		utility::vector1< Real > const & val
	);

	/// @brief Time point vector setter
	/// @param time_pts: Time points for the DEER trace
	void
	time_pts(
		utility::vector1< Real > const & val
	);

	/// @brief Modulation depth setter
	/// @param depth: Modulation depth
	void
	depth(
		Real const & val
	);

	/// @brief Background slope setter
	/// @param slope: Background slope
	void
	slope(
		Real const & val
	);

	/// @brief Dimensionality setter
	/// @param dim: Background dimensionality
	void
	dim(
		Real const & val
	);

private:

	/// @brief Final DEER trace
	utility::vector1< Real > deer_trace_;

	/// @brief Intramolecular DEER trace (no background)
	utility::vector1< Real > deer_trace_intra_;

	/// @brief Time points for DEER traces
	utility::vector1< Real > time_pts_;

	/// @brief Modulation depth in final DEER trace
	Real depth_;

	/// @brief Background slope of final DEER trace
	Real slope_;

	/// @brief Background dimensionality of final DEER trace
	Real dim_;

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_epr_deer_Simulated4PDEERTrace_hh
