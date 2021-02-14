// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERData.hh
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  These classes contain the base and derived types for various DEER data containers.
///      The DEERData parent class stores generic information. The DEERDistanceBounds type
///      stores a distance value of interest and evaluates as a harmonic function. The
///      DEERDistanceDistribution type store the data as a probability distribution and
///      tries to maximize overlap. The DEERDecayData type stores the raw data and tries
///      to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERData_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERData_hh

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.fwd.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief Alias to save space when defining residues and spin labels
using PairSizeString = std::pair< Size, std::string >;

/// Base class
class DEERData : public utility::VirtualBase {
public:

	/// @brief  Print the simulated distance distribution.
	/// @param  sim_histr: Simulated distance distribution
	/// @param  pose_name: Name of the pose
	void
	print_histogram(
		std::map< Size, Real > const & sim_histr,
		std::string const & pose_name = "NO NAME GIVEN"
	) const;

	/// @brief  Function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @param  set_score: Boolean to save score or not
	/// @return Freshly computed score
	Real
	get_score(
		std::map< Size, Real > const & sim_histr,
		bool const & set_score
	);

	/// @brief  Convolute a distribution with a Gaussian of a specific width
	/// @param  distr: DEER Distribution
	/// @param  std: Width of Gaussian (mean=0)
	/// @return New distribution
	std::map< Size, Real >
	convolute(
		std::map< Size, Real > const & distr,
		Real const & std
	) const;

	/// @brief  Virtual function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @return Freshly computed score
	virtual
	Real
	get_score(
		std::map< Size, Real > const & sim_histr
	) = 0;

	/// @brief  Returns the residues involved in this data set.
	/// @return Vector of residues (ID and spin label type)
	/// @detail Residues are saved with two parameters: the residue ID, and
	///  the label type. Label type is set to "DEFAULT" by default. Other
	///  options include DEFAULT_FAST and CUSTOM
	utility::vector1< PairSizeString > const &
	residues() const;

	/// @brief  Returns bins per angstrom for distribution (default: 2)
	/// @return Bins per angstrom
	Size const &
	bins_per_a() const;

	/// @brief  Returns the last computed score
	/// @return Score (0.0 if never set)
	Real
	score() const;

	/// @brief Returns the standard deviation of the distributions to generate
	/// @param  Standard deviation (in angstroms)
	/// @detail Function has a failsafe to avoid returning a nonzero value
	Real
	stdev() const;

	/// @brief Sets residue for data set
	/// @param residues: Vector of residues (index and spin label type)
	void
	residues(
		utility::vector1< PairSizeString > const & val
	);

	/// @brief Set the number of bins per angstrom for the data set.
	/// @param bins_per_a: Bins per angstrom
	void
	bins_per_a(
		Size const & val
	);

	/// @brief Set the score of the data set
	/// @param val: Score to save
	void
	score(
		Real const & val
	);

	/// @brief Set the standard deviation of the distributions to generate
	/// @param  val: Set the standard deviation to this value
	void
	stdev(
		Real const & val
	);

	/// @brief  Returns the map of distance values used for custom distributions
	/// @return Map of values (indeces to distances in Angstroms)
	/// @detail Only works when bins_per_a set to zero!
	std::map< Size, Real >
	dist_map() const;

	/// @brief  Append distance ID to custom distance map
	/// @param  dist_id: Unique distance ID used in distance map
	/// @param  dist: Distance value in angstroms
	/// @detail Only works when bins_per_a set to zero!
	void
	append_dist_id(
		Size dist_id,
		Real dist
	);

	/// @brief  Computes average distance of distribution for local functions
	/// @param  histogram: Distance distribution
	/// @return Vector with two items: Average and Standard Deviation
	utility::vector1< Real >
	avg_stdev(
		std::map< Size, Real > const & histogram
	);

protected:

	/// @brief Bins per angstrom for distributions stored and scored
	Size bins_per_a_ = 2;

private:

	/// @brief Residues used in this dataset
	utility::vector1< PairSizeString > residues_;

	/// @brief Last computed score
	Real score_ = 0.0;

	/// @brief Map if distance indeces to distance values (NOT amplitudes!)
	std::map< Size, Real > dist_map_ = {};

	/// @brief Standard deviation of distribution
	Real stdev_ = 1.0;

};

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
