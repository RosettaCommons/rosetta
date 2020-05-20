// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/DEERData.hh
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  These classes contain the base and derived types for various DEER data containers.
///      The DEERData parent class stores generic information. The DEERDistanceBounds type
///      stores a distance value of interest and evaluates as a harmonic function. The
///      DEERDistanceDistribution type store the data as a probability distribution and
///      tries to maximize overlap. The DEERDecayData type stores the raw data and tries
///      to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_DEERData_hh
#define INCLUDED_core_scoring_epr_deer_DEERData_hh

// Unit headers
#include <core/scoring/epr_deer/DEERData.fwd.hh>
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

/// Base class
class DEERData : public utility::VirtualBase {
public:

	/// @brief Print the simulated distance distribution.
	void
	print_histogram(
		std::map< Size, Real > const & simulated_histogram
	) const;

	/// @brief  Virtual scoring method. Returns value of zero for the base class.
	/// @detail Takes as input a simulated histogram and an option to set the private
	///     score value to the result (default: false. This saves some code)
	virtual
	Real
	get_score(
		std::map< Size, Real > const &,
		bool const & set_score = false
	);

	/// @brief Returns the residues involved in this data set.
	/// @detail Residues are saves with two parameters: the residue ID, and the label type
	///     Label type is set to "DEFAULT" by default (duh). Other options include
	///     DEFAULT_FAST and CUSTOM
	utility::vector1< std::pair< Size, std::string > > const &
	residues() const;

	/// @brief Returns the granularity of the distance distribution. Default is 2 bins per A
	Size const &
	bins_per_angstrom() const;

	/// @brief Returns the last computed score. Obtained using get_score() or manually set
	Real
	score() const;

	/// @brief Returns the relative weight. Default: 1. Can be lowered when the data is less trustworthy
	Real
	relative_weight() const;
	/*
	/// @brief Returns to condition ID to which this data belongs
	Size
	condition() const;
	*/
	/// @brief Sets residue for data set; info for each "residue" consists of the index and the spin label type
	void
	residues(
		utility::vector1< std::pair< Size, std::string > > const & val
	);

	/// @brief Set the number of bins per angstrom for the data set.
	void
	bins_per_angstrom(
		Size const & val
	);

	/// @brief Set the score of the data set
	void
	score(
		Real const & val
	);

	/// @brief Set the relative weight of the data set
	void
	relative_weight(
		Real const & val
	);
	/*
	/// @brief Returns to condition ID to which this data belongs
	void
	condition(
	Size const & val
	);
	*/
protected:

	/// @brief Computes average distance of distribution for local functions
	Real
	compute_avg_dist(
		std::map< Size, Real > const & simulated_histogram
	) const;

protected:

	utility::vector1< std::pair< Size, std::string > > residues_;
	Size bins_per_angstrom_ = 2;
	Real score_;
	Real rel_weight_ = 1.0;
	//Size condition_ = 1;

};

/// @brief  Derived class for storing the data as a bounded function.
/// @detail Contains the upper and lower bounds, as well as steepness.
///     If the average simulated distance falls within these bounds the score is zero
///     If the average falls outside, it is evaluated via the steepness
///     For example: ( ( lb - d ) / s ) ^2 OR ( ( d  - ub ) / s ) ^2
///     lb: Lower bound
///     ub: Upper bound
///     d: Avg distance
///     s: Steepness
class DEERDistanceBounds : public DEERData {
public:

	/// @brief Computes average distance and score based on bounds.
	Real
	get_score(
		std::map< Size, Real > const & simulated_histogram,
		bool const & set_score = false
	) override;

	/// @brief Returns the lower and upper distance bounds
	std::pair< Real, Real > const &
	bounds() const;

	/// @brief Returns the step / steepness
	Real const &
	step() const;

	/// @brief Sets the lower and upper bounds
	void
	bounds(
		Real const & lower,
		Real const & upper
	);

	/// @brief Sets the step
	void
	step(
		Real const & step
	);

private:

	std::pair< Real, Real > bounds_ = std::make_pair( 0.0, 0.0 );
	Real step_ = 1.0;

};

/// Derived class that stores the entire distance distribution
/// Score is evaluated using the cross-entropy of the simulated from the experimental
class DEERDistanceDistribution : public DEERData {
public:

	/// @brief  Computes cross-entropy of simulated distribution from experimental.
	/// @detail Cross-entropy corresponds to the negative log-likelihood that the
	///     experimental distribution could have given rise to the simulated. This
	///     allows boltzmann weighting and/or Bayesian statistical inference from score
	///     Note that although confidence bands can be received as input, they are
	///     not currently used for this purpose. If you know an information-theoretic
	///     approach to using them, please contact me. I would live to incorporate
	///     that information here.
	Real
	get_score(
		std::map< Size, Real > const & simulated_histogram,
		bool const & set_score = false
	) override;

	/// @brief Returns the lower bound/confidence band for the distance distribution
	std::map< Size, Real > const &
	lower_bound() const;

	/// @brief Returns the line of best fit. Used to calculate cross-entropy
	std::map< Size, Real > const &
	best_fit() const;

	/// @brief Returns the upper bound/confidence band for the distance distribution
	std::map< Size, Real > const &
	upper_bound() const;

	/// @brief Sets the lower bound/confidence band for the distance distribution
	void
	lower_bound(
		std::map< Size, Real > const & val
	);

	/// @brief Sets the line of best fit. Used to calculate cross-entropy
	void
	best_fit(
		std::map< Size, Real > const & val
	);

	/// @brief Sets the upper bound/confidence band for the distance distribution
	void
	upper_bound(
		std::map< Size, Real > const & val
	);

	/// @brief Adds data to map at a specific distance. Overwrites if distance is occupied
	void
	append(
		Size const & dist_bin,
		Real const & lower,
		Real const & upper
	);

private:

	std::map< Size, Real > lower_bound_ = {};
	std::map< Size, Real > best_fit_ = {};
	std::map< Size, Real > upper_bound_ = {};

};

/// @brief Struct for fitting DEER traces. Copies are stored in each DEERDecayData object
struct FittingInfo {
	Real last_slope_ = 0.0;
	Real last_mod_depth_ = 0.0;
	Real last_dim_ = 3.0;
	Real time_pts_sqd_ = 0.0;
	std::map< Size, std::map< Real, Real > > spin_map_; // distance outer, time inner
	std::pair< Real, Real > mod_depth_bounds_ = std::make_pair( 0.02, 0.75 );
	std::string bckg_ = "3D";
};

/// @brief Derived class that stores DEER decay data, either raw or background-corrected
class DEERDecayData : public DEERData {
public:

	/// @brief Computes sum-of-squares of experimental to simulated DEER trace
	/// @detail A DEER trace is calculated from the simulated histogram. For details,
	///     please read del Alamo et al Biophysical Journal 2019. The score is
	///     the average sum-of-squares residuals between the two. Note that there is
	///     a slight noise-dependence to this; noisier data will always return higher
	///     scores.
	Real
	get_score(
		std::map< Size, Real > const & simulated_histogram,
		bool const & set_score = false
	) override;

	/// @brief Returns a gaussian distribution with certain avg and stdev.
	/// @detail Used for generating distributions when the standard deviation is used
	///     as a fitting parameter.
	std::map< Size, Real >
	mod_distr(
		Real const & avg,
		Real const & stdev
	) const;

	/// @brief Returns the experimental DEER trace. Key is time point in microseconds.
	std::map< Real, Real > const &
	trace() const;

	/// @brief Returns the last value obtained when fitting background slope (stored in fitting info)
	Real const &
	k_fit() const;

	/// @brief Returns the last value obtained when fitting modulation depth (stored in fitting info)
	Real const &
	modulation_depth_fit() const;

	/// @brief Returns the time points squared, used to fit slope and avoids unnecessary calculation.
	Real const &
	time_points_sqd() const;

	/// @brief Returns the noise from the imaginary component, provided manually as an option
	Real const &
	noise() const;

	/// @brief Returns the spin value for a given time point and distance. Computed at the front end.
	Real const &
	spin_val(
		Size const & dist_bin,
		Real const & time_pt
	) const;

	/// @brief Returns the upper and lower bounds of the modulation depth
	std::pair< Real, Real > const &
	mod_depth_bounds() const;

	/// @brief Returns whether or not the standard deviation is used as a fitting parameter
	bool const &
	fit_stdev() const;

	/// @brief Returns whether the dimensionality is used as a fitting parameter - typical for membrane proteins
	std::string const &
	bckg() const;

	/// @brief Sets the experimental DEER trace
	void
	trace(
		std::map< Real, Real > const & trace
	);

	/// @brief Appends time point and signal values for the experimental DEER trace
	void
	append_trace_data(
		Real const & time,
		Real const & signal
	);

	/// @brief Sets the value obtained when fitting background slope (stored in fitting info)
	void
	k_fit(
		Real const & k
	);

	/// @brief Sets the value obtained when fitting modulation depth (stored in fitting info)
	void
	modulation_depth_fit(
		Real const & modulation_depth
	);

	/// @brief Sets the value for time points squared, used to fit slope and avoids unnecessary calculation.
	void
	time_pts_sqd(
		Real const & time_pts_sqd
	);

	/// @brief Sets the noise from the imaginary component, provided manually as an option
	void
	noise(
		Real const & noise
	);

	/// @brief Sets the "spin map" that contains trace values at all time points and distances
	void
	spin_map(
		std::map< Size, std::map< Real, Real > > const & spin_map
	);

	/// @brief Sets lower and upper bounds for modulation depth
	void
	mod_depth_bounds(
		Real const & lower,
		Real const & upper
	);

	/// @brief Sets whether the standard deviation is a parameter when fitting DEER traces
	void
	fit_stdev(
		bool const & fit_stdev
	);

	/// @brief Sets what type of background is used for fitting: options are 3D, NON-3D, NONE
	void
	bckg(
		std::string const & bckg
	);

	/// @brief Appends data to experimental decay trace and values to relevant variables
	void
	append_trace_data_and_calculate(
		Real const & time,
		Real const & signal,
		Size const & max_bin
	);

	/// @brief Calculates the DEER trace given a time point (microseconds) and distance (A)
	Real
	calculate_decay_2_spin(
		Real const & t,
		Real const & r
	) const;

private:

	// TRACE INFO
	std::map< Real, Real > trace_;
	FittingInfo fit_info_;

	// DATA FOR PROCESSING
	Real noise_ = 1.0;
	bool fit_stdev_ = false;
	utility::vector0< Real > theta_terms_;

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
