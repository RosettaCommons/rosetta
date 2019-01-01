// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PRESingleSet.hh
/// @brief   class that stores and handles data for one single PRE dataset/experiment (i.e. for one spinlabel at on tagging site)
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pre_PRESingleSet_HH
#define INCLUDED_core_scoring_nmr_pre_PRESingleSet_HH

// Unit headers
#include <core/scoring/nmr/pre/PRESingleSet.fwd.hh>

// Package headers
#include <core/io/nmr/AtomSelection.fwd.hh>
#include <core/scoring/nmr/pre/PRESingle.fwd.hh>
#include <core/scoring/nmr/pre/PREMultiSet.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <string>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

class PRESingleSet {

	friend class PREMultiSet; // We have to access the PRESingle Vector through the class PREMultiSet
	// but I don't want to make that data accessible for any other class
public: // Methods

	/// @brief construct PRESingleSet from data file, use mutator methods to set experimental conditions
	PRESingleSet(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief construct PRESingleSet from data file and set computation, weighting and rate type
	PRESingleSet(
		std::string const & filename,
		pose::Pose const & pose,
		Real const weight,
		std::string rate = "R2",
		std::string single_pre_weigting = "CONST"
	);

	/// @brief copy constructor
	PRESingleSet(PRESingleSet const & other);

	/// @brief assignment operator
	PRESingleSet &
	operator=(PRESingleSet const & rhs);

	/// @brief destructor
	~PRESingleSet();

	/// @brief return gyromagnetic ratio of the nuclear spin in rad/(s*T) (dimension is 10^6)
	Real calc_gamma_I() const;

	/// @brief calculate nuclear spin frequency at given field strength in rad/s
	Real calc_omega_I() const;

	/// Getters
	std::string get_dataset_name() const { return dataset_name_; }
	Size get_number_pre() const { return number_pre_; }
	Real get_weight() const { return weight_; }
	Real get_scaling_factor() const { return scaling_factor_; }
	utility::vector1< PRESingle > const & get_pre_single_vec() const { return pre_single_vec_; }
	SINGLE_NMR_VALUE_WEIGHTING get_single_pre_weighting_scheme() const { return single_pre_weighting_scheme_; }
	PRE_RATE_TYPE get_pre_rate_type() const { return rate_type_; }
	std::string pre_rate_type_to_string() const { return rate_type_ == R2_PARA ? "R2" : "R1"; }
	Real get_field_strength() const { return field_strength_; }
	bool normalized_data() const { return normalized_data_; }

	// Setters
	void set_weight(Real weight) { weight_ = weight; }
	void set_field_strength(Real field) { runtime_assert(field > 0.0); field_strength_ = field; }
	void set_single_pre_weighting_scheme(std::string const & weighting_scheme);

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	PRESingleSet();

	/// @brief utility function to initialize PRESingleSet from filedata
	void
	init_from_filedata(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief register options
	void register_options();
	void init_from_cml();

private: // Data

	std::string dataset_name_;
	utility::vector1<PRESingle> pre_single_vec_;
	Real weight_;      // weight the experiment score by this factor
	Real scaling_factor_;    // scale PREs by experiment StdDev
	Size number_pre_;
	PRE_RATE_TYPE rate_type_;
	SINGLE_NMR_VALUE_WEIGHTING single_pre_weighting_scheme_;

	// Experimental conditions
	Real field_strength_;  // in MHz
	// Normalization by experiment StdDev
	bool normalized_data_;

};

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pre_PRESingleSet_HH
