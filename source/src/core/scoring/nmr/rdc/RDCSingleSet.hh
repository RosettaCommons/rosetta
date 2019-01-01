// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCSingleSet.hh
/// @brief   class that stores and handles data for one single RDC dataset (i.e. of one type of dipolar coupling in one alignment medium)
/// @details last Modified: 07/27/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_rdc_RDCSingleSet_HH
#define INCLUDED_core_scoring_nmr_rdc_RDCSingleSet_HH

// Unit headers
#include <core/scoring/nmr/rdc/RDCSingleSet.fwd.hh>

// Package headers
#include <core/scoring/nmr/rdc/RDCSingle.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

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
namespace rdc {

class RDCSingleSet {

	friend class RDCMultiSet; // We have to access the RDCSingle Vector through the class RDCMultiSet
	// but I don't want to make that data accessible for any other class
public: // Methods

	/// @brief construct from experiment file
	///        set default values for RDCSingleSet weight and single_rdc_weighting_scheme
	RDCSingleSet(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief constructor with full argument list
	RDCSingleSet(
		std::string const & filename,
		pose::Pose const & pose,
		Real const weight,
		std::string single_rdc_weigting = "CONST"
	);

	/// @brief copy constructor
	RDCSingleSet(RDCSingleSet const & other);

	/// @brief assignment operator
	RDCSingleSet &
	operator=(RDCSingleSet const & rhs);

	/// @brief destructor
	~RDCSingleSet();

	// Getters
	Real get_weight() const { return weight_; }
	Size get_number_rdc() const { return number_rdc_; }
	std::string get_dataset_name() const { return dataset_name_; }
	utility::vector1<RDCSingle> const & get_single_rdc_vec() const { return rdc_single_vec_; }
	SINGLE_NMR_VALUE_WEIGHTING get_single_rdc_weighting_scheme() const { return single_rdc_weighting_scheme_; }
	RDC_TYPE get_rdc_type() const { return rdc_type_; }

	// Setters
	void set_weight(Real weight) { weight_ = weight; }
	void set_single_rdc_weighting_scheme(std::string const & weighting_scheme);

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	RDCSingleSet();

	/// @brief utility function used in constructor to initialize RDCSingelSet object from data file.
	void
	init_from_rdc_filedata(
		std::string const & filename,
		pose::Pose const & pose
	);

private: // Data

	std::string dataset_name_;
	utility::vector1<RDCSingle> rdc_single_vec_;
	Real weight_;    // weight the experiment score by this factor
	Size number_rdc_;
	SINGLE_NMR_VALUE_WEIGHTING single_rdc_weighting_scheme_;
	RDC_TYPE rdc_type_;
};

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_rdc_RDCSingleSet_HH
