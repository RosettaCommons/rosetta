// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSData.hh
/// @brief   class that stores and handles all PCS data for all tagging sites and all lanthanides
/// @details last Modified: 06/30/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pcs_PCSData_HH
#define INCLUDED_core_scoring_nmr_pcs_PCSData_HH

// Unit headers
#include <core/scoring/nmr/pcs/PCSData.fwd.hh>

// Package headers
#include <core/scoring/nmr/NMRDataFactory.fwd.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.fwd.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.fwd.hh>
#include <core/scoring/nmr/pcs/PCSTensor.fwd.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <string>
#include <algorithm>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

class PCSData : public basic::datacache::CacheableData {

public: // Methods

	/// @brief construct with filename
	PCSData(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief copy constructor
	PCSData(PCSData const & other);

	/// @brief assignment operator
	PCSData &
	operator=(PCSData const & rhs);

	/// @brief destructor
	~PCSData() override;

	basic::datacache::CacheableDataOP clone() const override;

	/// @brief compute the overall PCS score and individual scores for each tagging site
	Real
	compute_score_all_tags(
		pose::Pose & pose,
		utility::vector1<Real> & scores_all_tags,
		utility::vector1< utility::vector1< PCSTensorCOP > > & tensors_all_lanthanides
	);

	// Getters
	utility::vector1< PCSMultiSetOP > & get_pcs_multiset_vec() { return pcs_multiset_vec_; }
	utility::vector1< PCSMultiSetOP > const & get_pcs_multiset_vec() const { return pcs_multiset_vec_; }
	Size get_number_tags() const { return number_tags_; }
	bool optimize_tensors() const { return optimize_tensors_; }
	Size get_total_number_pcs() const;

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	PCSData();

	/// @brief register options
	void register_options();
	void init_from_cml();

	/// @brief utility function used during construction of PCSData object
	void
	init_pcs_data_from_file(
		std::string const & filename,
		pose::Pose const & /*pose*/
	);

private: // Data

	utility::vector1< PCSMultiSetOP > pcs_multiset_vec_;
	Size number_tags_;
	bool optimize_tensors_;

};

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pcs_PCSData_HH
