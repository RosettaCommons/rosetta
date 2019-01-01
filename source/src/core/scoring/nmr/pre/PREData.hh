// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREData.hh
/// @brief   class that stores all PRE data for all spinlabel sites
/// @details last Modified: 10/12/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pre_PREData_HH
#define INCLUDED_core_scoring_nmr_pre_PREData_HH

// Unit headers
#include <core/scoring/nmr/pre/PREData.fwd.hh>

// Package headers
#include <core/scoring/nmr/NMRDataFactory.fwd.hh>
#include <core/scoring/nmr/pre/PREMultiSet.fwd.hh>
#include <core/scoring/nmr/pre/PRESingleSet.fwd.hh>

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
namespace pre {

class PREData : public basic::datacache::CacheableData {

public: // Methods

	/// @brief constructor with filename
	PREData(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief copy constructor
	PREData(PREData const & other);

	/// @brief assignment operator
	PREData &
	operator=(PREData const & rhs);

	/// @brief destructor
	~PREData() override;

	basic::datacache::CacheableDataOP clone() const override;

	/// @brief compute the overall PRE score and individual scores for each spinlabel site
	Real
	compute_score_all_spinlabel(
		pose::Pose & pose,
		utility::vector1<Real> & individual_scores
	);

	// Getters
	utility::vector1< PREMultiSetOP > & get_pre_multiset_vec() { return pre_multiset_vec_; }
	utility::vector1< PREMultiSetOP > const & get_pre_multiset_vec() const { return pre_multiset_vec_; }
	Size get_number_spinlabel_sites() const { return number_spinlabel_sites_; }
	Size get_total_number_pre() const;

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	PREData();

	/// @brief register options
	void register_options();

	/// @brief utility function used during construction of PREData object
	void
	init_pre_data_from_file(
		std::string const & filename,
		pose::Pose const & pose
	);

private: // Data

	utility::vector1< PREMultiSetOP > pre_multiset_vec_;
	Size number_spinlabel_sites_;

};

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pre_PREData_HH
