// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCData.hh
/// @brief   class that stores and handles all RDC data for all alignment media and experiments
/// @details last Modified: 07/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_rdc_RDCData_HH
#define INCLUDED_core_scoring_nmr_rdc_RDCData_HH

// Unit headers
#include <core/scoring/nmr/rdc/RDCData.fwd.hh>

// Package headers
#include <core/scoring/nmr/NMRDataFactory.fwd.hh>
#include <core/scoring/nmr/rdc/RDCMultiSet.fwd.hh>
#include <core/scoring/nmr/rdc/RDCSingleSet.fwd.hh>
#include <core/scoring/nmr/rdc/RDCTensor.fwd.hh>

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
namespace rdc {

class RDCData : public basic::datacache::CacheableData {

public: // Methods

	/// @brief construct from filename
	RDCData(
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief copy constructor
	RDCData(RDCData const & other);

	/// @brief assignment operator
	RDCData & operator=(RDCData const & rhs);

	/// @brief destructor
	~RDCData() override;

	basic::datacache::CacheableDataOP clone() const override;

	/// @brief compute the overall RDC score and scores for the individual alignment media
	Real
	compute_score_all_media(
		pose::Pose const & pose,
		utility::vector1<Real> & scores_all_media,
		utility::vector1< RDCTensorCOP > & tensors_all_media
	);

	// Getters
	utility::vector1< RDCMultiSetOP > & get_rdc_multiset_vec() { return rdc_multiset_vec_; }
	utility::vector1< RDCMultiSetOP > const & get_rdc_multiset_vec() const { return rdc_multiset_vec_; }
	Size get_number_alignment_media() const { return number_alignment_media_; }
	Size get_total_number_rdc() const;

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	RDCData();

	/// @brief register options
	void register_options();

	/// @brief utility function used during construction of RDCData object
	void
	init_rdc_data_from_file(
		std::string const & filename,
		pose::Pose const & pose
	);

private: // Data

	utility::vector1< RDCMultiSetOP > rdc_multiset_vec_;
	Size number_alignment_media_;

};

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_rdc_RDCData_HH
