// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_filters/AbinitioBaseFilter.hh
/// @brief header file for AbinitioBaseFilter.cc
/// @details
///
///
///
/// @author Robert Vernon and Oliver Lange

#ifndef INCLUDED_protocols_simple_filters_AbinitioBaseFilter_hh
#define INCLUDED_protocols_simple_filters_AbinitioBaseFilter_hh

// Unit Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers
//#include <utility/pointer/ReferenceCount.hh>

//// C++ headers

namespace protocols {
namespace simple_filters {

class AbinitioBaseFilter : public protocols::filters::Filter {
public:

	/// c-tor and d-tor
	AbinitioBaseFilter();
	~AbinitioBaseFilter() override = default;

	filters::FilterOP clone() const override = 0;

	filters::FilterOP fresh_instance() const override = 0;


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	bool apply( core::pose::Pose const & pose ) const override = 0;

	std::string get_protein_sstype( core::pose::Pose const & pose ) const;


	std::string name() const override { return "AbinitioBaseFilter"; };

protected:
	mutable int max_helix_length_;
	mutable core::Real max_helix_fraction_;
	mutable int beta_;
	mutable float beta_ratio_;
	mutable std::string sstype_;
};

} // filters
} // protocols

#endif
