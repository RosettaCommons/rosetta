// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/filters/COFilters.hh
/// @brief header file for SheetCOFilter class.
/// @details
/// @author Robert Vernon
/// @author James Thompson

#ifndef INCLUDED_protocols_simple_filters_COFilter_hh
#define INCLUDED_protocols_simple_filters_COFilter_hh

// Unit Headers
#include <protocols/simple_filters/AbinitioBaseFilter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace simple_filters {

class COFilter : public AbinitioBaseFilter {
public:
	/// c-tor and d-tor
	COFilter() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
	}
	~COFilter() override = default;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new COFilter( *this ) ); }

	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new COFilter() ); }


	/// @brief Returns true if the given pose passes the filter, false otherwise.

	bool apply( core::pose::Pose const & pose ) const override;


	std::string name() const override {
		return "Contact-Order Filter";
	}
};

} // filters
} // protocols

#endif
