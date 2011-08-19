// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/filters/COFilters.hh
/// @brief header file for SheetCOFilter class.
/// @detailed
/// @author Robert Vernon
/// @author James Thompson

#ifndef INCLUDED_protocols_filters_COFilter_hh
#define INCLUDED_protocols_filters_COFilter_hh

// Unit Headers
#include <protocols/filters/AbinitioBaseFilter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <iostream>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace filters {

class COFilter : public AbinitioBaseFilter {
public:
	/// c-tor and d-tor
	COFilter() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
	}
	virtual ~COFilter() {}

	FilterOP clone() const {
		return new COFilter( *this ); }

	FilterOP fresh_instance() const{
		return new COFilter(); }


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	virtual
	bool apply( core::pose::Pose const & pose ) const;

	virtual
	std::string name() const {
		return "Contact-Order Filter";
	}
};

} // filters
} // protocols

#endif
