// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/filters/HelixKinkFilter.hh
/// @brief header file for HelixKinkFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_HelixKinkFilter_hh
#define INCLUDED_protocols_fldsgn_filters_HelixKinkFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/HelixKinkFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace fldsgn {
namespace filters {

class HelixKinkFilter : public protocols::filters::Filter {
public:


	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;

	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	HelixKinkFilter();

	virtual ~HelixKinkFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new HelixKinkFilter( *this ) ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const { return FilterOP( new HelixKinkFilter() ); }


public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "HelixKinkFilter"; }


public:// parser


	virtual void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;


private:


	/// @brief
	Real bend_angle_;

	/// @brief
	String secstruct_;

	/// @brief
	String string_resnums_;
	bool select_resnums_;
	bool select_range_;
	core::Size helix_start_;
	core::Size helix_end_;

};

} // filters
} // fldsgn
} // protocols

#endif
