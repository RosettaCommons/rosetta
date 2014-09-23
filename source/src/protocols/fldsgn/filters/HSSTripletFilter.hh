// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/HSSTripletFilter.hh
/// @brief header file for HSSTripletFilter class.
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_HSSTripletFilter_hh
#define INCLUDED_protocols_fldsgn_filters_HSSTripletFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/HSSTripletFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>

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

class HSSTripletFilter : public protocols::filters::Filter {
public:


	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::HSSTriplet  HSSTriplet;
	typedef protocols::fldsgn::topology::HSSTriplets  HSSTriplets;
	typedef protocols::fldsgn::topology::HSSTripletSet  HSSTripletSet;
	typedef protocols::fldsgn::topology::HSSTripletSetOP  HSSTripletSetOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	HSSTripletFilter();

	// @brief constructor with arguments
	HSSTripletFilter( HSSTriplets const & hss3s );

	// @brief constructor with arguments
	HSSTripletFilter( String const & hss3s );

	// @brief copy constructor
	HSSTripletFilter( HSSTripletFilter const & rval );

	virtual ~HSSTripletFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new HSSTripletFilter( *this ) ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const {	return FilterOP( new HSSTripletFilter() ); }


public:// mutator


	// @brief add hsstriplets for filtering
	void add_hsstriplets( HSSTriplets const & hss3s );

	// @brief set secondary strucure elements
	void secstruct( String const & ss );

	// @brief minimum distance for filtering
	void filter_min_dist( Real const r );

	// @brief maximum distance for filtering
	void filter_max_dist( Real const r );

	/// @brief miniimum angle for filtering
	void filter_min_angle( Real const r );

	/// @brief maximum angle for filtering
	void filter_max_angle( Real const r );

	/// @brief set output id
	void output_id( Size const i );

	/// @brief set output type
	void output_type( String const & s );


public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "HSSTripletFilter"; }


public:// parser


	virtual void parse_my_tag( TagCOP tag,
														 basic::datacache::DataMap &,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );


public:// virtual main operation


	/// @brief
	Real report_sm( Pose const & pose ) const;

	/// @brief
	Real compute( Pose const & pose ) const;


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;


private:


	/// @brief hsstriplet
	HSSTripletSetOP hss3set_;

	/// @brief if value is empty, dssp will run for ss definition ( default is emptry )
	mutable String secstruct_;

	/// @brief filtered min distance between helix and sheet
	Real filter_min_dist_;

	/// @brief filtered max distance between helix and sheet
	Real filter_max_dist_;

	/// @brief filtered min angle between helix and sheet
	Real filter_min_angle_;

	/// @brief filtered max angle between helix and sheet
	Real filter_max_angle_;

	/// @brief output id of HSSTriplet
	Size output_id_;

	/// @brief output type, dist or angle
	String output_type_;

	/// @brief output value, the result of filterring calculation
	mutable Real output_value_;

	/// @brief if set, the filter will ignore the direction of the helix and return an angle between -90 and 90 instead of -180 and 180
	bool ignore_helix_direction_;
};

} // filters
} // fldsgn
} // protocols

#endif
