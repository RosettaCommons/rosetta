// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/filters/HSSTripletFilter.hh
/// @brief header file for HSSTripletFilter class.
/// @details
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
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::HSSTriplet  HSSTriplet;
	typedef protocols::fldsgn::topology::HSSTriplets  HSSTriplets;
	typedef protocols::fldsgn::topology::HSSTripletSet  HSSTripletSet;
	typedef protocols::fldsgn::topology::HSSTripletSetOP  HSSTripletSetOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	/// @brief default constructor
	HSSTripletFilter();

	/// @brief constructor with arguments
	HSSTripletFilter( HSSTriplets const & hss3s );

	/// @brief constructor with arguments
	HSSTripletFilter( String const & hss3s );

	/// @brief copy constructor -- required because we clone the hss3set_ pointer
	HSSTripletFilter( HSSTripletFilter const & rval );

	virtual ~HSSTripletFilter(){}


public:// virtual constructor


	// @brief make clone
	FilterOP clone() const override { return FilterOP( new HSSTripletFilter( *this ) ); }

	// @brief make fresh instance
	FilterOP fresh_instance() const override { return FilterOP( new HSSTripletFilter() ); }


public:// mutator

	/// @brief add filtered HSSTriplets from string
	void add_hsstriplets( String const & hss3s );

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
	// XRW TEMP  virtual std::string name() const { return "HSSTripletFilter"; }


public:// parser


	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & ) override;


public:// virtual main operation


	/// @brief
	Real report_sm( Pose const & pose ) const override;

	/// @brief
	Real compute( Pose const & pose ) const;


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	/// @brief    computes and returns secondary structure string to use for this filter.
	/// @details  if secstruct_ is non-empty, returns that
	///           if use_dssp_ is true, use DSSP to compute secstruct
	///           otherwise, use secstruct stored in the pose
	String
	get_secstruct( Pose const & pose ) const;

	/// @brief    returns HSSTriplets object to use for this filter
	/// @details  if hss triplets are given prior to apply time, returns those
	///           otherwise, look for HSS info in the pose's StructureData and return that
	HSSTriplets
	get_hss3s( Pose const & pose ) const;

	/// @brief    checks secondary structure elements in the triplet, returns false if invalid
	/// @details  If the pose doesn't contain the helix or strands, returns false
	///           If the length of the helix is < 5, returns false
	///           If the length of either strand is < 2, returns false
	bool
	check_elements( HSSTriplet const & hss, SS_Info2 const & ss_info ) const;

	/// @brief given an HSS triplet, compute the distance from helix to sheet
	Real
	compute_dist( HSSTriplet const & hss ) const;

	/// @brief given an hs-angle, return a valid angle accounting for ignore_helix_direction_
	///        if ignore_helix_direction_ is true, this basically makes angle periodic from
	///        -90 to 90
	Real
	compute_angle( Real const angle ) const;

private:

	/// @brief hsstriplet
	HSSTripletSetOP hss3set_;

	/// @brief if value is empty, dssp will run for ss definition ( default is emptry )
	String secstruct_;

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

	/// @brief if set, the filter will ignore the direction of the helix and return an angle between -90 and 90 instead of -180 and 180
	bool ignore_helix_direction_;

	/// @brief if set, and secstruct_ is empty, secondary structure will be determined by dssp.  Otherwise, it will be taken from the pose.
	bool use_dssp_;
};

} // filters
} // fldsgn
} // protocols

#endif
