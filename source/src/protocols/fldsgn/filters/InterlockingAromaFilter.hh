// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/filters/InterlockingAromaFilter.hh
/// @brief header file for InterlockingAromaFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_InterlockingAromaFilter_hh
#define INCLUDED_protocols_fldsgn_filters_InterlockingAromaFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/InterlockingAromaFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

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

class InterlockingAromaFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::string String;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	InterlockingAromaFilter();

	// @brief constructor with arguments
	InterlockingAromaFilter( String const & type );

	// @brief copy constructor
	InterlockingAromaFilter( InterlockingAromaFilter const & rval );

	virtual ~InterlockingAromaFilter(){}


public:// virtual constructor


	// @brief make clone
	FilterOP clone() const override { return FilterOP( new InterlockingAromaFilter( *this ) ); }

	// @brief make fresh instance
	FilterOP fresh_instance() const override { return FilterOP( new InterlockingAromaFilter() ); }


public:// mutator


	// @brief
	void filter_value( Real const value );

	// @brief
	void contact_distance( Real const value );

	// @brief
	void verbose( bool const b );


public:// accessor


	// @brief get name of this filter
	// XRW TEMP  virtual std::string name() const { return "InterlockingAromaFilter"; }


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & pose ) override;


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	/// @brief
	Real report_sm( Pose const & pose ) const override;

	/// @brief used to report score
	void report( std::ostream & out, Pose const & pose ) const override;

	/// @brief compute this filter
	Real compute( Pose const & pose ) const;

	/// @brief compute this filter for a give residue
	bool compute( Size const & res, Pose const & pose, SS_Info2_COP const ssinfo ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:


	Real filter_value_;
	Real contact_dist2_;
	String input_ss_;
	bool verbose_;


};

} // filters
} // fldsgn
} // protocols

#endif
