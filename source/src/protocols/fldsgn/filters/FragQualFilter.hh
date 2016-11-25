// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/filters/FragQualFilter.hh
/// @brief header file for FragQualFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_FragQualFilter_hh
#define INCLUDED_protocols_fldsgn_filters_FragQualFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/FragQualFilter.fwd.hh>

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

class FragQualFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::string String;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	FragQualFilter();

	// @brief constructor with arguments
	//FragQualFilter( Real const & ss );

	// @brief copy constructor
	FragQualFilter( FragQualFilter const & rval );

	virtual ~FragQualFilter(){}


public:// virtual constructor


	// @brief make clone
	FilterOP clone() const override { return FilterOP( new FragQualFilter( *this ) ); }

	// @brief make fresh instance
	FilterOP fresh_instance() const override { return FilterOP( new FragQualFilter() ); }


public:// mutator


	// @brief
	void filtered_value( Real const & ss );

	// @brief
	void filtered_type( String const & ss );


public:// accessor


	// @brief get name of this filter
	// XRW TEMP  virtual std::string name() const { return "FragQualFilter"; }


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose ) override;


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	/// @brief
	Real report_sm( Pose const & pose ) const override;

	/// @brief used to report score
	void report( std::ostream & out, Pose const & pose ) const override;

	/// @brief
	Real compute( Pose const & pose ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:


	String filtered_type_;

	Real filtered_value_;

	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Real rmsd_cutoff_;


};

} // filters
} // fldsgn
} // protocols

#endif
