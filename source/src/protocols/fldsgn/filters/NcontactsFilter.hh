// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/NcontactsFilter.hh
/// @brief header file for NcontactsFilter class.
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_NcontactsFilter_hh
#define INCLUDED_protocols_fldsgn_filters_NcontactsFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/NcontactsFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

//// C++ headers
#include <map>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace filters {

class NcontactsFilter : public protocols::filters::Filter {
public:


	typedef protocols::filters::Filter Super;
	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	/// @brief default constructor
	NcontactsFilter();

	/// @brief default constructor
	NcontactsFilter( String const & report_type, Real const filter_value );

	/// @brief copy constructor
	NcontactsFilter( NcontactsFilter const & rval );

	/// @brief destructor
	virtual ~NcontactsFilter();


public:// virtual constructor


	/// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new NcontactsFilter( *this ) ); }

	/// @brief make fresh instance
	virtual FilterOP fresh_instance() const {	return FilterOP( new NcontactsFilter() ); }


public:// accessor


	/// @brief get name of this filter
	virtual std::string name() const { return "Ncontacts"; }


public:// main calculator


	/// @brief compute number of contacts
	Real compute( Pose const & pose ) const;


public:// virtual main operations


	/// @brief used to report score
	virtual Real report_sm( Pose const & pose ) const;

	/// @brief used to report score
	virtual void report( std::ostream & out, Pose const & pose ) const;

	/// @brief returns true if the given pose passes the filter, false otherwise.
	/// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;


public:// parser


	virtual void parse_my_tag( TagCOP tag,
														 basic::datacache::DataMap &,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );


private:

	/// @brief
	String report_type_;

	/// @brief
	Real filter_value_;

}; //NontactFilter


} // filters
} // fldsgn
} // protocols

#endif
