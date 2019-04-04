// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/SSElementLengthFilter.hh
/// @brief filter structures by longest,shortest or avg length of a given secondary structure type

/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_SSElementLengthFilter_hh
#define INCLUDED_protocols_simple_filters_SSElementLengthFilter_hh

// Unit Headers
#include <protocols/simple_filters/SSElementLengthFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
// Utility headers

// Parser headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>
#include <set>

//// C++ headers

namespace protocols {
namespace simple_filters {

class SSElementLengthFilter : public protocols::filters::Filter{
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
	SSElementLengthFilter();

	// @brief copy constructor
	SSElementLengthFilter( SSElementLengthFilter const & rval );

	virtual ~SSElementLengthFilter();


public:// virtual constructor


	// @brief make clone
	filters::FilterOP clone() const override { return filters::FilterOP(utility::pointer::make_shared<SSElementLengthFilter>(*this));}
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return filters::FilterOP(utility::pointer::make_shared<SSElementLengthFilter>());}


public:// accessor


	// @brief get name of this filter


public:// virtual main operation


	Real report_sm(const Pose & pose ) const override;
	void report( std::ostream & out,const Pose & pose ) const override;
	protocols::loops::Loops get_ss_elements(const Pose & pose) const;
	Real compute( const Pose & pose ) const;
	bool apply(const Pose & pose ) const override;


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		Movers_map const & ,
		Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	Real threshold_;
	bool report_avg_;
	bool report_longest_;
	bool report_shortest_;
	core::select::residue_selector::ResidueSelectorCOP selector_;

};

} // filters
} // protocols

#endif
