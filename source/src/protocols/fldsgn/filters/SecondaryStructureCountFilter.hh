// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/SecondaryStructureCountFilter.hh
/// @brief header file for SecondaryStructureCountFilter class.
/// @detailed
/// @author Lei Shi ( shilei@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_SecondaryStructureCountFilter_hh
#define INCLUDED_protocols_fldsgn_filters_SecondaryStructureCountFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureCountFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

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

class SecondaryStructureCountFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	SecondaryStructureCountFilter();

	virtual ~SecondaryStructureCountFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return new SecondaryStructureCountFilter( *this ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const {	return new SecondaryStructureCountFilter(); }

public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "SecondaryStructureCountFilter"; }


public:// parser

	virtual void parse_my_tag( TagCOP tag,
														 basic::datacache::DataMap &,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );


public:// main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const;
	core::Real report_sm( Pose const & pose ) const;
	void report( std::ostream & out,  Pose const & pose ) const;
	core::Size compute( Pose const & pose ) const;
private:

	mutable core::Size num_helix_pose_;
	mutable core::Size num_sheet_pose_;
	mutable core::Size num_loop_pose_;
	core::Size num_helix_;
	core::Size num_sheet_;
	core::Size num_loop_;
	core::Size num_helix_sheet_;

	core::Size min_helix_length_;
	core::Size max_helix_length_;
	core::Size min_sheet_length_;
	core::Size max_sheet_length_;
	core::Size min_loop_length_;
	core::Size max_loop_length_;

	bool filter_helix_;
	bool filter_sheet_;
	bool filter_loop_;
	bool filter_helix_sheet_;

};

} // filters
} // fldsgn
} // protocols

#endif
