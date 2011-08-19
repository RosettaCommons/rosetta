// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/filters/ScFilter.hh
/// @brief  header file for ScFilter class
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)


#ifndef INCLUDED_protocols_filters_ScFilter_hh
#define INCLUDED_protocols_filters_ScFilter_hh

// Unit Headers
#include <protocols/filters/ScFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Parser headers
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

//// C++ headers

namespace protocols {
namespace filters {

class ScFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagPtr TagPtr;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	ScFilter();

	// @brief constructor with arguments
	ScFilter( Real const & filtered_sc, Size const & jump_id, Size const & quick, Size const & verbose );

	// @brief copy constructor
	ScFilter( ScFilter const & rval );

	virtual ~ScFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return new ScFilter( *this ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const {	return new ScFilter(); }


public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "ScFilter"; }


public:// mutator

	void filtered_sc( Real const & filtered_sc );
	void jump_id( Size const & jump_id );
	void quick( Size const & quick );
	void verbose( Size const & verbose );

public:// parser

	virtual void parse_my_tag( TagPtr const tag,
		DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;

	/// @brief
	virtual Real report_sm( Pose const & pose ) const;

	/// @brief calc packstat score
	Real compute( Pose const & pose ) const;


private:

	Real filtered_sc_;
	Size verbose_;
	Size quick_;
	Size jump_id_;

};

} // filters
} // protocols

#endif
