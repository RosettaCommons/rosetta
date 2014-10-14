// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/filters/SecondaryStructureHasResidueFilter.hh
/// @brief header file for SecondaryStructureHasResidueFilter class.
/// @detailed
/// @author Chris King ( chrisk1@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_filters_SecondaryStructureHasResidueFilter_hh
#define INCLUDED_protocols_fldsgn_filters_SecondaryStructureHasResidueFilter_hh

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureHasResidueFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector0.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace fldsgn {
namespace filters {

class SecondaryStructureHasResidueFilter : public protocols::filters::Filter {
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
	SecondaryStructureHasResidueFilter();

	virtual ~SecondaryStructureHasResidueFilter(){}


public:// virtual constructor


	// @brief make clone
	virtual FilterOP clone() const { return FilterOP( new SecondaryStructureHasResidueFilter( *this ) ); }

	// @brief make fresh instance
	virtual FilterOP fresh_instance() const {	return FilterOP( new SecondaryStructureHasResidueFilter() ); }

public:// accessor


	// @brief get name of this filter
	virtual std::string name() const { return "SecondaryStructureHasResidueFilter"; }


public:// parser

	virtual void parse_my_tag( TagCOP tag,
														 basic::datacache::DataMap &,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );


public:// main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	bool apply( Pose const & pose ) const;
	core::Real report_sm( Pose const & pose ) const;
	void report( std::ostream & out,  Pose const & pose ) const;
	core::Real compute( Pose const & pose ) const;
	core::Size n_req_res_in_seq( std::string const &, utility::vector0< bool > const & is_checked ) const;
private:

	core::Size min_helix_length_;
	core::Size max_helix_length_;
	core::Size min_sheet_length_;
	core::Size max_sheet_length_;
	core::Size min_loop_length_;
	core::Size max_loop_length_;

	core::Real threshold_;
	std::string req_residue_str_;
	core::Size nres_req_per_secstruct_;
	core::pack::task::TaskFactoryOP res_check_task_factory_;
	core::pack::task::TaskFactoryOP ss_select_task_factory_;

	bool filter_helix_;
	bool filter_sheet_;
	bool filter_loop_;

};

} // filters
} // fldsgn
} // protocols

#endif
	// In this case, the test is whether the give pose is the topology we want.
