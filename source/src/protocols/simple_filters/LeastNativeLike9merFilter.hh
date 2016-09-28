// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/LeastNativeLike9merFilter.hh
/// @brief  Find the rmsd to the worst 9mer in the chain
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_LeastNativeLike9merFilter_hh
#define INCLUDED_protocols_simple_filters_LeastNativeLike9merFilter_hh

// Unit Headers
#include <protocols/simple_filters/LeastNativeLike9merFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/indexed_structure_store/SSHashedFragmentStore.hh>
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
namespace simple_filters {

class LeastNativeLike9merFilter : public protocols::filters::Filter{
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
	LeastNativeLike9merFilter();


	~LeastNativeLike9merFilter() override;


public:// virtual constructor

	// @brief make clone
	filters::FilterOP clone() const override { return filters::FilterOP(new LeastNativeLike9merFilter(*this));}
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return filters::FilterOP(new LeastNativeLike9merFilter());}

	// @brief get name of this filter
	std::string name() const override { return "LeastNativeLike9merFilter"; }
	Real report_sm(const Pose & pose ) const override;
	void report( std::ostream & out,const Pose & pose ) const override;
	Real compute( const Pose & pose ) const;
	bool apply(const Pose & pose ) const override;
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		Movers_map const &,
		Pose const & ) override;


private:

	Real filtered_value_;
	Real rmsd_lookup_thresh_;
	core::indexed_structure_store::SSHashedFragmentStore * SSHashedFragmentStore_;
	bool ignore_terminal_res_;
	bool only_helices_;

};

} // filters
} // protocols

#endif

