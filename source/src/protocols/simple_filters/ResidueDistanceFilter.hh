// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueDistanceFilter.hh
/// @brief definition of filter class ResidueDistanceFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_ResidueDistanceFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueDistanceFilter_hh

#include <protocols/simple_filters/ResidueDistanceFilter.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.hh>

namespace protocols {
namespace simple_filters {

class ResidueDistanceFilter : public filters::Filter
{
public:
	ResidueDistanceFilter() : filters::Filter( "ResidueDistance"  ) {}
	ResidueDistanceFilter( core::Size const res1, core::Size const res2, core::Real const distance_threshold ) :
		Filter( "ResidueDistance" ), res1_( res1 ), res2_( res2 ), distance_threshold_( distance_threshold ) {}
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new ResidueDistanceFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new ResidueDistanceFilter() );
	}

	~ResidueDistanceFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	core::Size res1_, res2_;
	core::Real distance_threshold_;
};

}
}

#endif
