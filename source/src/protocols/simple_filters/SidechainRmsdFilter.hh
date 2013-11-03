// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SidechainRmsdFilter.hh
/// @brief A filter based on automorphic sidechain RMSD
/// @author Noah Ollikainen

#ifndef INCLUDED_protocols_simple_filters_SidechainRmsdFilter_hh
#define INCLUDED_protocols_simple_filters_SidechainRmsdFilter_hh

#include <protocols/simple_filters/SidechainRmsdFilter.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.hh>

namespace protocols {
namespace simple_filters {

class SidechainRmsdFilter : public filters::Filter
{
public:
	SidechainRmsdFilter();
	SidechainRmsdFilter( core::Size const res1, core::Size const res2, core::Real const rmsd_threshold );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const;
	filters::FilterOP fresh_instance() const;

	virtual ~SidechainRmsdFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &);
private:
	core::Size res1_, res2_;
	core::Real rmsd_threshold_;
	core::pose::PoseOP reference_pose_;
	bool include_backbone_;
};

}
}

#endif
