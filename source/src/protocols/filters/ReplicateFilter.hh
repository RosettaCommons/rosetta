// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/ReplicateFilter.hh
/// @brief Repeat a subfilter multiple times, and pass a value based on the aggregate results
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_filters_ReplicateFilter_hh
#define INCLUDED_protocols_filters_ReplicateFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Unit headers

namespace protocols {
namespace filters {

class ReplicateFilter : public protocols::filters::Filter
{
public:
	ReplicateFilter();
	ReplicateFilter(protocols::filters::FilterOP subfilter, core::Size replicates=1, core::Size upper_trim=0, core::Size lower_trim=0);
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
  core::Real compute( core::pose::Pose const & pose ) const;

	virtual protocols::filters::FilterOP clone() const {
		return protocols::filters::FilterOP( new ReplicateFilter( *this ) );
	}
	virtual protocols::filters::FilterOP fresh_instance() const{
		return protocols::filters::FilterOP( new ReplicateFilter() );
	}

	virtual ~ReplicateFilter(){};
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void subfilter(protocols::filters::FilterOP subfilter) { subfilter_ = subfilter; }
	void median(bool median) { median_ = median; }
	void threshold( core::Real threshold) { threshold_ = threshold; }

private:
	protocols::filters::FilterOP subfilter_;
	core::Size replicates_;
	core::Size upper_trim_;
	core::Size lower_trim_;
  bool median_;
	core::Real threshold_;

};

} // filters
} // protocols

#endif //INCLUDED_protocols_Filters_ReplicateFilter_HH_
