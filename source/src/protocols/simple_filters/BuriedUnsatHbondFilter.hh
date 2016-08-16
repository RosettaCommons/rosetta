// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BuriedUnsatHbondFilter.hh
/// @brief definition of filter class BuriedUnsatHbondFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_BuriedUnsatHbondFilter_hh
#define INCLUDED_protocols_simple_filters_BuriedUnsatHbondFilter_hh

#include <protocols/simple_filters/BuriedUnsatHbondFilter.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

namespace protocols {
namespace simple_filters {

/// @brief filters based on an upper bound # of buried unsatisfied polar residues
class BuriedUnsatHbondFilter : public filters::Filter
{
public:
	BuriedUnsatHbondFilter();
	BuriedUnsatHbondFilter( core::Size const upper_threshold, core::Size const jump_num );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const;
	filters::FilterOP fresh_instance() const;

	virtual ~BuriedUnsatHbondFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::pack::task::TaskFactoryOP task_factory() const;
private:
	core::scoring::ScoreFunctionCOP sfxn_;
	core::Size upper_threshold_;
	core::Size jump_num_;
	core::pack::task::TaskFactoryOP task_factory_; // dflt NULL; only residues defined as packable by the taskoperations will be tested for burial
};

}
}
#endif
