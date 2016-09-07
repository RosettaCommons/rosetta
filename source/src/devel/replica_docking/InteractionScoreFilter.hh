// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_InteractionScoreFilter_hh
#define INCLUDED_devel_replica_docking_InteractionScoreFilter_hh

#include <devel/replica_docking/InteractionScoreFilter.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

namespace devel {
namespace replica_docking {

class InteractionScoreFilter : public protocols::filters::Filter
{
public:
	InteractionScoreFilter();
	InteractionScoreFilter( std::string const scorefxn_name, core::Size const rb_jump=1, core::Real const lower_threshold=-30, core::Real const upper_threshold=0.0 );
	InteractionScoreFilter( core::scoring::ScoreFunctionCOP scorefxn, core::Size const rb_jump=1, core::Real const lower_threshold=-30, core::Real const upper_threshold=0.0 );
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const; // which residue numbers are neighbors
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override;

	~InteractionScoreFilter() override;
	void jump( core::Size const jump );
	core::Size jump() const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	core::Real lower_threshold_;
	core::Real upper_threshold_;
	core::Size jump_; // dflt 1; across which jump to compute sasa
	core::scoring::ScoreFunctionOP scorefxn_;
};

}
}

#endif
