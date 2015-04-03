// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueIEFilter.hh
/// @brief definition of filter class ResidueIEFilter.
/// @author Sagar Khare (khares@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_ResidueIEFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueIEFilter_hh

#include <protocols/simple_filters/ResidueIEFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class ResidueIEFilter : public filters::Filter
{
public:
	ResidueIEFilter() : filters::Filter( "ResidueIE" ) {}

	ResidueIEFilter( utility::vector1<core::Size> const resnums, std::string const restype, core::scoring::ScoreFunctionCOP scorefxn,
	   core::scoring::ScoreType const score_type, core::Real const threshold,
	   bool const whole_pose = false, bool const whole_interface = false, core::Size const rb_jump = 1,
	   core::Real const interface_distance_cutoff =  8.0 , core::Real max_penalty = 0.0, core::Real penalty_factor = 1.0 );

	ResidueIEFilter( ResidueIEFilter const &init );
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new ResidueIEFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ResidueIEFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~ResidueIEFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	utility::vector1<core::Size> resnums() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::scoring::ScoreType score_type() const;
	core::Real threshold() const;
	void resnums( utility::vector1<core::Size> const & rn );
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void score_type( core::scoring::ScoreType score_type );
	void threshold( core::Real const th );
private:
	mutable utility::vector1<core::Size> resnums_;
	std::string restype_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_;
	core::Real threshold_;
	bool whole_pose_;
	bool whole_interface_;
	core::Size rb_jump_;
	core::Real interface_distance_cutoff_;
	core::Real max_penalty_;
	core::Real penalty_factor_;
	bool use_resE_;
	core::pack::task::residue_selector::ResidueSelectorCOP selector_;

};

}
}

#endif
