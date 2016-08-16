// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueSetChainEnergyFilter.hh
/// @brief definition of filter class ResidueSetChainEnergyFilter.
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_ResidueSetChainEnergyFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueSetChainEnergyFilter_hh

//unit headers
#include <protocols/simple_filters/ResidueSetChainEnergyFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>


namespace protocols {
namespace simple_filters {

class ResidueSetChainEnergyFilter : public filters::Filter
{
public:
	//default ctor
	ResidueSetChainEnergyFilter();
	//full ctor
	ResidueSetChainEnergyFilter( core::scoring::ScoreFunctionCOP scorefxn,
		core::scoring::ScoreType const score_type,
		core::Real const score_type_threshold,
		std::string string_resnums,
		core::Size chain
	);
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new ResidueSetChainEnergyFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ResidueSetChainEnergyFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	virtual ~ResidueSetChainEnergyFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real score_type_threshold_;
	core::scoring::ScoreType score_type_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::string string_resnums_;
	core::Size chain_;
};

}
}

#endif
