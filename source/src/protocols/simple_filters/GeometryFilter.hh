// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/GeometryFilter.hh
/// @brief definition of filter class GeometryFilter.
/// @author Possu Huang (possu@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_GeometryFilter_hh
#define INCLUDED_protocols_simple_filters_GeometryFilter_hh

#include <protocols/simple_filters/GeometryFilter.fwd.hh>


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class GeometryFilter : public filters::Filter
{
public:
	GeometryFilter();
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new GeometryFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new GeometryFilter() );
	}

	~GeometryFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	core::Real omega_cutoff_;
	core::Real cart_bonded_cutoff_;
	std::string filename_;
	core::Real cst_cutoff_;
	core::Size start_;
	core::Size end_;
	bool count_bad_residues_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
};

}
}
#endif
