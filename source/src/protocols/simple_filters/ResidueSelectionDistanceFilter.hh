// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueSelectionDistanceFilter.hh
/// @brief loops over the residues in a residue selector and reports the average distance
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_simple_filters_ResidueSelectionDistanceFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueSelectionDistanceFilter_hh

#include <protocols/simple_filters/ResidueSelectionDistanceFilter.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace protocols {
namespace simple_filters {

class ResidueSelectionDistanceFilter : public filters::Filter
{
public:
	ResidueSelectionDistanceFilter() : filters::Filter( "ResidueDistance"  ) {}
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return filters::FilterOP( new ResidueSelectionDistanceFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ResidueSelectionDistanceFilter() );
	}

	virtual ~ResidueSelectionDistanceFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real distance_threshold_;
	std::string atom_to_measure_;
	core::select::residue_selector::ResidueSelectorCOP selector_;

};

}
}

#endif
