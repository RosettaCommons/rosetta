// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SymUnsatHbondFilter.hh
/// @brief definition of filter class SymUnsatHbondFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_devel_matdes_SymUnsatHbondFilter_hh
#define INCLUDED_devel_matdes_SymUnsatHbondFilter_hh

#include <devel/matdes/SymUnsatHbondFilter.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace devel {
namespace matdes {

/// @brief filters based on an upper bound # of buried unsatisfied polar residues
class SymUnsatHbondFilter : public protocols::filters::Filter
{
public:
	SymUnsatHbondFilter() : protocols::filters::Filter( "SymUnsatHbonds" ) {}
	SymUnsatHbondFilter( core::Size const upper_threshold, core::Size const jump_num );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const {
		return new SymUnsatHbondFilter( *this );
	}
	protocols::filters::FilterOP fresh_instance() const{
		return new SymUnsatHbondFilter();
	}

	virtual ~SymUnsatHbondFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
	core::Real compute( core::pose::Pose const & pose ) const;
	core::Size upper_threshold_;
	core::Size jump_num_;
};

}
}
#endif
