// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/DdgFilter.hh
/// @brief definition of filter class DdgFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_DdgFilter_hh
#define INCLUDED_protocols_simple_filters_DdgFilter_hh

#include <protocols/simple_filters/DdgFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>



namespace protocols {
namespace simple_filters {

class DdgFilter : public filters::Filter
{
public:
	DdgFilter();
	DdgFilter( core::Real const ddg_threshold, core::scoring::ScoreFunctionCOP scorefxn, core::Size const rb_jump=1, core::Size const repeats=1, bool const symmetry=false );
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return new DdgFilter( *this );
	}
	filters::FilterOP fresh_instance() const{
		return new DdgFilter();
	}

	void repack( bool const repack );
	bool repack() const;
	void repeats( core::Size const repeats );
	core::Size repeats() const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~DdgFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real ddg_threshold_; //dflt -15
	core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL/score12 in cstrctr/rosettascripts
	core::Size rb_jump_; // dflt 1
	core::Size repeats_;//average of how many repeats? defaults to 1
	bool symmetry_; //dflt false
	bool repack_; //dflt true; Do you want to repack in the bound and unbound states (ddG) or merely compute the dG
};

	
}
}
#endif
