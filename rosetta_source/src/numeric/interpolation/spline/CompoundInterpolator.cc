// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/CompoundInterpolator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///


#include <numeric/interpolation/spline/CompoundInterpolator.hh>

#include <algorithm>

namespace numeric {
namespace interpolation {
namespace spline {

struct compare_interp_range {
  bool operator()( interp_range const & a, interp_range const & b ) {
		return a.lb < b.ub;
  }
};

void
CompoundInterpolator::add_range(
	InterpolatorOP interp,
	Real lb,
	Real ub
) {
	interp_range ir;
	ir.lb = lb;
	ir.ub = ub;
	ir.interp = interp;
	interpolators_.push_back( ir );
	std::sort( interpolators_.begin(), interpolators_.end(), compare_interp_range() );
	for( size_t i = 1; i < interpolators_.size(); ++i ) {
		assert( interpolators_[i].ub <= interpolators_[i+1].lb );
	}
}


void
CompoundInterpolator::interpolate(
	Real x,
	Real & y,
	Real & dy
) {

	if(has_lb_function() && x < get_lb_function_cutoff())
	{
		return compute_lb_function_solution(x,y);
	}
	if(has_ub_function() && x > get_ub_function_cutoff())
	{
		return compute_ub_function_solution(x,y);
	}

	for( size_t i = 1; i <= interpolators_.size(); ++i ) {
		if( interpolators_[i].lb <= x && x <= interpolators_[i].ub ) {
			return interpolators_[i].interp->interpolate(x,y,dy);
		}
	}
	assert(false);
}


} // end namespace spline
} // end namespace interpolation
} // end namespace numeric
