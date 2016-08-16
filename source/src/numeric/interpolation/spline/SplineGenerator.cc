// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/spline/SplineGenerator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/interpolation/spline/CompoundInterpolator.hh>

#include <algorithm>
#include <iostream>

namespace numeric {
namespace interpolation {
namespace spline {


struct OrderPoint {
	bool operator() ( Point const & a, Point const & b ) {
		return a.x < b.x;
	}
};

SplineGenerator::SplineGenerator(
	Real lbx, Real lby, Real lbdy,
	Real ubx, Real uby, Real ubdy
) : lbx_ (lbx), lby_ (lby), lbdy_(lbdy),
	ubx_ (ubx), uby_ (uby), ubdy_(ubdy)
{

}

SplineGenerator::SplineGenerator():
	lbx_(0.0),lby_(0.0),lbdy_(0.0),
	ubx_(0.0),uby_(0.0),ubdy_(0.0)
{

}

SplineGenerator::~SplineGenerator() {}

void
SplineGenerator::add_known_value(
	Real x, Real y
) {
	assert( interpolator_ == (InterpolatorOP)NULL );
	assert( lbx_ < x && x < ubx_ );
	points_.push_back( Point(x,y) );
}

void
SplineGenerator::add_known_value(
	Real x, Real y, Real dy
) {
	assert( interpolator_ == (InterpolatorOP)NULL );
	assert( lbx_ < x && x < ubx_ );
	points_.push_back( Point(x,y,dy) );
}

void SplineGenerator::add_boundary_function(std::string const & tag, Real const & cutoff, Real const & slope, Real const & intercept)
{
	boundary_functions_[tag] = LinearFunction(cutoff,slope,intercept);
}

InterpolatorOP
SplineGenerator::get_interpolator()
{
	if ( interpolator_ == (InterpolatorOP)NULL ) {
		std::sort( points_.begin(), points_.end(), OrderPoint() );
		bool compound = false;
		for ( size_t i = 1; i <= points_.size(); ++i ) {
			compound = compound || points_[i].has_dy;
		}
		utility::vector1<Real> x,y;
		x.push_back( lbx_ );
		y.push_back( lby_ );
		if ( !compound ) {
			for ( size_t i = 1; i <= points_.size(); ++i ) {
				x.push_back(points_[i].x);
				y.push_back(points_[i].y);
			}
			x.push_back( ubx_ );
			y.push_back( uby_ );
			interpolator_ = InterpolatorOP( new SimpleInterpolator(x,y,lbdy_,ubdy_) );
		} else {
			// std::cerr << "compound!" << std::endl;
			CompoundInterpolatorOP interp( new CompoundInterpolator() );
			points_.push_back( Point(ubx_,uby_,ubdy_));
			Real lbx  = lbx_;
			//Real lby  = lby_;
			Real lbdy = lbdy_;
			for ( size_t i = 1; i <= points_.size(); i++ ) {
				Point & p( points_[i] );
				if ( p.has_dy ) {
					x.push_back(p.x);
					y.push_back(p.y);
					// std::cerr << "add range " << lbx << " " << p.x << std::endl;
					interp->add_range( InterpolatorOP( new SimpleInterpolator(x,y,lbdy,p.dy) ), lbx, p.x );
					lbx  = p.x;
					//lby  = p.y;  // set, but never used ~Labonte
					lbdy = p.dy;
					x.clear();
					y.clear();
				}
				x.push_back( p.x );
				y.push_back( p.y );
			}
			interpolator_ = InterpolatorOP(interp);
		}
		std::map<std::string, LinearFunction>::iterator lower_bound(boundary_functions_.find("lb_function"));
		std::map<std::string, LinearFunction>::iterator upper_bound(boundary_functions_.find("ub_function"));
		if ( lower_bound != boundary_functions_.end() ) {
			interpolator_->set_lb_function(lower_bound->second.cutoff,lower_bound->second.slope,lower_bound->second.intercept);
		}

		if ( upper_bound != boundary_functions_.end() ) {
			interpolator_->set_ub_function(upper_bound->second.cutoff,upper_bound->second.slope,upper_bound->second.intercept);
		}
	}
	return interpolator_;
}

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric
