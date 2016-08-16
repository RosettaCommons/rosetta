// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief interpolate n_replica number

#include <devel/replica_docking/TempInterpolator.hh>
#include <core/types.hh>
#include <cmath>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer tr( "devel.replica_docking.TempInterpolator" );

namespace devel {
namespace replica_docking {

TempFixValue::TempFixValue( core::Real value ) {
	value_ = value;
}

// Core::Real
// TempFixValue::get_value( core::Size temp_level ) {
//   return value_;
// }

TempFixValue::~TempFixValue(){}

TempInterpolator::TempInterpolator( core::Size n_levels, core::Real start, core::Real end, std::string curve ) {

	start_ = start;
	end_ = end;
	//  ascend_ = ascend;
	n_levels_ = n_levels;
	curve_ = curve;
	calculated_ = false;
	interpolate();
}

TempInterpolator::TempInterpolator( utility::tag::TagCOP tag, core::Size n_levels ) {
	if ( tag->hasOption( "start" ) && tag->hasOption( "end" ) ) {
		start_ = tag->getOption< core::Real >( "start" );
		end_ = tag->getOption< core::Real >( "end" );
	} else {
		utility_exit_with_message( "TempInterpolator requires at least start and end values" ) ;
	}
	curve_ = tag->getOption< std::string >( "curve", "expotential" );
	n_levels_ = n_levels;
	calculated_ = false;
	interpolate();
}
TempInterpolator::~TempInterpolator(){}

void TempInterpolator::interpolate() {
	interpolated_nums_.clear();
	//  core::Size const n_levels = tempering_->n_temp_levels();

	runtime_assert( n_levels_ >= 2 );
	tr.Info << "interpolate temp_level dependent parameters from " << start_ << " to " << end_ << " with " << n_levels_
		<< " levels using " << curve_ << " interpolation. " << std::endl;

	if ( curve_ == "linear" ) {
		core::Real const step ( ( end_ - start_)/( n_levels_ - 1 ) );
		for ( Size ct=0; ct<n_levels_; ++ct ) {
			interpolated_nums_.push_back( start_ + ct*step );
		}
		calculated_ = true;
	} else if ( curve_ == "exponential" ) {
		core::Real next( start_ );
		core::Real factor( pow(end_/start_, 1/core::Real(n_levels_ - 1 )) );
		for ( core::Size ct=1; ct<n_levels_; ++ct ) {
			interpolated_nums_.push_back( next );
			next *= factor;
		}
		interpolated_nums_.push_back( end_ );
		calculated_ = true;
	} else {
		utility_exit_with_message( "we don't have " + curve_ + " interpolation!" );
	}
}


core::Real TempInterpolator::get_value ( core::Size temp_level ) {
	if ( !calculated_ ) interpolate();
	return interpolated_nums_[ temp_level ];
}

}
}
