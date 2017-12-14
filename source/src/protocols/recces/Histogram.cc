// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/Histogram.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/Histogram.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.Histogram" );

namespace protocols {
namespace recces {

//Constructor
Histogram::Histogram( core::Real const min, core::Real const max, core::Real const spacing ):
	min_( min ),
	max_( max ),
	spacing_( spacing )
{
	runtime_assert( max > min );
	n_elem_ = static_cast<core::Size>( ( max - min ) / spacing ) + 1;
	for ( core::Size i = 1; i <= n_elem_; ++i ) hist_.push_back( 0 );
}

//Destructor
Histogram::~Histogram() = default;

void
Histogram::add( float const value, core::Size const n_items ) {
	core::Size bin_index;
	if ( value <= min_ ) {
		bin_index = 1;
	} else if ( value >= max_ ) {
		bin_index = n_elem_;
	} else {
		bin_index = static_cast<core::Size>( ( value - min_ ) / spacing_ ) + 1;
	}
	hist_[bin_index] += n_items;
}

void
Histogram::clear() {
	for ( core::Size i = 0; i <= n_elem_; ++i ) hist_[i] = 0;
}

utility::vector1<core::Real>
Histogram::get_scores() const {
	utility::vector1<core::Real> scores;
	for ( core::Size i = 1; i <= n_elem_; ++i ) {
		scores.push_back( min_ + spacing_ * ( i - 0.5 ) );
	}
	return scores;
}


} //recces
} //protocols
