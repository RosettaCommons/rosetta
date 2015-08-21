// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/model_quality/rms_obj.cc
/// @brief  Small object to represent rms data.
/// @author James Thompson
/// @date   Thu Jan 10 06:55:37 2008


#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/RmsData.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>

namespace numeric {
namespace model_quality {

// begin singleton class stuff!
RmsData* RmsData::pinstance_ = 0;// initialize pointer
RmsData* RmsData::instance ()
{
	if ( pinstance_ == 0 ) {  // is it the first call?
		pinstance_ = new RmsData; // create sole instance
	}
	return pinstance_; // address of sole instance
}

/// @brief set up RmsData with default values
RmsData::RmsData()
{
	clear_rms();
}

void RmsData::add_rms(
	int i,
	ObjexxFCL::FArray2A< double > xp,
	ObjexxFCL::FArray2A< double > xe
) {

	xp.dimension( 3, i );
	xe.dimension( 3, i );


	++count_;
	for ( int j = 1; j <= 3; ++j ) {
		for ( int k = 1; k <= 3; ++k ) {
			xm_(k,j) += xp(j,i) * xe(k,i); // flopped
		}
	}

	xre_ += ( xe(1,i) * xe(1,i) ) + ( xe(2,i) * xe(2,i) ) + ( xe(3,i) * xe(3,i) );
	xrp_ += ( xp(1,i) * xp(1,i) ) + ( xp(2,i) * xp(2,i) ) + ( xp(3,i) * xp(3,i) );

	xse_(1) += xe(1,i);
	xse_(2) += xe(2,i);
	xse_(3) += xe(3,i);
	xsp_(1) += xp(1,i);
	xsp_(2) += xp(2,i);
	xsp_(3) += xp(3,i);


} // add_rms

void
RmsData::clear_rms()
{

	xm_ .dimension( 3, 3 ); // note zero origin
	xse_.dimension( 3 );
	xsp_.dimension( 3 );


	// init array we will making into a running sum
	for ( int j = 1; j <= 3; ++j ) {
		for ( int k = 1; k <= 3; ++k ) {
			xm_(j,k) = 0.0;
		}
	}
	xre_ = 0.0;
	xrp_ = 0.0;

	xse_(1) = 0.0;
	xse_(2) = 0.0;
	xse_(3) = 0.0;
	xsp_(1) = 0.0;
	xsp_(2) = 0.0;
	xsp_(3) = 0.0;
	count_  = 0;
}

} // end namespace model_quality
} // end namespace numeric
