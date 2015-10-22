// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/model_quality/maxsub.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details Routines for calculating maxsub-based structural quality scores. Based on code originally
/// written by Charlie Strauss for rosetta++, ported over by James Thompson.
///
/// @author James Thompson

#ifndef INCLUDED_numeric_model_quality_maxsub_hh
#define INCLUDED_numeric_model_quality_maxsub_hh


// ObjexxFCL Headers
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/RmsData.hh>

#include <utility/vector1.hh>

namespace numeric {
namespace model_quality {


/* void
maxsub_native(
FArray3_float const & x,
int & nali,
float & rms,
float & logeval
); */

/*
void
maxsub_partial(
const int nres,
FArray3A_float x1,
FArray3A_float x2,
FArray1A_bool occ1,
FArray1A_bool occ2,
int & nali,
float & rms,
float & logeval
); */


void
maxsub(
	int & nsup,
	ObjexxFCL::FArray1A_double xe,
	ObjexxFCL::FArray1A_double xp,
	double & rms,
	double & psi,
	int & nali,
	double & zscore,
	double & evalue,
	double & score,
	double rsmtol = 4.0,
	double distance_tolerance = 7.0
);


double
erfcc( double x );


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
//    Calculate the center of geometry for the selected atoms ---
///
/// @details
///
/// @param  C - [in/out]? -
/// @param  WT - [in/out]? -
/// @param  NAT - [in/out]? -
/// @param  XC - [in/out]? -
/// @param  YC - [in/out]? -
/// @param  ZC - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////

void
COMAS(
	ObjexxFCL::FArray1A< double > C,
	ObjexxFCL::FArray1A< double > WT,
	int NAT,
	double & XC,
	double & YC,
	double & ZC
);

} // namespace model_quality
} // namespace numeric

#endif
