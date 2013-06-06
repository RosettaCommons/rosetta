// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/SplineInterp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>

#include <iostream>

// option key includes
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/patterson.OptionKeys.gen.hh>

#include <utility/vector1.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace electron_density {

/// @brief read density weights from the cmd line into the scorefunction
void add_dens_scores_from_cmdline_to_scorefxn( core::scoring::ScoreFunction &scorefxn ) {
	using namespace basic::options;

	if ( option[ OptionKeys::patterson::weight ].user() ) {
		scorefxn.set_weight( core::scoring::patterson_cc,
		                     option[ OptionKeys::patterson::weight ]() );
	}
	if ( option[ OptionKeys::edensity::fastdens_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_fast,
		                     option[ OptionKeys::edensity::fastdens_wt ]() );
	}
	if ( option[ OptionKeys::edensity::sliding_window_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_window,
		                     option[ OptionKeys::edensity::sliding_window_wt ]() );
	}
	if ( option[ OptionKeys::edensity::whole_structure_ca_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_whole_structure_ca,
		                     option[ OptionKeys::edensity::whole_structure_ca_wt ]() );
	}
	if ( option[ OptionKeys::edensity::whole_structure_allatom_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_whole_structure_allatom,
		                     option[ OptionKeys::edensity::whole_structure_allatom_wt ]() );
	}
}


/// @brief spline interpolation with periodic boundaries
core::Real interp_spline(
                         ObjexxFCL::FArray3D< double > & coeffs ,
                         numeric::xyzVector< core::Real > const & idxX ) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real retval = SplineInterp::interp3(&coeffs[0], dims, pt);
	return retval;
}

/// @brief spline interpolation with periodic boundaries
numeric::xyzVector<core::Real> interp_dspline(
                         ObjexxFCL::FArray3D< double > & coeffs ,
                         numeric::xyzVector< core::Real > const & idxX ) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real grad[3] = { 0,0,0 };
	SplineInterp::grad3(&grad[0], &coeffs[0], dims, pt);
	return numeric::xyzVector<core::Real>(grad[2],grad[1],grad[0]);
}

void spline_coeffs(
           ObjexxFCL::FArray3D< double > & data ,
           ObjexxFCL::FArray3D< double > & coeffs) {
	int dims[3] = { data.u3(), data.u2(), data.u1() };
	coeffs = data;
	SplineInterp::compute_coefficients3( &coeffs[0] , dims );
}

void spline_coeffs(
                         ObjexxFCL::FArray3D< float > & data ,
                         ObjexxFCL::FArray3D< double > & coeffs) {
	int N = data.u3()*data.u2()*data.u1();
	ObjexxFCL::FArray3D< double > data_d(data.u1(),data.u2(),data.u3()) ;
	for (int i=0; i<N; ++i)
		data_d[i] = (double)data[i];
	spline_coeffs( data_d, coeffs );
}


} // namespace constraints
} // namespace scoring
} // namespace core
