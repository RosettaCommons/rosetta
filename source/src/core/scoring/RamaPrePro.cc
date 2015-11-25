// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RamaPrePro.cc
/// @brief
/// @author

// Unit Headers
#include <core/scoring/RamaPrePro.hh>

// Package Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <basic/basic.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/MathMatrix.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/algorithm/string.hpp>


namespace core {
namespace scoring {


RamaPrePro::RamaPrePro() {
	using namespace basic::options;
	read_rpp_tables( );
}


///////////////////////////////////////////////////////////////////////////////
///
void
RamaPrePro::eval_rpp_rama_score(
	AA const res_aa1,
	AA const res_aa2,
	Real const phi,
	Real const psi,
	Real & score_rama,
	Real & denergy_dphi,
	Real & denergy_dpsi
) const
{
	if (res_aa2 == core::chemical::aa_pro) {
		score_rama = rama_pp_splines_[res_aa1].F(phi,psi);
		denergy_dphi = rama_pp_splines_[res_aa1].dFdx(phi,psi);
		denergy_dpsi = rama_pp_splines_[res_aa1].dFdy(phi,psi);
	} else {
		score_rama = rama_splines_[res_aa1].F(phi,psi);
		denergy_dphi = rama_splines_[res_aa1].dFdx(phi,psi);
		denergy_dpsi = rama_splines_[res_aa1].dFdy(phi,psi);
	}
}

/// load tables

void
RamaPrePro::read_rpp_tables( ) {
	rama_splines_.resize(20);
	rama_pp_splines_.resize(20);

	// allocate space for raw data
	utility::vector1<  ObjexxFCL::FArray2D< Real > > data(20);
	for (int i=1; i<=20; ++i) {
		data[i].dimension(36,36);
	}

	///fpd hardcode for now
	read_rama_map_file_shapovalov("scoring/score_functions/rama/fd/all.ramaProb", data);
	for (int i=1; i<=20; ++i) {
		setup_interpolation( data[i], rama_splines_[i]);
	}

	read_rama_map_file_shapovalov("scoring/score_functions/rama/fd/prepro.ramaProb", data);
	for (int i=1; i<=20; ++i) {
		setup_interpolation( data[i], rama_pp_splines_[i]);
	}
}


// adapted from Max's code, ramachandran.cc
void
RamaPrePro::read_rama_map_file_shapovalov (
	std::string filename,
	utility::vector1<  ObjexxFCL::FArray2D< Real > > &data
) {
	//ala     -180.0  -180.0  5.436719e-004   7.517165e+000
	utility::io::izstream  iunit;

	// search in the local directory first
	iunit.open( filename );

	if ( !iunit.good() ) {
		iunit.close();
		if ( !basic::database::open( iunit, filename ) ) {
			std::stringstream err_msg;
			err_msg << "Unable to open Ramachandran map '" << filename << "'.";
			utility_exit_with_message(err_msg.str());
		}
	}


	char line[256];

	int i; // aa index
	int j, k; //phi and psi indices
	char aa[4];
	std::string aaStr;
	double phi, psi, prob, minusLogProb;
	utility::vector1< core::Real > entropy(20,0);

	// read the file
	do {
		iunit.getline( line, 255 );

		if ( iunit.eof() ) break;
		if ( line[0]=='#' || line[0]=='\n' || line[0]=='\r' ) continue;

		std::sscanf( line, "%3s%lf%lf%lf%lf", aa, &phi, &psi, &prob, &minusLogProb );
		std::string prevAAStr = aaStr;
		//std::cout << line << " :: " << aaStr << std::endl;
		aaStr = std::string(aa);
		boost::to_upper(aaStr);
		i = core::chemical::aa_from_name(aaStr);

		if ( phi < 0 ) phi += 360;
		if ( psi < 0 ) psi += 360;

		j = (int) ceil(phi / 10.0 - 0.5) + 1;
		k = (int) ceil(psi / 10.0 - 0.5) + 1;

		data[i](j,k) = prob;
		entropy[i] += prob * (-minusLogProb);
	} while (true);

	// correct
	for (int i=1; i<=20; ++i) {
		for (int j=1; j<=36; ++j) {
			for (int k=1; k<=36; ++k) {
				data[i](j,k) = -std::log( data[i](j,k) ) + entropy[i];
			}
		}
	}
}

void
RamaPrePro::setup_interpolation(
	ObjexxFCL::FArray2D< Real > & x,
	numeric::interpolation::spline::BicubicSpline  &sx
) {
	using namespace numeric;
	using namespace numeric::interpolation::spline;

	// just interpolate phi/psi
	MathMatrix< Real > energy_vals( 36, 36 );
	for ( Size jj = 0; jj < 36; ++jj ) {
		for ( Size kk = 0; kk < 36; ++kk ) {
			energy_vals( jj, kk ) = x(jj+1,kk+1);
		}
	}
	BorderFlag periodic_boundary[2] = { e_Periodic, e_Periodic };
	Real start_vals[2] = {0.0, 0.0}; // fpd: shapovalov has 0 degree shift
	Real deltas[2] = {10.0, 10.0};   // grid is 10 degrees wide
	bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
	std::pair< Real, Real > unused[2];
	unused[0] = std::make_pair( 0.0, 0.0 );
	unused[1] = std::make_pair( 0.0, 0.0 );
	sx.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
}

}
}
