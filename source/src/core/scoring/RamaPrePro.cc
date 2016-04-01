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
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) Feb 2016 -- made this compatible with canonical D-amino acids; returns 0 for noncanonicals.

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

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/algorithm/string.hpp>


namespace core {
namespace scoring {

static basic::Tracer TR("core.scoring.RamaPrePro");

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
) const {
	bool const is_d( core::chemical::is_canonical_D_aa(res_aa1) );
	core::Real const d_multiplier( is_d ? -1.0 : 1.0 );

	//If this is neither a canonical D-amino acid, nor a canonical L-amino acid, nor glycine return 0:
	if ( !core::chemical::is_canonical_L_aa( res_aa1 ) && !is_d && res_aa1 != core::chemical::aa_gly ) {
		score_rama = 0.0;
		denergy_dphi = 0.0;
		denergy_dpsi = 0.0;
	}

	//Get the L-equivalent if this is a canonical D-residue:
	core::chemical::AA const res_aa1_copy( is_d ? core::chemical::get_L_equivalent(res_aa1) : res_aa1 );
	core::Real const phi_copy( is_d ? -1.0*phi : phi );
	core::Real const psi_copy( is_d ? -1.0*psi : psi );

	if ( res_aa1_copy > core::chemical::num_canonical_aas ) { //Noncanonical case: return 0.
		score_rama = 0.0;
		denergy_dphi = 0.0;
		denergy_dpsi = 0.0;
	} else { //Canonical case: return something
		if ( res_aa2 == core::chemical::aa_pro || res_aa2 == core::chemical::aa_dpr ) { //VKM -- crude approximation: this residue is considered "pre-pro" if it precedes an L- or D-proline.  (The N and CD are achiral).
			score_rama = rama_pp_splines_[res_aa1_copy].F(phi_copy,psi_copy);
			denergy_dphi = d_multiplier * rama_pp_splines_[res_aa1_copy].dFdx(phi_copy,psi_copy);
			denergy_dpsi = d_multiplier * rama_pp_splines_[res_aa1_copy].dFdy(phi_copy,psi_copy);
		} else {
			score_rama = rama_splines_[res_aa1_copy].F(phi_copy,psi_copy);
			denergy_dphi = d_multiplier * rama_splines_[res_aa1_copy].dFdx(phi_copy,psi_copy);
			denergy_dpsi = d_multiplier * rama_splines_[res_aa1_copy].dFdy(phi_copy,psi_copy);
		}
	}
}

/// load tables

void
RamaPrePro::read_rpp_tables( ) {
	bool const symmetrize_gly( basic::options::option[ basic::options::OptionKeys::score::symmetric_gly_tables ]() ); //Should the gly tables be symmetrized?

	rama_splines_.resize(20);
	rama_pp_splines_.resize(20);

	// allocate space for raw data
	utility::vector1<  ObjexxFCL::FArray2D< Real > > data(20);
	for ( int i=1; i<=20; ++i ) {
		data[i].dimension(36,36);
	}

	///fpd hardcode for now
	read_rama_map_file_shapovalov("scoring/score_functions/rama/fd/all.ramaProb", data, symmetrize_gly);
	for ( int i=1; i<=20; ++i ) {
		setup_interpolation( data[i], rama_splines_[i], (i == static_cast<int>( core::chemical::aa_gly ) && symmetrize_gly) ? 5.0 : 0.0  ); //VKM: To symmetrize the gly tables, I need a five degree offset.  Question: should everything be offset by five degrees?
	}

	read_rama_map_file_shapovalov("scoring/score_functions/rama/fd/prepro.ramaProb", data, symmetrize_gly);
	for ( int i=1; i<=20; ++i ) {
		setup_interpolation( data[i], rama_pp_splines_[i], (i == static_cast<int>( core::chemical::aa_gly ) && symmetrize_gly) ? 5.0 : 0.0  ); //VKM: To symmetrize the gly tables, I need a five degree offset.  Question: should everything be offset by five degrees?
	}
}


/// @brief Adapted from Max's code, ramachandran.cc
/// @details If symmetrize_gly is true, the plot for glycine is made symmetric.
void
RamaPrePro::read_rama_map_file_shapovalov (
	std::string const &filename,
	utility::vector1<  ObjexxFCL::FArray2D< Real > > &data,
	bool const symmetrize_gly
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

		j = static_cast<int>( ceil(phi / 10.0 - 0.5) + 1 );
		k = static_cast<int>( ceil(psi / 10.0 - 0.5) + 1 );

		data[i](j,k) = prob;
		entropy[i] += prob * (-minusLogProb);
	} while (true);

	// Symmetrize the gly table, if we should:
	if ( symmetrize_gly ) {
		symmetrize_gly_table( data[static_cast<int>(core::chemical::aa_gly)], entropy[static_cast<int>(core::chemical::aa_gly)] );
	}

	// correct
	for ( int i=1; i<=20; ++i ) { //loop through all amino acids
		for ( int j=1; j<=36; ++j ) {
			for ( int k=1; k<=36; ++k ) {
				data[i](j,k) = -std::log( data[i](j,k) ) + entropy[i];
			}
		}
	}
}

void
RamaPrePro::setup_interpolation(
	ObjexxFCL::FArray2D< Real > & x,
	numeric::interpolation::spline::BicubicSpline  &sx,
	core::Real const &offset
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
	Real start_vals[2] = {offset, offset}; // fpd: shapovalov has 0 degree shift
	Real deltas[2] = {10.0, 10.0};   // grid is 10 degrees wide
	bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
	std::pair< Real, Real > unused[2];
	unused[0] = std::make_pair( 0.0, 0.0 );
	unused[1] = std::make_pair( 0.0, 0.0 );
	sx.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
}

/// @brief If the -symmetric_gly_tables option is used, symmetrize the aa_gly table.
/// @details By default, the gly table is asymmetric because it is based on statistics from the PDB (which disproportionately put glycine
/// in the D-amino acid region of Ramachandran space).  However, the intrinsic propensities of glycine make it equally inclined to favour
/// right- or left-handed conformation.  (Glycine is achrial, and can't have a preference.)  Must be called AFTER gly table load, but prior
/// to bicubic interpolation setup.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
RamaPrePro::symmetrize_gly_table(
	ObjexxFCL::FArray2D< core::Real > & data,
	core::Real &entropy
) const {
	TR << "Symmetrizing glycine RamaPrePro table." << std::endl;
	
	//The following is for debugging only:
	/*TR << "MATRIX_BEFORE:" << std::endl;
	for(core::Size j=1; j<=36; ++j) {
	for(core::Size k=1; k<=36; ++k) {
	TR << data(j,k) << "\t";
	}
	TR << std::endl;
	}
	TR << std::endl;*/

	entropy = 0;
	for ( core::Size j=1; j<=35; ++j ) {
		for ( core::Size k=j+1; k<=36; ++k ) { //Loop though one triangle of the data matrix
			core::Real const avg( (data(j,k) + data(37-j,37-k) ) / 2 );
			data(j,k) = avg;
			data(37-j,37-k) = avg;
			entropy += avg * (-2.0) * std::log(avg);
		}
	}

	for ( core::Size j=1; j<=18; ++j ) { //Loop through half of the diagonal
		core::Real const avg( ( data(j,j) + data(37-j,37-j) )/2 );
		data(j,j)=avg;
		data(37-j,37-j) = avg;
		entropy += avg * (-2.0) * std::log(avg);
	}

	//The following is for debugging only:
	/*TR << "MATRIX_AFTER:" << std::endl;
	for(core::Size j=1; j<=36; ++j) {
	for(core::Size k=1; k<=36; ++k) {
	TR << data(j,k) << "\t";
	}
	TR << std::endl;
	}*/
}

}
}
