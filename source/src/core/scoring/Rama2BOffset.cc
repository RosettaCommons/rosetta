// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Rama2BOffset.cc
/// @brief  
/// @author

// Unit Headers
#include <core/scoring/Rama2BOffset.hh>

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


namespace core {
namespace scoring {


Rama2BOffset::Rama2BOffset() {
	using namespace basic::options;
 	read_r2bo_tables( );
 	read_paapp_tables( );
}


///////////////////////////////////////////////////////////////////////////////
///
void
Rama2BOffset::eval_r2bo_rama_score(
	AA const res_aa1,
	AA const res_aa2,
	Real const psi1,
	Real const omega2,
	Real const phi2,
	Real & score,
	Real & dscore_dpsi1,
	Real & dscore_domega2,
	Real & dscore_dphi2
) const {
	using basic::subtract_degree_angles;

	core::Real omega_p = omega2;
	while( omega_p <  -90.0 ) omega_p += 360.0;
	while( omega_p >  270.0 ) omega_p -= 360.0;
	bool is_cis = (omega_p < 90.0);
	Size table = aapair_to_table_index(res_aa1, res_aa2, is_cis);

	// (1) phi-psi
	score = rama_[table].F(psi1,phi2);
	dscore_dpsi1 = rama_[table].dFdx(psi1,phi2);
	dscore_domega2 = 0.0;
	dscore_dphi2 = rama_[table].dFdy(psi1,phi2);
}

void
Rama2BOffset::eval_r2bo_omega_score(
	AA const res_aa1,
	AA const res_aa2,
	Real const psi1,
	Real const omega2,
	Real const phi2,
	Real & score,
	Real & dscore_dpsi1,
	Real & dscore_domega2,
	Real & dscore_dphi2
) const {
	using basic::subtract_degree_angles;

	core::Real omega_p = omega2;
	while( omega_p <  -90.0 ) omega_p += 360.0;
	while( omega_p >  270.0 ) omega_p -= 360.0;
	bool is_cis = (omega_p < 90.0);
	Size table = aapair_to_table_index(res_aa1, res_aa2, is_cis);

	core::Real mu, sigma, dmu_dphi2, dmu_dpsi1, dsigma_dphi2, dsigma_dpsi1;
	mu = omega_mu_[table].F(psi1,phi2);
	dmu_dpsi1 = omega_mu_[table].dFdx(psi1,phi2);
	dmu_dphi2 = omega_mu_[table].dFdy(psi1,phi2);
	sigma = omega_sig_[table].F(psi1,phi2);
	dsigma_dpsi1 = omega_sig_[table].dFdx(psi1,phi2);
	dsigma_dphi2 = omega_sig_[table].dFdy(psi1,phi2);

	core::Real normalization = log( 1/ (6* sqrt(2*numeric::constants::d::pi) ) );
	core::Real entropy = -log( 1/ (sigma* sqrt(2*numeric::constants::d::pi) ) );
	core::Real offset = subtract_degree_angles(omega_p, mu);
	core::Real logprob = offset*offset / (2*sigma*sigma);

	// omega score/deriv
	score = normalization + entropy + logprob;
	dscore_domega2 = offset / (sigma*sigma);

	// bb-dep omega score/deriv
	core::Real dscore_dmu = -dscore_domega2;
	core::Real dscore_dsigma = 1/(sigma) - 2*logprob/sigma;
	dscore_dpsi1 = dscore_dmu*dmu_dpsi1 + dscore_dsigma*dsigma_dpsi1;
	dscore_dphi2 = dscore_dmu*dmu_dphi2 + dscore_dsigma*dsigma_dphi2;
}


/// @brief Probability energies from P(aa|phi,psi): Low level calculation for non-terminus position
void
Rama2BOffset::eval_p_aa_pp_score(
	chemical::AA aa1, chemical::AA aa2, Real const psi1, Real const phi2, Real &score, Real &dE_dpsi1, Real &dE_dphi2 ) const
{
	core::Real subscore1 = paapp1_[(Size)aa1].F(psi1,phi2);
	core::Real sub1_dscore_dpsi1 = paapp1_[(Size)aa1].dFdx(psi1,phi2);
	core::Real sub1_dscore_dphi2 = paapp1_[(Size)aa1].dFdy(psi1,phi2);

	core::Real subscore2 = paapp2_[(Size)aa2].F(psi1,phi2);
	core::Real sub2_dscore_dpsi1 = paapp2_[(Size)aa2].dFdx(psi1,phi2);
	core::Real sub2_dscore_dphi2 = paapp2_[(Size)aa2].dFdy(psi1,phi2);

	score = 0.5*( subscore1+subscore2 );
	dE_dpsi1 = 0.5*(sub1_dscore_dpsi1+sub2_dscore_dpsi1);
	dE_dphi2 = 0.5*(sub1_dscore_dphi2+sub2_dscore_dphi2);

	//if (std::fabs(score) > 10000) {
	//	std::cerr << "PAAPP: " << score << " " << subscore1 << " " << subscore2 << " " << (Size)aa1 << " " << (Size)aa2 << " " << psi1 << " " << phi2 << std::endl;
	//}
}

Size
Rama2BOffset::aapair_to_table_index( chemical::AA const res_aa1, chemical::AA const res_aa2, bool cis ) const {
	using namespace core::chemical;
	if (cis) {
		if (res_aa2 == aa_pro) return CIS_XP;
		else return CIS_XX;
	}

	if (res_aa1 == aa_pro) {
		if (res_aa2 == aa_pro) return TRANS_PP;
		else if (res_aa2 == aa_gly) return TRANS_PG;
		else if (res_aa2 == aa_val || res_aa2 == aa_ile) return TRANS_PV;
		else return TRANS_PX;
	} else if (res_aa1 == aa_gly) {
		if (res_aa2 == aa_pro) return TRANS_GP;
		else if (res_aa2 == aa_gly) return TRANS_GG;
		else if (res_aa2 == aa_val || res_aa2 == aa_ile) return TRANS_GV;
		else return TRANS_GX;
	} else if (res_aa1 == aa_val || res_aa1 == aa_ile) {
		if (res_aa2 == aa_pro) return TRANS_VP;
		else if (res_aa2 == aa_gly) return TRANS_VG;
		else if (res_aa2 == aa_val || res_aa2 == aa_ile) return TRANS_VV;
		else return TRANS_VX;
	} else {
		if (res_aa2 == aa_pro) return TRANS_XP;
		else if (res_aa2 == aa_gly) return TRANS_XG;
		else if (res_aa2 == aa_val || res_aa2 == aa_ile) return TRANS_XV;
		else return TRANS_XX;
	}
}


/// load tables
void
Rama2BOffset::read_r2bo_tables( ) {
	rama_.resize(NRAMATABLES);
	omega_mu_.resize(NRAMATABLES);
	omega_sig_.resize(NRAMATABLES);

	ObjexxFCL::FArray2D< core::Real > rama_i, omega_mu_i, omega_sig_i;

	for (int i=1; i<=NRAMATABLES; ++i) {
		utility::io::izstream stream;
		switch (i) {
		case TRANS_XX:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_XX.txt"); break;
		case TRANS_XG:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_XG.txt"); break;
		case TRANS_XP:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_XP.txt"); break;
		case TRANS_XV:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_XV.txt"); break;
		case TRANS_GX:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_GX.txt"); break;
		case TRANS_GG:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_GG.txt"); break;
		case TRANS_GP:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_GP.txt"); break;
		case TRANS_GV:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_GV.txt"); break;
		case TRANS_PX:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_PX.txt"); break;
		case TRANS_PG:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_PG.txt"); break;
		case TRANS_PP:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_PP.txt"); break;
		case TRANS_PV:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_PV.txt"); break;
		case TRANS_VX:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_VX.txt"); break;
		case TRANS_VG:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_VG.txt"); break;
		case TRANS_VP:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_VP.txt"); break;
		case TRANS_VV:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_VV.txt"); break;
		case CIS_XP:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_cis_XP.txt"); break;
		case CIS_XX:
			basic::database::open( stream, "scoring/score_functions/rama/offset/rama_cis_XX.txt"); break;
		}

		read_table_from_stream( stream, rama_i, omega_mu_i, omega_sig_i );
		setup_interpolation( rama_i, rama_[i] );
		setup_interpolation( omega_mu_i, omega_mu_[i] );
		setup_interpolation( omega_sig_i, omega_sig_[i] );
	}
}

void
Rama2BOffset::read_paapp_tables()
{
	paapp1_.resize(chemical::num_canonical_aas);
	paapp2_.resize(chemical::num_canonical_aas);

	{
		utility::vector1 <ObjexxFCL::FArray2D<Real> > paas;
		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/P_AA_pp/P_AA_pp_offset1");
		read_paapp_table_from_stream( stream, paas );
		for (int i=1; i<=chemical::num_canonical_aas; ++i) {
			setup_interpolation( paas[i], paapp1_[i] );
		}
	}

	{
		utility::vector1 <ObjexxFCL::FArray2D<Real> > paas;
		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/P_AA_pp/P_AA_pp_offset2");
		read_paapp_table_from_stream( stream, paas );
		for (int i=1; i<=chemical::num_canonical_aas; ++i) {
			setup_interpolation( paas[i], paapp2_[i] );
		}
	}
}


void
Rama2BOffset::setup_interpolation(
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
	Real start_vals[2] = {5.0, 5.0}; // grid is shifted by five degrees.
	Real deltas[2] = {10.0, 10.0}; // grid is 10 degrees wide
	bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
	std::pair< Real, Real > unused[2];
	unused[0] = std::make_pair( 0.0, 0.0 );
	unused[1] = std::make_pair( 0.0, 0.0 );
	sx.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
}

void
Rama2BOffset::read_table_from_stream(
	utility::io::izstream &stream,
	ObjexxFCL::FArray2D< Real > &ramas,
	ObjexxFCL::FArray2D< Real > &om_mus,
	ObjexxFCL::FArray2D< Real > &om_sigmas
) {
	core::Size i, j;
	core::Real phi, psi, rama, mu,sig;
	std::string line;

	ramas.dimension(36,36);
	om_mus.dimension(36,36);
	om_sigmas.dimension(36,36);

	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		l >> i >> j >> psi >> phi >> rama >> mu >> sig;
		ramas(i+1,j+1) = rama;
		om_mus(i+1,j+1) = mu;
		om_sigmas(i+1,j+1) = sig;
	}
}

void
Rama2BOffset::read_paapp_table_from_stream(
	utility::io::izstream &stream,
	utility::vector1 <ObjexxFCL::FArray2D< Real > > &paas
) {
	core::Size i, j;
	core::Real prob_n, phi,psi;
	std::string line;

	paas.resize( chemical::num_canonical_aas );
	for (core::Size z=1; z<=paas.size(); ++z) {
		paas[z].dimension(36,36);
	}

	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		l >> i >> j >> psi >> phi;
		for (core::Size z=1; z<=paas.size(); ++z ) {
			l >> prob_n;
			paas[z](i+1,j+1) = prob_n;
		}
	}
}


}
}
