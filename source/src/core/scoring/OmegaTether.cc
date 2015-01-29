// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/OmegaTether.cc
/// @brief  OmegaTether potential class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/scoring/OmegaTether.hh>

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


OmegaTether::OmegaTether() {
	using namespace basic::options;
	use_phipsi_dep_ = option[ OptionKeys::corrections::score::bbdep_omega ]();
	if (use_phipsi_dep_)
		read_omega_tables( );
}


///////////////////////////////////////////////////////////////////////////////
/// @brief evaluate omega score for each (protein) residue and store that score
/// in the pose.energies() object
void
OmegaTether::eval_omega_score_all(
	pose::Pose & pose,
	ScoreFunction const & scorefxn
) const
{
	if ( scorefxn.has_zero_weight( omega ) ) return; // unnecessary, righ?

	int const total_residue = pose.total_residue();

	Energies & pose_energies( pose.energies() );

	for ( int ii = 1; ii <= total_residue; ++ii )
	{
		if ( pose.residue(ii).is_protein()  && ! pose.residue(ii).is_terminus() && ! pose.residue(ii).is_virtual_residue()  )
		{
			Real omega_score,dscore_domega, dscore_dphi, dscore_dpsi ;
			eval_omega_score_residue(pose.residue(ii),omega_score,dscore_domega, dscore_dphi, dscore_dpsi );
			pose_energies.onebody_energies( ii )[omega] = omega_score;
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
///
Real
OmegaTether::eval_omega_score_residue(
	AA const res_aa,
	Real const omega,
	Real const phi,
	Real const psi
) const
{
	Real score, dscore_domega, dscore_dphi, dscore_dpsi ;
	eval_omega_score_residue( res_aa, omega, phi, psi, score, dscore_domega, dscore_dphi, dscore_dpsi  );
	return score;
}


///////////////////////////////////////////////////////////////////////////////
void
OmegaTether::eval_omega_score_residue(
	conformation::Residue const & rsd,
	Real & score,
	Real & dscore_domega,
	Real & dscore_dphi,
	Real & dscore_dpsi
) const
{
	using namespace numeric;

debug_assert( rsd.is_protein() );

	Real const phi_angle
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi_angle
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));
	Real const omega_angle
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(3 + ((rsd.has("CM") && rsd.mainchain_torsions().size()==4) ? 1 : 0) ))); //Use backbone torsion angle 4 for omega if this is a beta-amino acid.  Test for this by looking for a CM atom AND by checking the size of the mainchain_torsions vector (since methylated lysine has a CM).

	if ( rsd.is_upper_terminus() || rsd.is_virtual_residue() ) { // begin or end of chain
		score = 0.0;
		dscore_domega = dscore_dphi = dscore_dpsi = 0.0;
		return;
	}

	eval_omega_score_residue( rsd.aa(), omega_angle, phi_angle, psi_angle, score, dscore_domega, dscore_dphi, dscore_dpsi );
}


///////////////////////////////////////////////////////////////////////////////
///
void
OmegaTether::eval_omega_score_residue(
	AA const aa,
	Real const omega,
	Real const phi,
	Real const psi,
	Real & score,
	Real & dscore_domega,
	Real & dscore_dphi,
	Real & dscore_dpsi
) const
{
	using basic::subtract_degree_angles;

	core::Real omega_p = omega;
	while( omega_p <  -90.0 ) omega_p += 360.0;
	while( omega_p >  270.0 ) omega_p -= 360.0;
	
	if ( !use_phipsi_dep_ || omega_p < 90.0 ) {  // use standard form for cis omegas
		// standard form
		core::Real dangle;
		core::Real weight = 0.01;  // This is 1 in rosetta but divided by the number of residues oddly. we'll just assume N=100 here such
															 // that omega can be calculated on a per residue basis
		
		if( omega_p >= 90.0 ){
			// trans
			dangle = subtract_degree_angles(omega_p, 180);
		} else {
			// cis
			dangle = subtract_degree_angles(omega_p, 0);
		}

		score = weight*dangle*dangle;
		dscore_domega = weight*2*dangle;
		dscore_dphi = dscore_dpsi = 0.0;  // stdard form is independent of phi and psi
	} else {
		// figure out which table to use
		core::Size table=1;
		if (aa == core::chemical::aa_gly)
			table = 2;
		if (aa == core::chemical::aa_pro)
			table = 3;
		if (aa == core::chemical::aa_ile || aa == core::chemical::aa_val)
			table = 4;

		//fpd note: preproline is not yet implemented as this would have to be a 2b energy term

		// phi-psi dependent
		// interpolate phi/psi
		core::Real mu, sigma, dmu_dphi, dmu_dpsi, dsigma_dphi, dsigma_dpsi;

		mu = omega_mus_all_splines_[table].F(phi,psi);
		dmu_dphi = omega_mus_all_splines_[table].dFdx(phi,psi);
		dmu_dpsi = omega_mus_all_splines_[table].dFdy(phi,psi);

		sigma = omega_sigmas_all_splines_[table].F(phi,psi);
		dsigma_dphi = omega_sigmas_all_splines_[table].dFdx(phi,psi);
		dsigma_dpsi = omega_sigmas_all_splines_[table].dFdy(phi,psi);

		// energy
		core::Real normalization = log( 1/ (6* sqrt(2*numeric::constants::f::pi) ) );  // shift scores to match default
		core::Real entropy = -log( 1/ (sigma* sqrt(2*numeric::constants::f::pi) ) );
		core::Real offset = subtract_degree_angles(omega_p, mu);
		core::Real logprob = offset*offset / (2*sigma*sigma);
		score = normalization + entropy + logprob;
		//score = logprob;

		// derivatives
		dscore_domega = offset / (sigma*sigma);

		// derivatives w.r.t mu,sigma
		core::Real dscore_dmu = -dscore_domega;
		core::Real dscore_dsigma = 1/(sigma) - 2*logprob/sigma;

		// derivatives w.r.t. phi/psi
		dscore_dphi = dscore_dmu*dmu_dphi + dscore_dsigma*dsigma_dphi;
		dscore_dpsi = dscore_dmu*dmu_dpsi + dscore_dsigma*dsigma_dpsi;
	}
}

///
/// load bb-dep omega tables
void
OmegaTether::read_omega_tables( ) {
	omega_mus_all_.resize(5);
	omega_sigmas_all_.resize(5);
	omega_mus_all_splines_.resize(5);
	omega_sigmas_all_splines_.resize(5);

	for (int i=1;  i<=5; ++i) {
		utility::io::izstream stream;
		if (i==1)      basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.all.txt");
		else if (i==2) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.gly.txt");
		else if (i==3) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.pro.txt");
		else if (i==4) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.valile.txt");
		//else if (i==5) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.prepro.txt");
		else if (i==5) continue;  // prepro not yet used

		read_table_from_stream( stream, omega_mus_all_[i], omega_sigmas_all_[i] );
		setup_interpolation( omega_mus_all_[i], omega_mus_all_splines_[i] );
		setup_interpolation( omega_sigmas_all_[i], omega_sigmas_all_splines_[i] );
	}
}

void
OmegaTether::setup_interpolation(
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
OmegaTether::read_table_from_stream(
	utility::io::izstream &stream,
	ObjexxFCL::FArray2D< Real > &mus, 
	ObjexxFCL::FArray2D< Real > &sigmas
) {
	core::Size i, j;
	core::Real phi, psi,  mu,sig;
	std::string line;

	mus.dimension(36,36);
	sigmas.dimension(36,36);

	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		l >> i >> j >> phi >> psi >> mu >> sig;
		mus(i+1,j+1) = mu;
		sigmas(i+1,j+1) = sig;
	}
}


}
}
