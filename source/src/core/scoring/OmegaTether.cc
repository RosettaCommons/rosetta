// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <numeric/interpolation/spline/Cubic_spline.hh>

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
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>


namespace core {
namespace scoring {


OmegaTether::OmegaTether() {
	using namespace basic::options;
	use_phipsi_dep_ = option[ OptionKeys::corrections::score::bbdep_omega ]();
	if ( use_phipsi_dep_ ) {
		read_omega_tables( );
	}
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

	for ( int ii = 1; ii <= total_residue; ++ii ) {
		if ( !pose.residue(ii).is_protein()  || pose.residue(ii).is_terminus() || pose.residue(ii).is_virtual_residue() ) continue;

		Real omega_score,dscore_domega, dscore_dphi, dscore_dpsi ;
		eval_omega_score_residue(pose.residue(ii),omega_score,dscore_domega, dscore_dphi, dscore_dpsi );
		pose_energies.onebody_energies( ii )[omega] = omega_score;
	}
}

/// @brief Returns the mainchain torsion index corresponding to "phi".
/// @details Generally 1.  Set to 2 for beta-amino acids so that derivatives are calculated
/// for the dihedral two spaces before the peptide bond.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Size OmegaTether::phi_index( core::conformation::Residue const &rsd ) const
{
	if ( rsd.type().is_beta_aa() ) return 2; //Special case.
	return 1; //Default for alpha-amino acids.
}

/// @brief Returns the mainchain torsion index corresponding to "psi".
/// @details Generally 2 (alpha-amino acids) or 3 (beta-amino acids).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Size OmegaTether::psi_index( core::conformation::Residue const &rsd ) const
{
	if ( rsd.type().is_beta_aa() ) return 3; //Special case.
	return 2; //Default for alpha-amino acids.
}

/// @brief Returns the mainchain torsion index corresponding to "omega".
/// @details Should be 3 for alpha amino acids, 4 for beta amino acids.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::Size OmegaTether::omega_index( core::conformation::Residue const &rsd ) const
{
	if ( rsd.type().is_beta_aa() ) return 4;
	return 3; //Default for alpha-amino acids.
}


///////////////////////////////////////////////////////////////////////////////
///
Real
OmegaTether::eval_omega_score_residue(
	AA const res_aa,
	Real const omega,
	Real const phi,
	Real const psi
) const {
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
) const {
	using namespace numeric;

	debug_assert( rsd.is_protein() );

	bool const is_d( rsd.type().is_d_aa() );
	core::Real const d_multiplier( is_d ? -1.0 : 1.0 );
	core::chemical::AA the_aa( rsd.aa() );
	if ( !core::chemical::is_canonical_L_aa( the_aa ) /*includes gly*/ && !core::chemical::is_canonical_D_aa( the_aa ) ) {
		the_aa = core::chemical::aa_gly;
	} else if ( is_d ) {
		the_aa = core::chemical::get_L_equivalent( the_aa );
	}

	//Use backbone torsion angle 4 for omega if this is a beta-amino acid.
	Real const phi_angle
		( d_multiplier * nonnegative_principal_angle_degrees( rsd.mainchain_torsion( phi_index(rsd) )));
	Real const psi_angle
		( d_multiplier * nonnegative_principal_angle_degrees( rsd.mainchain_torsion( psi_index(rsd) )));
	Real const omega_angle
		( d_multiplier * nonnegative_principal_angle_degrees( rsd.mainchain_torsion( omega_index(rsd) )));
	if ( rsd.is_upper_terminus() || rsd.is_virtual_residue() ) { // begin or end of chain
		score = 0.0;
		dscore_domega = dscore_dphi = dscore_dpsi = 0.0;
		return;
	}

	eval_omega_score_residue( the_aa, omega_angle, phi_angle, psi_angle, score, dscore_domega, dscore_dphi, dscore_dpsi );

	dscore_dphi *= d_multiplier;
	dscore_dpsi *= d_multiplier;
	dscore_domega *= d_multiplier;
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
) const {
	using basic::subtract_degree_angles;

	core::Real omega_p = omega;
	while ( omega_p <  -90.0 ) omega_p += 360.0;
	while ( omega_p >  270.0 ) omega_p -= 360.0;

	if ( !use_phipsi_dep_ || omega_p < 90.0 ) {  // use standard form for cis omegas
		// standard form
		core::Real dangle;
		core::Real weight = 0.01;  // This is 1 in rosetta but divided by the number of residues oddly. we'll just assume N=100 here such
		// that omega can be calculated on a per residue basis

		if ( omega_p >= 90.0 ) {
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
		if ( aa == core::chemical::aa_gly ) {
			table = 2;
		}
		if ( aa == core::chemical::aa_pro || aa ==  core::chemical::aa_dpr ) {
			table = 3;
		}
		if ( aa == core::chemical::aa_ile || aa == core::chemical::aa_val || aa == core::chemical::aa_dil || aa == core::chemical::aa_dva ) {
			table = 4;
		}

		// phi-psi dependent
		// interpolate phi/psi
		core::Real mu, sigma, dmu_dphi, dmu_dpsi, dsigma_dphi, dsigma_dpsi;

		// D AAs
		core::Real phi_eff=phi, psi_eff=psi, scale=1.0;
		if ( core::chemical::is_canonical_D_aa(aa) ) {
			phi_eff=360.0-phi;
			psi_eff=360.0-psi;
			scale=-1.0;
		}

		mu = scale*omega_mus_all_splines_[table].F(phi_eff,psi_eff);
		dmu_dphi = scale*omega_mus_all_splines_[table].dFdx(phi_eff,psi_eff);
		dmu_dpsi = scale*omega_mus_all_splines_[table].dFdy(phi_eff,psi_eff);

		sigma = omega_sigmas_all_splines_[table].F(phi_eff,psi_eff);
		dsigma_dphi = omega_sigmas_all_splines_[table].dFdx(phi_eff,psi_eff);
		dsigma_dpsi = omega_sigmas_all_splines_[table].dFdy(phi_eff,psi_eff);

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


/// load bb-dep omega tables
void
OmegaTether::read_omega_tables() {
	omega_mus_all_.resize(4);
	omega_sigmas_all_.resize(4);
	omega_mus_all_splines_.resize(4);
	omega_sigmas_all_splines_.resize(4);

	for ( int i=1;  i<=4; ++i ) {
		utility::io::izstream stream;
		if ( i==1 )      basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.all.txt");
		else if ( i==2 ) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.gly.txt");
		else if ( i==3 ) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.pro.txt");
		else if ( i==4 ) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.valile.txt");
		//else if (i==5) basic::database::open( stream, "scoring/score_functions/omega/omega_ppdep.prepro.txt");

		read_table_from_stream( stream, omega_mus_all_[i], omega_sigmas_all_[i] );
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::symmetric_gly_tables].user() ) {
		//fpd this is tricky because the function needs to be symmetric around both phi/psi and omega
		//fpd there are two ways to do this:
		//fpd   1) set everything to 180
		//fpd   2) symmetrize here and make the potential symmetric about omega=180
		//fpd we do the former
		for ( core::Size iphi = 1; iphi <= 36; ++iphi ) {
			for ( core::Size ipsi = 1; ipsi <= 36; ++ipsi ) {
				omega_mus_all_[2](iphi,ipsi) = 180.0;
				omega_sigmas_all_[2](iphi,ipsi) = 6.0;
			}
		}
	}

	for ( int i=1;  i<=4; ++i ) {
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
