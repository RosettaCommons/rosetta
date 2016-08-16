// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/dna/DNABFormPotential.cc
/// @brief  DNA B-form specific torsion potential class implementation
/// @author Jim Havranek

// Unit Headers
#include <core/scoring/dna/DNABFormPotential.hh>

// Package Headers
#include <core/scoring/ScoreFunction.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
//#include <core/io/database/open.hh>
//#include <basic/options/option.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer tr( "core.scoring.dna.DNABFormPotential" );

using namespace ObjexxFCL;

namespace core {
namespace scoring {
namespace dna {

///////////////////////////////////////////////////////////////////////////////
///
Real
TorsionFourierComponent::compute( Real const torsion_angle, Real & return_deriv ) const
{
	using numeric::conversions::radians;

	Real score( 0.0 );

	// the torsion_angle and the phase should be in degrees
	Real const trig_term( radians( ( periodicity()*torsion_angle ) - phase() ) );

	score = factor()*( 1.0 + cos( trig_term ) );
	return_deriv = -1.0*factor()*periodicity()*sin( trig_term );

	return score;
}

///////////////////////////////////////////////////////////////////////////////


DNABFormPotential::DNABFormPotential()
{
	dummy_ = 0;
	init_dna_bform_data();
}

///////////////////////////////////////////////////////////////////////////////


void
DNABFormPotential::init_dna_bform_data()
{
	// Parameters for DNA backbone torsions are taken from the Amber MM code -
	// Alpha and gamma are from the refitting of Perez et.al. in Biophys. J. (2007)
	// v92 pp. 3816-3829.
	// The rest are taken from the parm99 parameter set available with the Ambertools
	// online.
	// It is my understanding that these should also work for RNA.
	// -jjh

	utility::vector1< TorsionFourierComponentCOP > alpha_components;
	alpha_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 0.185181, 1.0,  31.79508 ) ) );
	alpha_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 1.256531, 2.0, 351.95960 ) ) );
	alpha_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 0.354858, 3.0, 357.24748 ) ) );

	utility::vector1< TorsionFourierComponentCOP > beta_components;
	beta_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 1.150000, 3.0,   0.00000 ) ) );

	utility::vector1< TorsionFourierComponentCOP > gamma_components;
	gamma_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 1.178040, 1.0, 190.97653 ) ) );
	gamma_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 0.092102, 2.0, 295.63279 ) ) );
	gamma_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 0.962830, 3.0, 348.09535 ) ) );

	utility::vector1< TorsionFourierComponentCOP > delta_components;
	delta_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 1.400000, 3.0,   0.00000 ) ) );

	utility::vector1< TorsionFourierComponentCOP > epsilon_components;
	epsilon_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 1.150000, 3.0,   0.00000 ) ) );

	utility::vector1< TorsionFourierComponentCOP > zeta_components;
	zeta_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 1.200000, 2.0,   0.00000 ) ) );
	zeta_components.push_back( TorsionFourierComponentOP( new TorsionFourierComponent( 0.250000, 3.0,   0.00000 ) ) );

	bb_fourier_data.push_back(   alpha_components );
	bb_fourier_data.push_back(    beta_components );
	bb_fourier_data.push_back(   gamma_components );
	bb_fourier_data.push_back(   delta_components );
	bb_fourier_data.push_back( epsilon_components );
	bb_fourier_data.push_back(    zeta_components );

}

///////////////////////////////////////////////////////////////////////////////
///
void
DNABFormPotential::eval_dna_bform_bb_torsion_score_residue(
	conformation::Residue const & rsd,
	Real & score,
	Real & dscore_dchi,
	Size const torsion_id
) const
{
	using namespace numeric;

	debug_assert( rsd.is_DNA() );

	// Get the correct set of dihedral score components
	utility::vector1< TorsionFourierComponentCOP > const & this_data( bb_fourier_data[ torsion_id ] );

	Real total_score( 0.0 );
	Real total_deriv( 0.0 );

	// This loop sums over all the Fourier components for a dihedral's total score
	for ( Size icomp( 1 ) ; icomp <= this_data.size() ; ++icomp ) {
		Real this_deriv( 0.0 );
		total_score += (this_data[icomp])->compute( rsd.mainchain_torsion( torsion_id ), this_deriv );
		total_deriv += this_deriv;
	}

	score = total_score;
	dscore_dchi = total_deriv;
}

///////////////////////////////////////////////////////////////////////////////
///
void
DNABFormPotential::eval_dna_bform_chi_torsion_score_residue(
	conformation::Residue const & rsd,
	Real & score,
	Real & dscore_dchi
) const
{
	using namespace numeric;

	debug_assert( rsd.is_DNA() );

	Real chi = basic::unsigned_periodic_range( (rsd.chi())[1], 360.0 );

	// tr << "Found chi of " << chi << " for " << rsd.seqpos() << std::endl;

	// Very hacky for now - decide whether syn or anti, then pick a value based
	// on whether a purine or pyrimidine base.

	bool syn( false );
	if ( chi > 0.0 && chi < 90.0 ) {
		//  tr << "Using syn data " << std::endl;
		syn = true;
	} else {
		//  tr << "Using anti data " << std::endl;
		syn = false;
	}

	Real pot_min;
	Real pot_width;
	Real pot_depth( 1.0 );
	if ( rsd.aa() == core::chemical::na_ade || rsd.aa() == core::chemical::na_gua ) {
		// purines
		if ( syn ) {
			pot_min = 45.0;
			pot_width = 8.0;
		} else {
			pot_min = 258.0;
			pot_width = 14.0;
		}
	} else {
		// pyrimidines
		if ( syn ) {
			pot_min = 45.0;
			pot_width = 8.0;
		} else {
			pot_min = 241.0;
			pot_width = 8.0;
		}
	}

	score = (pot_depth*( chi - pot_min )*( chi - pot_min ))/( pot_width * pot_width );
	dscore_dchi = -1.0*(pot_depth/(pot_width*pot_width))*( chi - pot_min );
}


}
}
}
