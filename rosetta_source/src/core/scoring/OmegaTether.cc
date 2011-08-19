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
// AUTO-REMOVED #include <basic/database/open.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/basic.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <numeric/interpolation/periodic_range/half/interpolation.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2A.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray4D.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <utility/options/keys/BooleanOptionKey.hh>


namespace core {
namespace scoring {


OmegaTether::OmegaTether()
{
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
		if ( pose.residue(ii).is_protein()  && ! pose.residue(ii).is_terminus()  )
		{
			Real omega_score,dscore_domega;
			eval_omega_score_residue(pose.residue(ii),omega_score,dscore_domega);
			//std::cout << "Rama: residue " << ii << " = " << omega_score << std::endl;
			pose_energies.onebody_energies( ii )[omega] = omega_score;
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
void
OmegaTether::eval_omega_score_residue(
	conformation::Residue const & rsd,
	Real & score,
	Real & dscore_domega
) const
{
	using namespace numeric;

	//assert( pose.residue(res).is_protein() );
	assert( rsd.is_protein() );

	Real const omega_angle
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(3)));

	if ( rsd.is_upper_terminus() ) { // begin or end of chain
		score = 0.0;
		dscore_domega = 0.0;
		return;
	}

	eval_omega_score_residue( rsd.aa(), omega_angle, score, dscore_domega );
}


///////////////////////////////////////////////////////////////////////////////
///
Real
OmegaTether::eval_omega_score_residue(
	AA const res_aa,
	Real const omega
) const
{

	Real score, dscore_domega;
	eval_omega_score_residue( res_aa, omega, score, dscore_domega );
	return score;
}

///////////////////////////////////////////////////////////////////////////////
///
void
OmegaTether::eval_omega_score_residue(
	AA const,
	Real const omega,
	Real & score,
	Real & dscore_domega
) const
{

	using basic::subtract_degree_angles;

	core::Real dangle;
	core::Real weight = 0.01;  // This is 1 in rosetta but divided by the number of residues oddly. we'll just assume N=100 here such
												// that omega can be calculated on a per residue basis

	core::Real omega_p = omega;

	while( omega_p <  -90.0 ) omega_p += 360.0;
	while( omega_p >  270.0 ) omega_p -= 360.0;

	if( omega_p >= 90.0 ){
		// trans
		dangle = subtract_degree_angles(omega_p, 180);
	}else{
		// cis
		dangle = subtract_degree_angles(omega_p, 0);
	}

	score = weight*dangle*dangle;
	//std::cout << "OMEGA  "  << omega_p << "  " <<  dangle << "   " <<  weight * ( dangle * dangle ) << "  " << score << std::endl;

	dscore_domega = weight*2*dangle;
}




}
}
