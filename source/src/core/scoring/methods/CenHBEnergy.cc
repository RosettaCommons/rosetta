// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/CenHBEnergy.cc
/// @brief  Smooth, differentiable version of centroid hbond term
/// @author Frank DiMaio


#include <core/scoring/methods/CenHBEnergy.hh>
#include <core/scoring/methods/CenHBEnergyCreator.hh>
#include <core/scoring/CenHBPotential.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


#include <numeric/xyz.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/numeric.functions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyMap.hh>

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the CenHBEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CenHBEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new CenHBEnergy( ) );
}

ScoreTypes
CenHBEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cen_hb );
	return sts;
}


/// @details  C-TOR with method options object
CenHBEnergy::CenHBEnergy( ):
	parent( methods::EnergyMethodCreatorOP( new CenHBEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_CenHBPotential( ) ),
	soft_( false )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ corrections::score::hb_cen_soft ].user() ) {
		soft_ = option[ corrections::score::hb_cen_soft ]();
	}
}

/// clone
EnergyMethodOP
CenHBEnergy::clone() const {
	return EnergyMethodOP( new CenHBEnergy( *this ) );
}

/// @details  copy c-tor
CenHBEnergy::CenHBEnergy( CenHBEnergy const & src ):
	parent( src ),
	potential_( src.potential_ ),
	soft_( src.soft_ )
{ }


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
CenHBEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	pose.update_residue_neighbors();
}


void
CenHBEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {
	pose.update_residue_neighbors();
}


void
CenHBEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	Real score(0.0);

	//std::cout << "enter soft, " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

	Size seqsep = rsd1.polymeric_sequence_distance( rsd2 );

	// is there a way to get these coords without string comparisons?
	Vector bbH, bbO, bbC, bbN;
	Real r, xd, xh;

	if ( soft_ ) {
		if ( (rsd1.aa() != core::chemical::aa_gly && rsd1.aa() != core::chemical::aa_pro ) &&
				(rsd2.aa() != core::chemical::aa_gly && rsd2.aa() != core::chemical::aa_pro ) &&
				seqsep > 7 ) {
			Vector a1 = rsd1.atom( rsd1.atom_index("CEN") ).xyz() - rsd1.atom( rsd1.atom_index("CA") ).xyz();
			Vector a2 = rsd2.atom( rsd2.atom_index("CEN") ).xyz() - rsd2.atom( rsd2.atom_index("CA") ).xyz();
			Vector b1 = rsd1.atom( rsd1.atom_index("C") ).xyz()   - rsd1.atom( rsd1.atom_index("N") ).xyz();
			Vector b2 = rsd2.atom( rsd2.atom_index("C") ).xyz()   - rsd2.atom( rsd2.atom_index("N") ).xyz();
			Vector dv  = rsd2.atom( rsd2.atom_index("CA") ).xyz()  - rsd1.atom( rsd1.atom_index("CA") ).xyz();
			//printf( "%3d %3d %6.2f %6.2f %6.2f %6.2f %6.2f\n", int(rsd1.seqpos()), int(rsd2.seqpos()),
			//    a1.length(), a2.length(), b1.length(), b2.length(), dv.length() );
			score = potential_.func_soft( a1, a2, b1, b2, dv );
		}

	} else {

		if ( rsd1.aa() != core::chemical::aa_pro ) {
			bbH = rsd1.atom( rsd1.atom_index("H") ).xyz();
			bbO = rsd2.atom( rsd2.atom_index("O") ).xyz();

			r =  bbH.distance( bbO );
			if ( r <= potential_.cutoff( seqsep ) ) {
				bbN = rsd1.atom( rsd1.atom_index("N") ).xyz();
				bbC = rsd2.atom( rsd2.atom_index("C") ).xyz();
				xd = numeric::angle_degrees( bbN,bbH,bbO );
				xh = numeric::angle_degrees( bbH,bbO,bbC );
				//std::cerr << " r,xd,xh = " << r << ","  << 180-xd << ","  << 180-xh << " -- " << potential_.func( seqsep, r,180-xd,180-xh ) << std::endl;
				score += potential_.func( seqsep, r,180-xd,180-xh );
			}
		}

		if ( rsd2.aa() != core::chemical::aa_pro ) {
			bbH = rsd2.atom( rsd2.atom_index("H") ).xyz();
			bbO = rsd1.atom( rsd1.atom_index("O") ).xyz();

			r =  bbH.distance( bbO );//, xd, xh;
			if ( r <= potential_.cutoff( seqsep ) ) {
				bbN = rsd2.atom( rsd2.atom_index("N") ).xyz();
				bbC = rsd1.atom( rsd1.atom_index("C") ).xyz();
				xd = numeric::angle_degrees( bbN,bbH,bbO );
				xh = numeric::angle_degrees( bbH,bbO,bbC );

				score += potential_.func( seqsep, r,180-xd,180-xh );
				//std::cerr << " r,xd,xh = " << r << ","  << 180-xd << ","  << 180-xh << " -- " << potential_.func( seqsep, r,180-xd,180-xh ) << std::endl;
			}
		}
	}

	emap[ cen_hb ] += score;

	//std::cout << " soft done! " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;
}


void
CenHBEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	// is there a way to get these coords without string comparisons?
	using numeric::constants::f::pi;
	using namespace numeric::deriv;
	Size idxH, idxO, idxC, idxN;
	Vector bbH, bbO, bbC, bbN;
	Real r, xd, xh;
	Size seqsep = rsd1.polymeric_sequence_distance( rsd2 );
	Real weight = weights[ cen_hb ];

	if ( soft_ ) {
		eval_residue_pair_derivatives_soft( rsd1, rsd2, weights,
			r1_atom_derivs, r2_atom_derivs );
		return;
	}

	if ( rsd1.aa() != core::chemical::aa_pro ) {
		idxH = rsd1.atom_index("H"); bbH = rsd1.atom( idxH ).xyz();
		idxO = rsd2.atom_index("O"); bbO = rsd2.atom( idxO ).xyz();
		r =  bbH.distance( bbO );
		if ( r <= potential_.cutoff( seqsep ) ) {
			idxN = rsd1.atom_index("N"); bbN = rsd1.atom( idxN ).xyz();
			idxC = rsd2.atom_index("C"); bbC = rsd2.atom( idxC ).xyz();
			xd = numeric::angle_degrees( bbN,bbH,bbO );
			xh = numeric::angle_degrees( bbH,bbO,bbC );

			Vector dfunc = potential_.dfunc( seqsep, r,180-xd,180-xh );

			// convert deg->radians + invert angular terms, scale by weight
			dfunc[0] *= weight;
			dfunc[1] *= -weight * 180/pi;
			dfunc[2] *= -weight * 180/pi;

			// compute gradients
			// distance
			Vector f1,f2;
			Real temp_dist, temp_ang;
			distance_f1_f2_deriv(bbH, bbO, temp_dist, f1, f2);
			r1_atom_derivs[ idxH ].f1() += dfunc[0] * f1;
			r1_atom_derivs[ idxH ].f2() += dfunc[0] * f2;
			r2_atom_derivs[ idxO ].f1() += -dfunc[0] * f1;
			r2_atom_derivs[ idxO ].f2() += -dfunc[0] * f2;

			/// N-H-O angle
			angle_p1_deriv(  bbO, bbH, bbN, temp_ang, f1, f2);
			r2_atom_derivs[ idxO ].f1() += dfunc[1] * f1;
			r2_atom_derivs[ idxO ].f2() += dfunc[1] * f2;
			angle_p1_deriv(  bbN, bbH, bbO, temp_ang, f1, f2);
			r1_atom_derivs[ idxN ].f1() += dfunc[1] * f1;
			r1_atom_derivs[ idxN ].f2() += dfunc[1] * f2;
			angle_p2_deriv(  bbN, bbH, bbO, temp_ang, f1, f2);
			r1_atom_derivs[ idxH ].f1() += dfunc[1] * f1;
			r1_atom_derivs[ idxH ].f2() += dfunc[1] * f2;

			/// H-O-C angle
			angle_p1_deriv(  bbH, bbO, bbC, temp_ang, f1, f2);
			r1_atom_derivs[ idxH ].f1() += dfunc[2] * f1;
			r1_atom_derivs[ idxH ].f2() += dfunc[2] * f2;
			angle_p1_deriv(  bbC, bbO, bbH, temp_ang, f1, f2);
			r2_atom_derivs[ idxC ].f1() += dfunc[2] * f1;
			r2_atom_derivs[ idxC ].f2() += dfunc[2] * f2;
			angle_p2_deriv(  bbC, bbO, bbH, temp_ang, f1, f2);
			r2_atom_derivs[ idxO ].f1() += dfunc[2] * f1;
			r2_atom_derivs[ idxO ].f2() += dfunc[2] * f2;
		}
	}

	if ( rsd2.aa() != core::chemical::aa_pro ) {
		idxH = rsd2.atom_index("H"); bbH = rsd2.atom( idxH ).xyz();
		idxO = rsd1.atom_index("O"); bbO = rsd1.atom( idxO ).xyz();
		r =  bbH.distance( bbO );
		if ( r <= potential_.cutoff( seqsep ) ) {
			idxN = rsd2.atom_index("N"); bbN = rsd2.atom( idxN ).xyz();
			idxC = rsd1.atom_index("C"); bbC = rsd1.atom( idxC ).xyz();
			xd = numeric::angle_degrees( bbN,bbH,bbO );
			xh = numeric::angle_degrees( bbH,bbO,bbC );

			Vector dfunc = potential_.dfunc( seqsep, r,180-xd,180-xh );
			// convert deg->radians + invert angular terms, scale by weight
			dfunc[0] *= weight;
			dfunc[1] *= -weight * 180/pi;
			dfunc[2] *= -weight * 180/pi;


			// compute gradients
			// distance
			Vector f1,f2;
			Real temp_dist, temp_ang;
			distance_f1_f2_deriv(bbH, bbO, temp_dist, f1, f2);
			r2_atom_derivs[ idxH ].f1() += dfunc[0] * f1;
			r2_atom_derivs[ idxH ].f2() += dfunc[0] * f2;
			r1_atom_derivs[ idxO ].f1() -= dfunc[0] * f1;
			r1_atom_derivs[ idxO ].f2() -= dfunc[0] * f2;

			/// N-H-O angle
			angle_p1_deriv(  bbO, bbH, bbN, temp_ang, f1, f2);
			r1_atom_derivs[ idxO ].f1() += dfunc[1] * f1;
			r1_atom_derivs[ idxO ].f2() += dfunc[1] * f2;
			angle_p1_deriv(  bbN, bbH, bbO, temp_ang, f1, f2);
			r2_atom_derivs[ idxN ].f1() += dfunc[1] * f1;
			r2_atom_derivs[ idxN ].f2() += dfunc[1] * f2;
			angle_p2_deriv(  bbN, bbH, bbO, temp_ang, f1, f2);
			r2_atom_derivs[ idxH ].f1() += dfunc[1] * f1;
			r2_atom_derivs[ idxH ].f2() += dfunc[1] * f2;

			/// H-O-C angle
			angle_p1_deriv(  bbH, bbO, bbC, temp_ang, f1, f2);
			r2_atom_derivs[ idxH ].f1() += dfunc[2] * f1;
			r2_atom_derivs[ idxH ].f2() += dfunc[2] * f2;
			angle_p1_deriv(  bbC, bbO, bbH, temp_ang, f1, f2);
			r1_atom_derivs[ idxC ].f1() += dfunc[2] * f1;
			r1_atom_derivs[ idxC ].f2() += dfunc[2] * f2;
			angle_p2_deriv(  bbC, bbO, bbH, temp_ang, f1, f2);
			r1_atom_derivs[ idxO ].f1() += dfunc[2] * f1;
			r1_atom_derivs[ idxO ].f2() += dfunc[2] * f2;

		}
	}
}

/// @brief CenHBEnergy distance cutoff
Distance
CenHBEnergy::atomic_interaction_cutoff() const {
	return 6.0;
}

/// @brief CenHBEnergy
void
CenHBEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const {}

core::Size
CenHBEnergy::version() const {
	return 1; // Initial versioning
}

void
CenHBEnergy::eval_residue_pair_derivatives_soft(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {


	Size seqsep = rsd1.polymeric_sequence_distance( rsd2 );
	Real weight = weights[ cen_hb ];

	if ( (rsd1.aa() != core::chemical::aa_gly && rsd1.aa() != core::chemical::aa_pro ) &&
			(rsd2.aa() != core::chemical::aa_gly && rsd2.aa() != core::chemical::aa_pro ) &&
			seqsep > 7 ) {
		Size idxN2 = rsd1.atom_index("N");
		Size idxC1 = rsd1.atom_index("C");
		Size idxA1 = rsd1.atom_index("CA");
		Size idxB1 = rsd1.atom_index("CEN");
		Size idxN1 = rsd2.atom_index("N");
		Size idxC2 = rsd2.atom_index("C");
		Size idxA2 = rsd2.atom_index("CA");
		Size idxB2 = rsd2.atom_index("CEN");

		Vector const &crdA1 = rsd1.atom( idxA1 ).xyz();
		Vector const &crdB1 = rsd1.atom( idxB1 ).xyz();
		Vector const &crdN1 = rsd1.atom( idxN1 ).xyz();
		Vector const &crdC1 = rsd1.atom( idxC1 ).xyz();
		Vector const &crdA2 = rsd2.atom( idxA2 ).xyz();
		Vector const &crdB2 = rsd2.atom( idxB2 ).xyz();
		Vector const &crdN2 = rsd2.atom( idxN2 ).xyz();
		Vector const &crdC2 = rsd2.atom( idxC2 ).xyz();

		Vector a1 = crdB1 - crdA1; //rsd1.atom( idxB1 ).xyz() - rsd1.atom( idxA1 ).xyz();
		Vector a2 = crdB2 - crdA2; //rsd2.atom( idxB2 ).xyz() - rsd2.atom( idxA2 ).xyz();
		Vector b1 = crdC1 - crdN1; //rsd1.atom( idxC1 ).xyz() - rsd1.atom( idxN1 ).xyz();
		Vector b2 = crdC2 - crdN2; //rsd2.atom( idxC2 ).xyz() - rsd2.atom( idxN2 ).xyz();
		Vector dv = crdA2 - crdA1; //rsd2.atom( idxA2 ).xyz() - rsd1.atom( idxA1 ).xyz();

		utility::vector1< Vector > df_dABNC_1( 4, Vector(0.0) ); // in order of A1, B1, N1, C1
		utility::vector1< Vector > df_dABNC_2( 4, Vector(0.0) ); // in order of A2, B2, N2, C2

		potential_.dfunc_soft( a1, a2, b1, b2, dv, df_dABNC_1, df_dABNC_2 );

		for ( Size i = 1; i <= 4; ++i ) df_dABNC_1[i] *= weight;
		for ( Size i = 1; i <= 4; ++i ) df_dABNC_2[i] *= weight;

		r1_atom_derivs[ idxA1 ].f2() += df_dABNC_1[1];
		r1_atom_derivs[ idxB1 ].f2() += df_dABNC_1[2];
		r1_atom_derivs[ idxN1 ].f2() += df_dABNC_1[3];
		r1_atom_derivs[ idxC1 ].f2() += df_dABNC_1[4];
		r2_atom_derivs[ idxA2 ].f2() += df_dABNC_2[1];
		r2_atom_derivs[ idxB2 ].f2() += df_dABNC_2[2];
		r2_atom_derivs[ idxN2 ].f2() += df_dABNC_2[3];
		r2_atom_derivs[ idxC2 ].f2() += df_dABNC_2[4];

		r1_atom_derivs[ idxA1 ].f1() += crdA1.cross( -df_dABNC_1[1] + crdA1 );
		r1_atom_derivs[ idxB1 ].f1() += crdB1.cross( -df_dABNC_1[2] + crdB1 );
		r1_atom_derivs[ idxN1 ].f1() += crdN1.cross( -df_dABNC_1[3] + crdN1 );
		r1_atom_derivs[ idxC1 ].f1() += crdC1.cross( -df_dABNC_1[4] + crdC1 );
		r2_atom_derivs[ idxA2 ].f1() += crdA2.cross( -df_dABNC_2[1] + crdA2 );
		r2_atom_derivs[ idxB2 ].f1() += crdB2.cross( -df_dABNC_2[2] + crdB2 );
		r2_atom_derivs[ idxN2 ].f1() += crdN2.cross( -df_dABNC_2[3] + crdN2 );
		r2_atom_derivs[ idxC2 ].f1() += crdC2.cross( -df_dABNC_2[4] + crdC2 );
	}

}



}
}
}
