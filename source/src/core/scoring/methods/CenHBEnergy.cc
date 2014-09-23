// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenHBEnergy.cc
/// @brief  Smooth, differentiable version of centroid hbond term
/// @author Frank DiMaio


#include <core/scoring/methods/CenHBEnergy.hh>
#include <core/scoring/methods/CenHBEnergyCreator.hh>
#include <core/scoring/CenHBPotential.hh>

// AUTO-REMOVED #include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairFunction.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/CountPairFactory.hh>
// AUTO-REMOVED #include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>


#include <numeric/xyz.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>
// AUTO-REMOVED #include <numeric/deriv/dihedral_deriv.hh>
#include <numeric/numeric.functions.hh>

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
	potential_( ScoringManager::get_instance()->get_CenHBPotential( ) ) { }


/// clone
EnergyMethodOP
CenHBEnergy::clone() const {
	return EnergyMethodOP( new CenHBEnergy( *this ) );
}

/// @details  copy c-tor
CenHBEnergy::CenHBEnergy( CenHBEnergy const & src ):
	parent( src ),
	potential_( src.potential_ ) { }


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
CenHBEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	pose.update_residue_neighbors();
}


///
void
CenHBEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {
	pose.update_residue_neighbors();
}


///
void
CenHBEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	Real score(0.0);

	Size seqsep = rsd1.polymeric_sequence_distance( rsd2 );

	// is there a way to get these coords without string comparisons?
	Vector bbH, bbO, bbC, bbN;
	Real r, xd, xh;

	if (rsd1.aa() != core::chemical::aa_pro) {
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

	if (rsd2.aa() != core::chemical::aa_pro) {
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

	emap[ cen_hb ] += score;
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

	if (rsd1.aa() != core::chemical::aa_pro) {
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

	if (rsd2.aa() != core::chemical::aa_pro) {
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



}
}
}
