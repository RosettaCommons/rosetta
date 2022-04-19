// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/energy_methods/FaMPSolvEnergy.cc
///
/// @brief  LK-Type Membrane Solvation Energy
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc
#define INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc

// Unit Headers
#include <core/energy_methods/FaMPSolvEnergy.hh>
#include <core/energy_methods/FaMPSolvEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>


// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>

// Utility headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyzVector.hh>

// C++ Headers

namespace core {
namespace energy_methods {

using namespace core::scoring::methods;

/// @brief Create Fresh Instance of the Energy Method
core::scoring::methods::EnergyMethodOP
FaMPSolvEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {

	return utility::pointer::make_shared< FaMPSolvEnergy >(
		( core::scoring::ScoringManager::get_instance()->etable( options ) ),
		( core::scoring::ScoringManager::get_instance()->memb_etable( options.etable_type() )),
		options.analytic_membetable_evaluation()
	);
}

core::scoring::ScoreTypes
FaMPSolvEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( FaMPSolv );
	return sts;
}


/// @brief Construct MP Solv energy from standard and membrane etable
FaMPSolvEnergy::FaMPSolvEnergy(
	core::scoring::etable::EtableCAP etable_in,
	core::scoring::etable::MembEtableCAP memb_etable_in,
	bool const analytic_membetable_evaluation
) :
	parent( utility::pointer::make_shared< FaMPSolvEnergyCreator >() ),
	etable_( etable_in ),
	memb_etable_( memb_etable_in ),
	// FIXME: move inside with lock() success check?
	solv1_( memb_etable_in.lock()->solv1() ),
	solv2_( memb_etable_in.lock()->solv2() ),
	dsolv1_( memb_etable_in.lock()->dsolv1() ),
	dsolv2_( memb_etable_in.lock()->dsolv2() ),
	dsolv_( etable_in.lock()->dsolv() ),
	memb_solv1_( memb_etable_in.lock()->memb_solv1() ),
	memb_solv2_( memb_etable_in.lock()->memb_solv2() ),
	memb_dsolv1_( memb_etable_in.lock()->memb_dsolv1() ),
	memb_dsolv2_( memb_etable_in.lock()->memb_dsolv2() ),
	//added to calculate fampsolv here instead of binning
	lk_dgfree_( memb_etable_in.lock()->lk_dgfree() ),
	memb_lk_dgfree_( memb_etable_in.lock()->memb_lk_dgfree() ),
	lj_radius_( memb_etable_in.lock()->lj_radius() ),
	lk_volume_( memb_etable_in.lock()->lk_volume() ),
	lk_lambda_( memb_etable_in.lock()->lk_lambda() ),
	//
	safe_max_dis2_( etable_in.lock()->get_safe_max_dis2() ),
	get_bins_per_A2_( etable_in.lock()->get_bins_per_A2() ),
	max_dis_(etable_.lock()->max_dis()),
	max_normal_dis_( max_dis_ - 1.5 ),
	analytic_etable_evaluation_(analytic_membetable_evaluation)
	//verbose_( false )
{
	// core::scoring::etable::MembEtableCOP memb_etable( memb_etable_in );
}


/// @brief Clone Energy Method
core::scoring::methods::EnergyMethodOP
FaMPSolvEnergy::clone() const {
	return utility::pointer::make_shared< FaMPSolvEnergy >( *this );
}

/// @brief Setup Energy Method for Derivatives
void
FaMPSolvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & scfxn
) const {

	init( pose );
	pose.update_residue_neighbors();
	fa_weight_ = scfxn.weights()[ core::scoring::FaMPSolv ];
}

/// @brief Evaluate Derivatives
/// @details Called during graident-based minimization inside dfunc
/// note: f1 and f2 are not zeroed - contributions are summed
void
FaMPSolvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {

	// Grab residue and atom numbers
	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	if ( m > rsd1.nheavyatoms() ) return;

	Vector const heavy_atom_i( rsd1.xyz( m ) );
	bool const pos1_fixed( domain_map( i ) != 0 );

	// Grab neighbor graphs from energies object
	core::scoring::Energies const & energies( pose.energies() );
	core::scoring::EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( utility::graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( (*iter)->get_other_ind( i ) );

		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue;

		// Grab second residue and check if residues are the same
		conformation::Residue const & rsd2( pose.residue( j ) );
		bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

		using namespace core::scoring::etable::count_pair;
		CountPairFunctionOP cpfxn( nullptr );

		if ( same_res ) {
			cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
		} else {
			cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
		}

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0; Size path_dist(0);

			if ( ! cpfxn->count(m, n, cp_weight, path_dist ) ) continue;
			// Grab projection on xy plane from both atoms
			Real proj_m = fa_proj_[ rsd1.seqpos() ][ m ];
			Real proj_n = fa_proj_[ rsd2.seqpos() ][ n ];

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();
			//Vector const d_ij_norm = d_ij.normalized();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;
			Vector f1_ij( 0.0 ), f2_ij( 0.0 );

			Real const dE_dR_over_r
				( eval_dE_dR_over_r( rsd1.atom(m), rsd2.atom(n), weights, f1_ij, f2_ij, proj_m, proj_n ) );

			//now calculate f1 and f2 with respect to distance from membrane center
			//F1 and F2 are defined in H. Abe, W. Braun, T. Noguti, N. Go, Computers & Chemistry 1984

			//needed for derivative with respect to distance from membrane center
			Real solvE1( 0.0 ), solvE2( 0.0 ), membsolvE1( 0.0 ), membsolvE2( 0.0 );
			solvationE( rsd1.atom(m), rsd2.atom(n), d2, solvE1, solvE2, membsolvE1, membsolvE2);

			Real diff = membsolvE1 - solvE1;

			Vector f1( 0.0 ), f2( 0.0 );
			f1 = (dE_dR_over_r * f1_ij) + fa_f1_[ rsd1.seqpos() ][ m ]*diff;
			f2 = (dE_dR_over_r * f2_ij) + fa_f2_[ rsd1.seqpos() ][ m ]*diff;

			if ( same_res ) {
				F1 += 0.5 * fa_weight_ * cp_weight * f1;
				F2 += 0.5 * fa_weight_ * cp_weight * f2;

			} else {
				F1 += fa_weight_ * cp_weight * f1;
				F2 += fa_weight_ * cp_weight * f2;

			}


		}
	}
}


/// @brief Compute Residue Pair Energy
void
FaMPSolvEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {

	Real fa_mbsolv_score( 0.0 );
	get_residue_pair_energy( rsd1, rsd2, pose, fa_mbsolv_score);
	emap[ core::scoring::FaMPSolv ] += fa_mbsolv_score;
}

/// @brief Evaluate Intra-Residue Energies
void
FaMPSolvEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {

	Real fa_mbsolv_score( 0.0 );
	get_residue_pair_energy( rsd, rsd, pose, fa_mbsolv_score);
	emap[ core::scoring::FaMPSolv ] += fa_mbsolv_score;
}

/// @brief Specify Interaction Cutoff for computing pair energies
Distance
FaMPSolvEnergy::atomic_interaction_cutoff() const {
	core::scoring::etable::EtableCOP etable( etable_ );
	return etable->max_dis();
}

/// @brief Provide context graphs
void
FaMPSolvEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Energy Method for Scoring
void
FaMPSolvEnergy::setup_for_scoring(
	pose::Pose & pose, core::scoring::ScoreFunction const &
) const {

	// Initialize fullatom data
	init( pose );
}


/// @brief Finalize method after computing totals
void
FaMPSolvEnergy::finalize_total_energy(
	pose::Pose &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & // totals
) const {}

///// Helper Methods //////////////////////

/// @brief Compute Residue Pair Energies
void
FaMPSolvEnergy::get_residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	Real & fa_mbsolv_score
) const {

	// Check if interacting residues are the same
	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	Real score (0.0);

	// Initialize Count Pair Function
	using namespace core::scoring::etable::count_pair;
	CountPairFunctionOP cpfxn( nullptr );

	if ( same_res ) {
		cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
	} else {
		cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}

	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const heavy_atom_i( rsd1.xyz( i ) );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0; Size path_dist( 0 );
			if ( ! cpfxn->count( i, j, cp_weight, path_dist ) ) continue;

			// Grab proj from both atoms
			Real proj_i = fa_proj_[ rsd1.seqpos() ][ i ];
			Real proj_j = fa_proj_[ rsd2.seqpos() ][ j ];

			Vector const heavy_atom_j( rsd2.xyz( j ) );

			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			bool debug( false );

			score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), d2, proj_i, proj_j, debug );
			if ( same_res ) score *= 0.5;

			fa_mbsolv_score += score;
		}
	}
}

/// @brief Evaluate LK Energy
Real
FaMPSolvEnergy::eval_lk(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const & d2,
	Real const & f1,
	Real const & f2,
	bool &
) const {


	if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) return 0.0;

	// Initialize Variables
	Real score( 0.0 );

	if ( !analytic_etable_evaluation_ ) {


		Real const d2_bin = d2 * get_bins_per_A2_;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real frac = d2_bin - ( disbin - 1 );
		int const l1 = solv1_.index( disbin, atom2.type(), atom1.type() );
		int const l2 = l1 + 1;

		// Membrane specific solvation
		// solvation of atom1 based on its distance from the membrane center on the membrane normal
		Real e11 = f1 * solv1_[ l1 ] + (1 - f1) * memb_solv1_[ l1 ];
		Real e12 = f1 * solv1_[ l2 ] + (1 - f1) * memb_solv1_[ l2 ];

		//pba solvation of atom2 based on its distance from the membrane center on the membrane normal
		Real e21 = f2 * solv2_[ l1 ] + (1 - f2) * memb_solv2_[ l1 ];
		Real e22 = f2 * solv2_[ l2 ] + (1 - f2) * memb_solv2_[ l2 ];

		Real e1 = e11 + e21;
		Real e2 = e12 + e22;

		score = e1 + frac * ( e2 - e1 );

	} else {
		Real solve1(0.0), solve2(0.0), solvE1( 0.0 ), solvE2( 0.0 ), membsolvE1( 0.0 ), membsolvE2( 0.0 );
		solvationE( atom1, atom2, d2, solvE1, solvE2, membsolvE1, membsolvE2);

		solve1 = f1 * solvE1 + (1 - f1) * membsolvE1;
		solve2 = f2 * solvE2 + (1 - f2) * membsolvE2;
		score = solve1 + solve2;
	}

	return score;
}

/// @brief Compute Change in Energy over distance (for minimization)
//F1 and F2 are vectors needed for the derivatives
//f1 and f2 are the values for the membrane transition function for atom i and j
Real
FaMPSolvEnergy::eval_dE_dR_over_r(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	core::scoring::EnergyMap const &,
	Vector & F1,
	Vector & F2,
	Real const & f1,
	Real const & f2
) const {

	F1 = atom1.xyz().cross( atom2.xyz() );
	F2 = atom1.xyz() - atom2.xyz();

	Real d2 = atom1.xyz().distance_squared( atom2.xyz() );
	if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real(0.0) ) ) return 0.0;

	Real deriv;

	if ( !analytic_etable_evaluation_ ) {
		// bin by distance:
		Real const d2_bin = d2 * get_bins_per_A2_;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real frac = d2_bin - ( disbin - 1 );

		int const l1 = dsolv1_.index( disbin, atom1.type(), atom2.type()),
			l2 = l1 + 1;
		//f1 is transition function for atom i
		//f2 is transition function for atom j
		Real e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
		Real e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
		Real e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
		Real e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
		Real e1 = e11 + e21;
		Real e2 = e12 + e22;

		deriv = ( e1 + frac * ( e2 - e1 ) );
	} else {

		Real dsolve1( 0.0 ), dsolve2( 0.0 ), dmembsolve1( 0.0 ), dmembsolve2( 0.0 );
		dsolvationE( atom1, atom2, d2, dsolve1, dsolve2, dmembsolve1, dmembsolve2 );
		Real de1 = f1 * dsolve1 + (1 - f1) * dmembsolve1;
		Real de2 = f2 * dsolve2 + (1 - f2) * dmembsolve2;
		deriv = de1 + de2;

	}

	return deriv / std::sqrt( d2 );
	//divides by the distance to incorporate 1/|rb - ra|, instead of dividing F1 and F2 by distance
}

/// @brief Versioning
core::Size
FaMPSolvEnergy::version() const {
	return (core::Size)2.0;
}

void
FaMPSolvEnergy::init( pose::Pose & pose ) const {

	// Alloc appropriate Farrays
	setup_for_fullatom( pose );

	core::conformation::Conformation const & conf( pose.conformation() );
	core::conformation::membrane::MembraneGeometryCOP mp_geometry( conf.membrane_info()->membrane_geometry() );

	// For convenience - grab nres
	Real nres = pose.size();

	for ( Size i = 1; i <= nres; ++i ) {
		for ( Size j = 1, j_end = pose.residue( i ).nheavyatoms(); j <= j_end; ++j ) {

			Vector const xyz( pose.residue( i ).xyz( j ) );

			//Compute transition function
			fa_proj_[i][j] = mp_geometry->f_transition( conf, i, j );

			// Compute derivatives
			fa_f1_[i][j] = mp_geometry->f_transition_f1( conf, i, j );
			fa_f2_[i][j] = mp_geometry->f_transition_f2( conf, i, j );

		}
	}
}

//solvation of atom i and j in water and chex
void
FaMPSolvEnergy::solvationE(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real dis2,
	Real &solvE1,
	Real &solvE2,
	Real &membsolvE1,
	Real &membsolvE2
) const {
	Real solv1, solv2;
	core::scoring::etable::EtableCOP etable = etable_.lock();
	core::scoring::etable::MembEtableCOP memb_etable = memb_etable_.lock();
	Real min_dis2 = etable->min_dis2(); //If distance is less than min_dis then make distance=min_dis
	//If distance is greater than max dis + epsilon then we don't need to calucate solvation
	Real max_dis2 = etable->max_dis2();
	Real epsilon = etable->epsilon();


	if ( dis2 > max_dis2 + epsilon ) {
		return;
	}
	if ( dis2 < min_dis2 ) {
		dis2 = min_dis2;
	}

	//distance between atom i and j
	Real dis = std::sqrt(dis2);


	//pulling sigma from MembEtable
	Real lk_min_dis2sigma = etable->lk_min_dis2sigma();
	Real sigma = memb_etable->lj_sigma(atom1.type(), atom2.type());

	Real thresh_dis = lk_min_dis2sigma * sigma;

	Real thresh_dis_min = thresh_dis/2.0;


	//If dis > max_dis solvation = 0 and <thresh_dis_min score should be
	//what it is at thresh_dis_min. A cubic spline is used to ensure that the function remains continuously differentiable.
	if ( dis > max_normal_dis_ && dis < max_dis_ ) {

		Real start_damping = max_normal_dis_;
		Real t = (dis - start_damping) / (max_dis_ - start_damping);
		//pk1 is the value of solvation for atom1 at the max_normal dis
		//mk1 is the derivative at max_normal_dis
		//pk2 and mk2 are the same excpet for atom2
		//the value of solvation and derivative at max_dis  for atom 1 and 2 are 0
		//so those terms drop out
		Real pk1 = solv( atom1.type(), atom2.type(), start_damping );
		Real mk1 = solv_deriv( atom1, start_damping ) * pk1;

		Real pk2 = solv( atom2.type(), atom1.type(), start_damping );
		Real mk2 = solv_deriv( atom2, start_damping ) * pk2;

		Real a1 = mk1*(max_dis_ - start_damping) - ( -pk1 );
		Real b1 = -pk1;

		Real a2 = mk2*(max_dis_ - start_damping) - ( -pk2 );
		Real b2 = -pk2;

		solv1 = (1-t)*pk1 + t*(1-t)*((1-t)*a1 + t*b1);
		solv2 = (1-t)*pk2 + t*(1-t)*((1-t)*a2 + t*b2);

	} else if ( dis < thresh_dis && dis > thresh_dis_min ) {
		Real t = ( dis - thresh_dis_min ) / ( thresh_dis - thresh_dis_min );
		//pk11 is the value of solvation at thresh_dis_min for atom 1
		//pk12 is the value of solvation at thresh_dis for atom 1
		//mk12 is the value of the derivative at thresh_dis for atom 1
		//the value of the slope at thresh_dis_min is 0
		//pk21, pk22, mk22 are the same except it is for atom2

		Real pk11 = solv( atom1.type(), atom2.type(), thresh_dis_min );
		Real pk12 = solv( atom1.type(), atom2.type(), thresh_dis );
		Real mk12 = solv_deriv( atom1, thresh_dis ) * pk12;

		Real pk21 = solv( atom2.type(), atom1.type(), thresh_dis_min );
		Real pk22 = solv( atom2.type(), atom1.type(), thresh_dis );
		Real mk22 = solv_deriv( atom2, thresh_dis ) * pk22;

		solv1 = ( 2*std::pow(t,3) - 3*std::pow(t,2) + 1 )*pk11 + (-2*std::pow(t,3) + 3*std::pow(t,2))*pk12 + ( std::pow(t,3) - std::pow(t,2) )*( thresh_dis - thresh_dis_min )*mk12;
		solv2 = ( 2*std::pow(t,3) - 3*std::pow(t,2) + 1 )*pk21 + (-2*std::pow(t,3) + 3*std::pow(t,2))*pk22 + ( std::pow(t,3) - std::pow(t,2) )*( thresh_dis - thresh_dis_min )*mk22;

	} else if ( dis <= thresh_dis_min ) {
		//constant value of solvation at distance less than thresh_dis_min
		solv1 = solv( atom1.type(), atom2.type(), thresh_dis_min );
		solv2 = solv( atom2.type(), atom1.type(), thresh_dis_min );

	} else if ( dis >= max_dis_ ) {
		//solvaiton at distance greater than max_dis is 0
		solv1 = 0;
		solv2 = 0;
	} else {
		//the solvation for distances between thresh_dis and max_normal_dis
		solv1 = solv( atom1.type(), atom2.type(), dis );
		solv2 = solv( atom2.type(), atom1.type(), dis );
	}

	solvE1 = lk_dgfree_[ atom1.type() ] * solv1;
	solvE2 = lk_dgfree_[ atom2.type() ] * solv2;
	membsolvE1 = memb_lk_dgfree_[ atom1.type() ] * solv1;
	membsolvE2 = memb_lk_dgfree_[ atom2.type() ] * solv2;

}

//solvation function used in solvationE that is independent of membrane depth of atom
Real
FaMPSolvEnergy::solv(
	int atom1type,
	int atom2type,
	Real dis
) const {
	Real const k = { -0.089793561062583294 }; //inv_neg2_tms_pi_sqrt_pi
	//k is a combination of the constant terms in the equation from P. Barth, J. Schonbrun, D. Baker PNAS 2007 equation 2 in the SI Materials and Methods.
	return k * lk_volume_[ atom2type ] * solv_piece( atom1type, dis );
}

//A portion of the solvation calculated in solv that is only dependent on one atom
Real
FaMPSolvEnergy::solv_piece(
	int atom_type,
	Real d
) const {
	Real lambda = lk_lambda_[ atom_type ];
	Real denom = d * d * lambda;
	Real dis_rad = d - lj_radius_[ atom_type ];
	Real dis_rad_lambda = dis_rad/lambda;
	Real expo = dis_rad_lambda * dis_rad_lambda;
	Real solv = std::exp(-expo) / denom;
	return solv;
}

//solvation partial derivative wrt distance from atom i and j
void
FaMPSolvEnergy::dsolvationE(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real dis2,
	Real &dsolvE1,
	Real &dsolvE2,
	Real &dmembsolvE1,
	Real &dmembsolvE2
) const {
	core::scoring::etable::EtableCOP etable = etable_.lock();
	core::scoring::etable::MembEtableCOP memb_etable = memb_etable_.lock();
	Real min_dis2 = etable->min_dis2();
	Real max_dis2 = etable->max_dis2();
	Real epsilon = etable->epsilon();


	if ( dis2 > max_dis2 + epsilon ) {
		return;
	}
	if ( dis2 < min_dis2 ) {
		dis2 = min_dis2;
	}

	Real dis = std::sqrt(dis2);

	//see comments above in SolvationE for notes on sigma and thresh_dis
	Real lk_min_dis2sigma = etable->lk_min_dis2sigma();
	Real sigma = memb_etable->lj_sigma(atom1.type(), atom2.type());


	Real thresh_dis = lk_min_dis2sigma * sigma;
	Real thresh_dis_min = thresh_dis/2.0;

	if ( dis > max_normal_dis_ && dis < max_dis_ ) {
		Real start_damping = max_normal_dis_;
		Real t = (dis - start_damping) / (max_dis_ - start_damping);

		Real pk1 = solv( atom1.type(), atom2.type(), start_damping );
		Real mk1 = solv_deriv( atom1, start_damping ) * pk1;

		Real pk2 = solv( atom2.type(), atom1.type(), start_damping );
		Real mk2 = solv_deriv( atom2, start_damping ) * pk2;

		Real a1 = mk1*(max_dis_ - start_damping) - ( -pk1 );
		Real b1 = -pk1;

		Real a2 = mk2*(max_dis_ - start_damping) - ( -pk2 );
		Real b2 = -pk2;

		Real deriv1 = (-pk1 + (1 - 2*t)*((1-t)*a1 + t*b1) + (t - t*t)*(-a1 + b1))*(1/(max_dis_ - start_damping));
		Real deriv2 = (-pk2 + (1 - 2*t)*((1-t)*a2 + t*b2) + (t - t*t)*(-a2 + b2))*(1/(max_dis_ - start_damping));

		dsolvE1 = lk_dgfree_[ atom1.type() ] * deriv1;
		dsolvE2 = lk_dgfree_[ atom2.type() ] * deriv2;
		dmembsolvE1 = memb_lk_dgfree_[ atom1.type() ] * deriv1;
		dmembsolvE2 = memb_lk_dgfree_[ atom2.type() ] * deriv2;

	} else if ( dis < thresh_dis && dis > thresh_dis_min ) {
		Real t = ( dis - thresh_dis_min ) / ( thresh_dis - thresh_dis_min );

		Real pk11 = solv( atom1.type(), atom2.type(), thresh_dis_min );
		Real pk12 = solv( atom1.type(), atom2.type(), thresh_dis );
		Real mk12 = solv_deriv( atom1, thresh_dis ) * pk12;

		Real pk21 = solv( atom2.type(), atom1.type(), thresh_dis_min );
		Real pk22 = solv( atom2.type(), atom1.type(), thresh_dis );
		Real mk22 = solv_deriv( atom2, thresh_dis ) * pk22;

		Real deriv1 = (( 6*std::pow(t,2) - 6*t )*pk11 + ( -6*std::pow(t,2) + 6*t )*pk12 + ( 3*std::pow(t,2) - 2*t )*( thresh_dis - thresh_dis_min )*mk12 )*(1/(thresh_dis-thresh_dis_min));
		Real deriv2 = (( 6*std::pow(t,2) - 6*t )*pk21 + ( -6*std::pow(t,2) + 6*t )*pk22 + ( 3*std::pow(t,2) - 2*t )*( thresh_dis - thresh_dis_min )*mk22 )*(1/(thresh_dis-thresh_dis_min));

		dsolvE1 = lk_dgfree_[ atom1.type() ] * deriv1;
		dsolvE2 = lk_dgfree_[ atom2.type() ] * deriv2;
		dmembsolvE1 = memb_lk_dgfree_[ atom1.type() ] * deriv1;
		dmembsolvE2 = memb_lk_dgfree_[ atom2.type() ] * deriv2;


	} else if ( dis <= thresh_dis_min || dis >= max_dis_ ) {
		dsolvE1 = 0.0;
		dsolvE1 = 0.0;
		dmembsolvE1 = 0.0;
		dmembsolvE1 = 0.0;

	} else {
		Real solv1 = solv( atom1.type(), atom2.type(), dis );
		Real solv2 = solv( atom2.type(), atom1.type(), dis );

		Real deriv1 = solv_deriv( atom1, dis );
		Real deriv2 = solv_deriv( atom2, dis );

		dsolvE1 = lk_dgfree_[ atom1.type() ] * solv1 * deriv1;
		dsolvE2 = lk_dgfree_[ atom2.type() ] * solv2 * deriv2;
		dmembsolvE1 = memb_lk_dgfree_[ atom1.type() ] * solv1 * deriv1;
		dmembsolvE2 = memb_lk_dgfree_[ atom2.type() ] * solv2 * deriv2;

	}
}

Real
FaMPSolvEnergy::solv_deriv(
	conformation::Atom const & atom,
	Real dis
) const {
	Real deriv = -2.0 * (((dis - lj_radius_[ atom.type() ]) * ( 1/std::pow(lk_lambda_[ atom.type() ], 2) )) + (1/dis));
	return deriv;
}

/// @brief Setup Data Members for Fullatom Info
void
FaMPSolvEnergy::setup_for_fullatom( pose::Pose & pose ) const {

	Real nres = pose.size();

	fa_proj_.resize( (Size)nres );
	fa_f1_.resize( (core::Size)nres );
	fa_f2_.resize( (core::Size)nres );


	static Size const MAX_AMINOACID_SIZE = 15;

	for ( Size i = 1; i <= nres; ++i ) {

		Size const max_size = std::max( MAX_AMINOACID_SIZE, pose.residue( i ).nheavyatoms() );

		fa_proj_[i].resize( max_size );
		fa_f1_[i].resize( max_size );
		fa_f2_[i].resize( max_size );

		for ( Size j = 1; j <= max_size; ++j ) {
			fa_proj_[i][j] = 0.0;
			fa_f1_[i][j].assign(0.0,0.0,0.0);
			fa_f2_[i][j].assign(0.0,0.0,0.0);


		}
	}
}

} // scoring
} // core

#endif // INCLUDED_core_energy_methods_FaMPSolvEnergy_cc
