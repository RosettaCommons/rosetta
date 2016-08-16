// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/FaMPSolvEnergy.cc
///
/// @brief  LK-Type Membrane Solvation Energy
/// @details Last Modified: 5/13/14
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc
#define INCLUDED_core_scoring_membrane_FaMPSolvEnergy_cc

// Unit Headers
#include <core/scoring/membrane/FaMPSolvEnergy.hh>
#include <core/scoring/membrane/FaMPSolvEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>


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

// Utility headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

using namespace core::scoring::methods;

/// @brief Create Fresh Instance of the Energy Method
methods::EnergyMethodOP
FaMPSolvEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {

	return methods::EnergyMethodOP( new FaMPSolvEnergy(
		( ScoringManager::get_instance()->etable( options ) ),
		( ScoringManager::get_instance()->memb_etable( options.etable_type() ))
		) );
}

ScoreTypes
FaMPSolvEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( FaMPSolv );
	return sts;
}


/// @brief Construct MP Solv energy from standard and membrane etable
FaMPSolvEnergy::FaMPSolvEnergy(
	etable::EtableCAP etable_in,
	etable::MembEtableCAP memb_etable_in
) :
	parent( EnergyMethodCreatorOP( new FaMPSolvEnergyCreator ) ),
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
	safe_max_dis2_( etable_in.lock()->get_safe_max_dis2() ),
	get_bins_per_A2_( etable_in.lock()->get_bins_per_A2() ),
	verbose_( false )
{
	// etable::MembEtableCOP memb_etable( memb_etable_in );
}


/// @brief Clone Energy Method
EnergyMethodOP
FaMPSolvEnergy::clone() const {
	return EnergyMethodOP( new FaMPSolvEnergy( *this ) );
}

/// @brief Setup Energy Method for Derivatives
void
FaMPSolvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const & scfxn
) const {

	init( pose );
	pose.update_residue_neighbors();
	fa_weight_ = scfxn.weights()[ FaMPSolv ];
}

/// @brief Evaluate Derivatives
/// @details Called during graident-based minimization inside dfunc
/// note: f1 and f2 are not zeroed - contributions are summed
void
FaMPSolvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
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
	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );


	for ( graph::Graph::EdgeListConstIter
			iter  = energy_graph.get_node( i )->const_edge_list_begin(),
			itere = energy_graph.get_node( i )->const_edge_list_end();
			iter != itere; ++iter ) {

		Size const j( (*iter)->get_other_ind( i ) );

		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue;

		// Grab second residue and check if residues are the same
		conformation::Residue const & rsd2( pose.residue( j ) );
		bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

		using namespace etable::count_pair;
		CountPairFunctionOP cpfxn( 0 );

		if ( same_res ) {
			cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
		} else {
			cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
		}

		for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

			Real cp_weight = 1.0; Size path_dist(0);

			if ( ! cpfxn->count(m, n, cp_weight, path_dist ) ) continue;

			// Grab proj from both atoms
			core::Real proj_m = fa_proj_[ rsd1.seqpos() ][ m ];
			core::Real proj_n = fa_proj_[ rsd2.seqpos() ][ n ];

			Vector const heavy_atom_j( rsd2.xyz( n ) );
			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();
			Vector const d_ij_norm = d_ij.normalized();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			Vector f1( 0.0 ), f2( 0.0 );

			Real const dE_dR_over_r
				( eval_dE_dR_over_r( rsd1.atom(m), rsd2.atom(n), weights, f1, f2, proj_m, proj_n ) );
			if ( dE_dR_over_r == 0.0 ) continue;

			if ( same_res ) {
				F1 += 0.5 * dE_dR_over_r * cp_weight * f1;
				F2 += 0.5 * dE_dR_over_r * cp_weight * f2;
			} else {
				F1 += dE_dR_over_r * cp_weight * f1;
				F2 += dE_dR_over_r * cp_weight * f2;
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
	ScoreFunction const &,
	EnergyMap & emap
) const {

	Real fa_mbsolv_score( 0.0 );
	get_residue_pair_energy( rsd1, rsd2, pose, fa_mbsolv_score);
	emap[ FaMPSolv ] += fa_mbsolv_score;
}

/// @brief Evaluate Intra-Residue Energies
void
FaMPSolvEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	Real fa_mbsolv_score( 0.0 );
	get_residue_pair_energy( rsd, rsd, pose, fa_mbsolv_score);
	emap[ FaMPSolv ] += fa_mbsolv_score;
}

/// @brief Specify Interaction Cutoff for computing pair energies
Distance
FaMPSolvEnergy::atomic_interaction_cutoff() const {
	etable::EtableCOP etable( etable_ );
	return etable->max_dis();
}

/// @brief Provide context graphs
void
FaMPSolvEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Energy Method for Scoring
void
FaMPSolvEnergy::setup_for_scoring(
	pose::Pose & pose, ScoreFunction const &
) const {

	// Initialize fullatom data
	init( pose );
}


/// @brief Finalize method after computing totals
void
FaMPSolvEnergy::finalize_total_energy(
	pose::Pose &,
	ScoreFunction const &,
	EnergyMap & // totals
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
	using namespace etable::count_pair;
	CountPairFunctionOP cpfxn( 0 );

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
			core::Real proj_i = fa_proj_[ rsd1.seqpos() ][ i ];
			core::Real proj_j = fa_proj_[ rsd2.seqpos() ][ j ];

			Vector const heavy_atom_j( rsd2.xyz( j ) );

			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			Real dummy_deriv( 0.0 );
			bool debug( false );

			score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), d2, dummy_deriv, proj_i, proj_j, debug );
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
	Real & deriv,
	Real const & f1,
	Real const & f2,
	bool &
) const {

	if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) return 0.0;

	// Initialize Variables
	Real temp_score( 0.0 );
	deriv = 0.0;
	bool const eval_deriv( true );

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

	temp_score = e1 + frac * ( e2 - e1 );

	// Always evaluate derivatives
	if ( eval_deriv ) {

		e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
		e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
		e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
		e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
		e1 = e11 + e21;
		e2 = e12 + e22;

		deriv = e1 + frac * ( e2 - e1 );
		deriv = deriv / std::sqrt( d2 );
	}
	return temp_score;
}

/// @brief Compute Change in Energy over distance (for minimization)
Real
FaMPSolvEnergy::eval_dE_dR_over_r(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	EnergyMap const &,
	Vector & F1,
	Vector & F2,
	Real const & f1,
	Real const & f2
) const {

	F1 = atom1.xyz().cross( atom2.xyz() );
	F2 = atom1.xyz() - atom2.xyz();
	Real d2 = atom1.xyz().distance_squared( atom2.xyz() );

	if ( ( d2 >= safe_max_dis2_ ) || ( d2 == Real(0.0) ) ) return 0.0;

	// bin by distance:
	Real const d2_bin = d2 * get_bins_per_A2_;
	int disbin = static_cast< int >( d2_bin ) + 1;
	Real frac = d2_bin - ( disbin - 1 );

	int const l1 = dsolv1_.index( disbin, atom1.type(), atom2.type()),
		l2 = l1 + 1;

	Real e11 = f1 * dsolv1_[ l1 ] + (1 - f1) * memb_dsolv1_[ l1 ];
	Real e12 = f1 * dsolv1_[ l2 ] + (1 - f1) * memb_dsolv1_[ l2 ];
	Real e21 = f2 * dsolv2_[ l1 ] + (1 - f2) * memb_dsolv2_[ l1 ];
	Real e22 = f2 * dsolv2_[ l2 ] + (1 - f2) * memb_dsolv2_[ l2 ];
	Real e1 = e11 + e21;
	Real e2 = e12 + e22;

	Real deriv = fa_weight_ * ( e1 + frac * ( e2 - e1 ) );

	return deriv / std::sqrt( d2 );
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

	core::Real thickness = pose.conformation().membrane_info()->membrane_thickness();
	core::Real steepness = pose.conformation().membrane_info()->membrane_steepness();

	// For convenience - grab nres
	Real nres = pose.total_residue();

	for ( Size i = 1; i <= nres; ++i ) {
		for ( Size j = 1, j_end = pose.residue( i ).nheavyatoms(); j <= j_end; ++j ) {

			Vector const xyz( pose.residue( i ).xyz( j ) );

			// Compute Standard Z Position
			fa_z_position_[i][j] = pose.conformation().membrane_info()->atom_z_position( pose.conformation(), i, j );

			// Compute Fa Projection
			fa_proj_[i][j] = compute_fa_proj( fa_z_position_[i][j], thickness, steepness );
		}
	}
}

/// @brief Helper Method - Compute Fa Proj
core::Real
FaMPSolvEnergy::compute_fa_proj(
	core::Real z_position,
	core::Real thickness,
	core::Real steepness
) const {

	Real internal_product(0), z(0), zn(0);
	internal_product = std::abs( z_position );
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	Real result = zn/(1 + zn);

	return result;
}

/// @brief Setup Data Members for Fullatom Info
void
FaMPSolvEnergy::setup_for_fullatom( pose::Pose & pose ) const {

	core::Real nres = pose.total_residue();

	fa_z_position_.resize( (core::Size)nres );
	fa_proj_.resize( (core::Size)nres );

	static Size const MAX_AMINOACID_SIZE = 15;

	for ( Size i = 1; i <= nres; ++i ) {

		Size const max_size = std::max( MAX_AMINOACID_SIZE, pose.residue( i ).nheavyatoms() );

		fa_z_position_[i].resize( max_size );
		fa_proj_[i].resize( max_size );

		for ( Size j = 1; j <= max_size; ++j ) {
			fa_z_position_[i][j] = 0.0;
			fa_proj_[i][j] = 0.0;
		}
	}
}

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaMPSolvEnergy_cc
