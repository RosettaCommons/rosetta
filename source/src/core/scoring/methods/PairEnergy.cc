// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/PairEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/PairEnergy.hh>
#include <core/scoring/methods/PairEnergyCreator.hh>

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/PairEPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the PairEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
PairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new PairEnergy );
}

ScoreTypes
PairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_pair );
	sts.push_back( fa_pair_aro_aro );
	sts.push_back( fa_pair_aro_pol );
	sts.push_back( fa_pair_pol_pol );
	return sts;
}


PairEnergy::PairEnergy() :
	parent( methods::EnergyMethodCreatorOP( new PairEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_PairEPotential() )
{}


/// clone
EnergyMethodOP
PairEnergy::clone() const
{
	return EnergyMethodOP( new PairEnergy() );
}


void
PairEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
}


void
PairEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


void
PairEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

void
PairEnergy::prepare_rotamers_for_packing(
	pose::Pose const & /*pose*/,
	conformation::RotamerSetBase & set
) const
{
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		set.nonconst_rotamer( ii )->update_actcoord(); // the Rotamer set does not take responsibility for this; why not? -- maybe Residue could?
	}
}

void
PairEnergy::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	pose.update_actcoord( resid );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
PairEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ||
			(!rsd1.is_polar() && !rsd1.is_aromatic() ) ||
			(!rsd2.is_polar() && !rsd2.is_aromatic() ) ) return;

	/// Enforce interaction cutoff
	Real const rsd1_reach( rsd1.nbr_radius() ), rsd2_reach( rsd2.nbr_radius() );
	Distance const intxn_dist = rsd1_reach + rsd2_reach + interaction_cutoff();
	DistanceSquared const intxn_dist2 = intxn_dist * intxn_dist;
	DistanceSquared const nbr_dist2 = rsd1.xyz( rsd1.nbr_atom() ).distance_squared( rsd2.xyz( rsd2.nbr_atom() ) );
	if ( nbr_dist2 > intxn_dist2 ) return;

	TenANeighborGraph const & tenA_neighbor_graph
		( pose.energies().tenA_neighbor_graph() );

	Real pairE = potential_.pair_term_energy(
		rsd1,
		tenA_neighbor_graph.get_node( rsd1.seqpos() )->
		num_neighbors_counting_self_static(),
		rsd2,
		tenA_neighbor_graph.get_node( rsd2.seqpos() )->
		num_neighbors_counting_self_static() );

	if ( rsd1.is_polar() && rsd2.is_polar() ) {
		emap[ fa_pair_pol_pol ] += pairE;
		emap[ fa_pair ] += pairE;
	} else if ( rsd1.is_aromatic() && rsd2.is_aromatic() ) {
		emap[ fa_pair_aro_aro ] += pairE;
	} else {
		emap[ fa_pair_aro_pol ] += pairE;
	}

	//if ( rsd1.actcoord_atoms().size() == 1 && rsd1.actcoord().distance( rsd1.xyz( rsd1.actcoord_atoms()[ 1 ] )) > 0.0001 ) {
	// std::cout << "Actcoord discrepancy!" << std::endl;
	//}

	//if ( rsd2.actcoord_atoms().size() == 1 && rsd2.actcoord().distance( rsd2.xyz( rsd2.actcoord_atoms()[ 1 ] )) > 0.0001 ) {
	// std::cout << "Actcoord discrepancy2!" << std::endl;
	//}


	//if ( pairE != 0.0 )
	// std::cout << "  pairE " << rsd1.seqpos() << " " << rsd2.seqpos() << " " << pairE << " " << std::sqrt( nbr_dist2 ) << std::endl;
	//std::cout << "  pairE: " << rsd1.seqpos() << " " << rsd2.seqpos() << " " << pairE << std::endl;
}

void
PairEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	using namespace conformation;
	using namespace numeric;

	EnergyMap emap;

	bool const any_aro( weights[ fa_pair_aro_pol ] != 0 || weights[ fa_pair_aro_aro ] != 0 );

	for ( Size ii = 1; ii <= set1.get_n_residue_types(); ++ii ) {
		if ( set1.get_n_rotamers_for_residue_type( ii ) == 0 ) continue;
		Size const ii_offset = set1.get_residue_type_begin( ii );
		Residue const & ii_example_rotamer( set1.rotamer_ref( ii_offset ));
		if ( ! ii_example_rotamer.is_protein() || ! ( ii_example_rotamer.is_polar() || (any_aro && ii_example_rotamer.is_aromatic()) ) ) continue;


		Vector const & ii_coord( ii_example_rotamer.atom( ii_example_rotamer.type().nbr_atom() ).xyz());
		Real const ii_radius( ii_example_rotamer.type().nbr_radius() );

		for ( Size jj = 1; jj <= set2.get_n_residue_types(); ++jj ) {
			if ( set2.get_n_rotamers_for_residue_type( jj ) == 0 ) continue;
			Size const jj_offset = set2.get_residue_type_begin( jj );
			Residue const & jj_example_rotamer( set2.rotamer_ref( jj_offset ));

			if ( ! jj_example_rotamer.is_protein() || ! ( jj_example_rotamer.is_polar() || (any_aro &&  jj_example_rotamer.is_aromatic()) ) ) continue;

			Vector const & jj_coord( jj_example_rotamer.atom( jj_example_rotamer.type().nbr_atom() ).xyz());
			Real const jj_radius( jj_example_rotamer.type().nbr_radius() );

			if ( ii_coord.distance_squared( jj_coord ) >= std::pow(ii_radius+jj_radius+atomic_interaction_cutoff(), 2 ) ) continue;

			for ( Size kk = 1, kke = set1.get_n_rotamers_for_residue_type( ii ); kk <= kke; ++kk ) {
				Size const kk_rot_id = ii_offset + kk - 1;
				for ( Size ll = 1, lle = set2.get_n_rotamers_for_residue_type( jj ); ll <= lle; ++ll ) {
					Size const ll_rot_id = jj_offset + ll - 1;

					//emap.zero();
					emap[ fa_pair ] = 0;
					emap[ fa_pair_aro_aro ] = 0;
					emap[ fa_pair_aro_pol ] = 0;
					emap[ fa_pair_pol_pol ] = 0;
					PairEnergy::residue_pair_energy( set1.rotamer_ref( kk_rot_id ), set2.rotamer_ref( ll_rot_id ), pose, sfxn, emap );
					energy_table( ll_rot_id, kk_rot_id ) += static_cast< core::PackerEnergy > (
						weights[ fa_pair ] * emap[ fa_pair ] +
						weights[ fa_pair_aro_aro ] * emap[ fa_pair_aro_aro ] +
						weights[ fa_pair_aro_pol ] * emap[ fa_pair_aro_pol ] +
						weights[ fa_pair_pol_pol ] * emap[ fa_pair_pol_pol ]
					);
				}
			}
		}
	}

}

void
PairEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	if ( ! residue.is_protein() || ! ( residue.is_polar() || residue.is_aromatic() ) ) return;
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		//emap.zero();
		emap[ fa_pair ] = 0;
		emap[ fa_pair_aro_aro ] = 0;
		emap[ fa_pair_aro_pol ] = 0;
		emap[ fa_pair_pol_pol ] = 0;
		if ( ! set.rotamer(ii)->is_protein() || ! ( set.rotamer(ii)->is_polar() || set.rotamer(ii)->is_aromatic() ) ) continue;
		PairEnergy::residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
		energy_vector[ ii ] += static_cast< core::PackerEnergy > (weights.dot( emap ) );
	}
}


void
PairEnergy::evaluate_rotamer_background_energy_maps(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & , // weights
	utility::vector1< EnergyMap > & emaps
) const
{
	if ( ! residue.is_protein() || ! residue.is_polar() || ! residue.is_aromatic()  ) return;
	EnergyMap emap;
	for ( Size ii = 1, ii_end = set.num_rotamers(); ii <= ii_end; ++ii ) {
		//emap.zero();
		emap[ fa_pair ] = 0;
		emap[ fa_pair_aro_aro ] = 0;
		emap[ fa_pair_aro_pol ] = 0;
		emap[ fa_pair_pol_pol ] = 0;
		if ( ! set.rotamer(ii)->is_protein() || ! ( set.rotamer(ii)->is_polar() || set.rotamer(ii)->is_aromatic() ) ) continue;
		PairEnergy::residue_pair_energy( *set.rotamer( ii ), residue, pose, sfxn, emap );
		emaps[ ii ] += emap;
	}
}


/// @details Returns false if !res_movign_wrt_eachother since the score function does not update
/// the neighbor counts for residues during minimization.  If two residues are not moving wrt
/// each other, their scores are not changing during minimization, even though this is a context
/// dependent score function.
bool
PairEnergy::defines_score_for_residue_pair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const
{
	if ( ! res_moving_wrt_eachother ) return false;
	if ( ! res1.is_protein() || ( ! res1.is_polar() && ! res1.is_aromatic())  ) return false;
	if ( ! res2.is_protein() || ( ! res2.is_polar() && ! res2.is_aromatic())  ) return false;
	return true;
}

void
PairEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	debug_assert( (rsd1.is_polar() || rsd1.is_aromatic()) && (rsd2.is_polar() || rsd2.is_aromatic() ) );

	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	int const nbr_count1( tenA_neighbor_graph.get_node( rsd1.seqpos() )->
		num_neighbors_counting_self_static() );
	int const nbr_count2( tenA_neighbor_graph.get_node( rsd2.seqpos() )->
		num_neighbors_counting_self_static() );

	Real weight(0.0); Size naros(0);
	if ( rsd1.is_aromatic() ) {
		++naros;
	}
	if ( rsd2.is_aromatic() ) {
		++naros;
	}
	switch (naros ) {
	case 0 :
		weight = std::max( weights[ fa_pair ], weights[ fa_pair_pol_pol ] );
		break;
	case 1 :
		weight = weights[ fa_pair_aro_pol ];
		break;
	case 2 :
		weight = weights[ fa_pair_aro_aro ];
		break;
	default :
		utility_exit_with_message( "ERROR in fa_pair derivaties, too many aromatics!!!");
		break;
	}

	if ( weight == 0.0 ) return;

	Real dpairE_dr( 0.0 );
	potential_.pair_term_energy_and_deriv( rsd1, nbr_count1, rsd2, nbr_count2, dpairE_dr );
	Vector f1( cross( rsd1.actcoord(), rsd2.actcoord() ) );
	Vector f2( rsd1.actcoord() - rsd2.actcoord() );
	Real const dis( f2.length() );

	if ( dis == Real(0.0 ) ) {
		utility_exit_with_message("dis==0 in pairtermderiv!");
	}
	dpairE_dr /= dis;
	f1 *= weight * dpairE_dr;
	f2 *= weight * dpairE_dr;

	/// APL: so where does this next section come from?  Well, previously the code assumed that
	/// all actcoord atoms were controlled by the same chi (or in arg's case, were moved synchronously
	/// by chi4 so that their motions canceled out, assuming the CN bond lengths were ideal)
	/// but this assumption does not hold in the case that bond angles or bond lengths are allowed
	/// to move.  The code below gives accurate derivatives by splitting up the derivative
	/// contribution to the atoms that are controlling the actcoord.
	Vector f1_1( f1 ), f2_1( f2 ), f1_2( f1 ), f2_2( f2 );
	Size const nact1 = rsd1.actcoord_atoms().size();
	Size const nact2 = rsd2.actcoord_atoms().size();
	Real inv_nact1( 1 / (Real) nact1 ), inv_nact2( -1 / (Real) nact2 );
	f1_1 *= inv_nact1;
	f2_1 *= inv_nact1;
	f1_2 *= inv_nact2;
	f2_2 *= inv_nact2;

	for ( Size ii = 1; ii <= nact1; ++ii ) {
		r1_atom_derivs[ rsd1.actcoord_atoms()[ ii ]].f1() += f1_1;
		r1_atom_derivs[ rsd1.actcoord_atoms()[ ii ]].f2() += f2_1;
	}
	for ( Size ii = 1; ii <= nact2; ++ii ) {
		r2_atom_derivs[ rsd2.actcoord_atoms()[ ii ]].f1() += f1_2;
		r2_atom_derivs[ rsd2.actcoord_atoms()[ ii ]].f2() += f2_2;
	}
}

/// @brief PairEnergy distance cutoff set to the same cutoff used by EtableEnergy, for now
Distance
PairEnergy::atomic_interaction_cutoff() const
{
	return interaction_cutoff();
}

/// @details non-virtual accessor for speed; assumption: PairEnergy is not inherrited from.
Distance
PairEnergy::interaction_cutoff() const
{
	//return 5.5; // TOO SHORT.  GOES OUT TO 7.5 A if nbins == 2; if nbins == 3, then it goes out to 9 A.
	return potential_.range();
}

/// @brief PairEnergy requires that Energies class maintains a TenANeighborGraph
void
PairEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ ten_A_neighbor_graph ] = true;
}

/// @brief PairEnergy does not define intraresidue interactions
bool
PairEnergy::defines_intrares_energy( EnergyMap const & /*weights*/ ) const
{
	return false;
}

void
PairEnergy::eval_intrares_energy(
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const {}
core::Size
PairEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
