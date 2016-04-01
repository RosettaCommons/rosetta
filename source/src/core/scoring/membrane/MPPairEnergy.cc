// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPPairEnergy.cc
///
/// @brief  Membrane Residue Pair Energy Term
/// @details Two Body Term - score residue-residue interactions in the membrane. Derived from Membrane
///    base potential and uses mpframework data
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPPairEnergy_cc
#define INCLUDED_core_scoring_membrane_MPPairEnergy_cc

// Unit Headers
#include <core/scoring/membrane/MPPairEnergy.hh>
#include <core/scoring/membrane/MPPairEnergyCreator.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/membrane/SpanningTopology.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

// Symmetry Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace membrane {

/// Creator methods ////////////////////

/// @brief Create a Fresh Instance of the MP Pair Energy Method
methods::EnergyMethodOP
MPPairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPPairEnergy );
}

/// @brief Return the relevant score type - MPPair
ScoreTypes
MPPairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPPair );
	return sts;
}

/// @brief Default Constructor
MPPairEnergy::MPPairEnergy() :
	parent( EnergyMethodCreatorOP( new MPPairEnergyCreator ) ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() ),
	no_interpolate_mpair_( true ) // temp
{}

/// @brief Clone Method
EnergyMethodOP
MPPairEnergy::clone() const {
	return EnergyMethodOP( new MPPairEnergy() );
}

/// Scoring Methods ///////////////////////

/// @brief Setup by updating residue neighbors and computing cen env
void
MPPairEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	pose.update_residue_neighbors();
	mpdata_.compute_centroid_environment( pose );
}

/// @brief Compute Reisdue Pair Energy in the membrane
void
MPPairEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	// Check Structure is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}

	// Setup Initial Score
	Real pair_score( 0.0 );

	// Setup Pair Positions
	Size seqpos1( rsd1.seqpos() );
	Size seqpos2( rsd2.seqpos() );

	// Pull amino acids
	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	// Grab Initial z_positions
	conformation::Conformation const & conf( pose.conformation() );
	Real const z_position1( conf.membrane_info()->residue_z_position( conf, rsd1.seqpos() ) );
	Real const z_position2( conf.membrane_info()->residue_z_position( conf, rsd2.seqpos() ) );

	// Compute Centroid Info
	conformation::Atom const & cen1 ( rsd1.atom( rsd1.nbr_atom() ) ), cen2 (rsd2.atom( rsd2.nbr_atom() ) );
	Real const cendist = cen1.xyz().distance_squared( cen2.xyz() );

	// Optimizations for Symmetry
	if ( core::pose::symmetry::is_symmetric( pose ) ) {

		using namespace core::conformation::symmetry;

		SymmetricConformation const & symm_conf (
			dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );


		// Ignore Virtual Symmetry Residues
		if ( symm_info->is_virtual( seqpos1 ) ) seqpos1 = 0;
		if ( symm_info->is_virtual( seqpos2 ) ) seqpos2 = 0;

		// Correct for bb_is_independent
		if ( !symm_info->bb_is_independent( seqpos1 ) ) {
			seqpos1 = symm_info->bb_follows( seqpos1 );
		}

		if ( !symm_info->bb_is_independent( seqpos2 ) ) {
			seqpos2 = symm_info->bb_follows( seqpos2 );
		}
	}

	// If pair is at a termini, return
	if ( seqpos1 == 0 || seqpos2 == 0 ) return;

	// Only score protein residues
	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	// Skip Pair Scoring if disulfide
	if ( aa1 == chemical::aa_cys && aa2 == chemical::aa_cys &&
			rsd1.is_bonded( rsd2 ) && rsd1.polymeric_sequence_distance( rsd2 ) > 1 &&
			rsd1.has_variant_type( chemical::DISULFIDE ) && rsd2.has_variant_type( chemical::DISULFIDE ) ) return;

	// no pair score for residues closer than 9 in sequence
	if ( rsd1.polymeric_sequence_distance( rsd2 ) <= 8 ) return;

	// Compute Pair Score
	pair_score = compute_mpair_score( z_position1, z_position2, aa1, aa2, cendist );

	// Set Score in the Energy Map
	emap[ MPPair ] += pair_score;
}

void
MPPairEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const {
	mpdata_.finalize( pose );
}

Distance
MPPairEnergy::atomic_interaction_cutoff() const {
	return 6.0; /// now subtracted off 6.0 from cutoffs in MembraneCentroid params files
}

core::Size
MPPairEnergy::version() const { return 2; }

/// @brief Compute Residue-Residue Energy from z_positions, aas and distance between Cbs (cendist)
core::Real
MPPairEnergy::compute_mpair_score(
	core::Real const z_position1,
	core::Real const z_position2,
	chemical::AA const & aa1,
	chemical::AA const & aa2,
	core::Real const cendist
) const {

	// Initialize Pair Score
	core::Real mpair_score = 0.0;

	// Initial Conditions
	int icon = 5;
	Real interp2( 0.0 );

	// Interpolate between layers
	int hydro_layer = 1;  // 1 = not_hydrophobic_core, 2 = hydrophobic core
	Real avg_z_position = ( z_position1 + z_position2 )/2;

	// Both residues must be in the hydrophobic core for this type of scoring
	if ( z_position1 > -12 && z_position1 < 12 && z_position2 > -12 && z_position2 < 12 ) hydro_layer = 2;

	// Score based on centroid-centroid distance between neighbors
	if ( cendist > mpdata_.cen_dist10_pad_plus() ) {
		icon = 4;
		interp2 = ( cendist + mpdata_.cen_dist12_pad_minus() ) * mpdata_.cen_dist12_pad_hinv();
	} else {
		if ( cendist > mpdata_.cen_dist7_pad_plus() ) {
			icon = 3;
			interp2 = ( cendist + mpdata_.cen_dist10_pad_minus() ) * mpdata_.cen_dist10_pad_hinv();
		} else {
			if ( cendist > mpdata_.cen_dist5_pad_plus() ) {
				icon = 2;
				interp2 = ( cendist + mpdata_.cen_dist7_pad_minus() ) * mpdata_.cen_dist7_pad_hinv();
			} else {
				icon = 1;
				interp2 = ( cendist + mpdata_.cen_dist5_pad_minus() ) * mpdata_.cen_dist5_pad_hinv();
			}
		}
	}

	if ( interp2 < 0.0 ) interp2 = 0.0;

	// note in theory this will never happen but in practice round off
	// error can cause problem
	if ( interp2 > 1.0 ) interp2 = 1.0;

	// handle last bin specially since icon+1 would be past array end
	Real f(0);

	if ( icon == 5 ) {
		mpair_score = ( 1.0f - interp2 ) * mpdata_.mem_pair_log()( hydro_layer, icon, aa1, aa2 );
	} else {
		// if interpolation is true (by default):
		if ( !no_interpolate_mpair_ ) {

			// bin 1
			if ( std::abs( avg_z_position - -12 ) < 4 ) {

				f = 1/(1+std::exp(1.5*(18-avg_z_position)));
				mpair_score = ( ( 1.0f - interp2 ) * ((1-f)*mpdata_.mem_pair_log()( 1, icon  , aa1, aa2 ) + f*mpdata_.mem_pair_log()( 2, icon  , aa1, aa2 )) +
					(        interp2 ) * ((1-f)*mpdata_.mem_pair_log()( 1, icon+1, aa1, aa2 ) + f*mpdata_.mem_pair_log()( 2, icon+1, aa1, aa2 )));
				// bin 2
			} else if ( std::abs( avg_z_position - 12 ) < 4 ) {

				f=1/(1+std::exp(1.5*(avg_z_position-42)));
				mpair_score = ( ( 1.0f - interp2 ) * ((1-f)*mpdata_.mem_pair_log()( 1, icon  , aa1, aa2 ) + f*mpdata_.mem_pair_log()( 2, icon  , aa1, aa2 )) +
					(        interp2 ) * ((1-f)*mpdata_.mem_pair_log()( 1, icon+1, aa1, aa2 ) + f*mpdata_.mem_pair_log()( 2, icon+1, aa1, aa2 )));

				// bin 3
			}  else {
				mpair_score = ( ( 1.0f - interp2 ) * mpdata_.mem_pair_log()( hydro_layer, icon  , aa1, aa2 ) +
					(        interp2 ) * mpdata_.mem_pair_log()( hydro_layer, icon+1, aa1, aa2 ));
			}

			// Otherwise, don't interpolate
		} else {
			mpair_score = ( ( 1.0f - interp2 ) * mpdata_.mem_pair_log()( hydro_layer, icon  , aa1, aa2 ) +
				(        interp2 ) * mpdata_.mem_pair_log()( hydro_layer, icon+1, aa1, aa2 ));
		}

	}

	mpair_score *= 2.019; // if this is a weight, it should be applied somewhere else

	// Return result
	return mpair_score;
}

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPPairEnergy_cc

