// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPEnvEnergy.cc
///
/// @brief  Membrane Environemnt Energy
/// @details One Body Term - score residue interaction with specific hydrophobic layer
///    derived from Membrane base potential and uses mpframework data
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPEnvEnergy_cc
#define INCLUDED_core_scoring_membrane_MPEnvEnergy_cc

// Unit Headers
#include <core/scoring/membrane/MPEnvEnergy.hh>
#include <core/scoring/membrane/MPEnvEnergyCreator.hh>

// Project Headers
#include <core/scoring/membrane/MembraneData.hh>
#include <core/scoring/ScoringManager.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/scoring/EnergyMap.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

using namespace core::scoring;
using namespace core::scoring::methods;

static thread_local basic::Tracer TR( "core.membrane.MPEnvEnergy" );

namespace core {
namespace scoring {
namespace membrane {

/// Creator Methods ////////////////////

/// @brief Return a Fresh Instance of the Energy Method
methods::EnergyMethodOP
MPEnvEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MPEnvEnergy );
}

/// @brief Return Appropriate Score Type Name
ScoreTypes
MPEnvEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( MPEnv );
	return sts;
}

/// @brief Default Constructor
MPEnvEnergy::MPEnvEnergy() :
	parent( EnergyMethodCreatorOP( new MPEnvEnergyCreator ) ),
	mpdata_( ScoringManager::get_instance()->get_MembraneData() )
	// penalties_( true )
{}

/// @brief Clone an Energy Method
EnergyMethodOP
MPEnvEnergy::clone() const
{
	return EnergyMethodOP( new MPEnvEnergy );
}

/// Scoring Methods //////////////////////

/// @brief Setup Method for Scoring
void
MPEnvEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	mpdata_.compute_centroid_environment( pose );
}

/// @brief Setup Energy Method for Derivatives
void
MPEnvEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const
{
	setup_for_scoring( pose, sf );
}

/// @brief Sore Residue-Environemnt Energy Based on the Membrane Layer and Residue Z_Position
void
MPEnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{

	// Check Structure is a membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message("Error: Cannot use mpframework energy term using a non membrane pose!");
	}

	// Skip Virtual Residues
	if ( rsd.aa() == core::chemical::aa_vrt ) return;

	// If residue is not a protein residue, set to 0 and return
	if ( !rsd.is_protein() ) {
		emap[ MPEnv ] = 0.0;
		return;
	}

	// Initialize Env Score
	Real env_score( 0.0 );

	// Initialize Info for Scoring
	CenListInfo const & cenlist = mpdata_.get_cenlist_from_pose( pose );
	core::Real const z_position = pose.conformation().membrane_info()->residue_z_position( rsd.seqpos() );
	core::Size const seqpos = rsd.seqpos();
	core::chemical::AA const & aa = rsd.aa();

	// Compute MPEnv Score
	env_score = compute_mpenv_score( cenlist, aa, z_position, seqpos );

	// Set the new score in the energies map
	emap[ MPEnv ] += env_score;

} // residue_energy

/// @brief Finalize Total Energy in the Poose
void
MPEnvEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const
{
	mpdata_.finalize( pose );
}

/// @brief Versioning
core::Size
MPEnvEnergy::version() const
{
	return 2;
}

/// @brief Evaluate Environemnt Energy for Residue based on z_position, aa, and layer
core::Real
MPEnvEnergy::compute_mpenv_score(
	CenListInfo const & cenlist,
	core::chemical::AA const & aa,
	core::Real const z_position,
	core::Size const seqpos
) const
{
	// Initialise score
	core::Real menv_score( 0.0 );

	// Setup Initial Conditions
	Real t2 = 2.0;
	Real t3 = 2.0;

	int s2 = 14;
	int s3 = 14;
	int layer1, layer2, layer;
	Real f, z, zn, low;

	// Setup Weights
	Real const env6_weight = 1.0;
	Real const env10_weight = 1.0;

	// Extract cenlist info
	Real fcen6  ( cenlist.fcen6( seqpos ) );
	Real fcen10 ( cenlist.fcen10( seqpos ) );

	// in rare cases, the density is over 15 within 6Å
	if ( fcen6 > 15 ) fcen6 = 15 ;
	if ( fcen10 > 40 ) fcen10 = 40;

	// Score in the Pure Water Layer
	if ( ( z_position < -19.0 ) || ( z_position > 19.0 ) ) {

		layer = 3;
		Real score6( env6_weight * mpdata_.mem_env_log6()( aa, layer, static_cast< int >( fcen6 ) ) );
		Real score10( env10_weight * mpdata_.mem_env_log10()( aa, layer, static_cast< int >( fcen10 ) ) );
		menv_score = score6 + score10;

		// Interpolate Between Water and Interface Phases
	} else if ( ( z_position >= -19.0 && z_position <= -17.0 ) || ( z_position >= 17.0 && z_position <= 19.0 ) ) {

		layer1 = 2; //interface layer
		layer2 = 3; //water layer

		if ( z_position <= -17.0 ) {
			low = -17.0;
		} else {
			low = 17.0;
		}

		z = 2*std::abs( (z_position - low) ) / t2;
		zn = std::pow( z, s2 );
		f = zn/(1 + zn);

		// Score By Layer
		Real score6_layer2( env6_weight * mpdata_.mem_env_log6()( aa, layer2, static_cast< int >( fcen6 ) ) );
		Real score10_layer2( env10_weight * mpdata_.mem_env_log10()( aa, layer2, static_cast< int >( fcen10 ) ) );
		Real score6_layer1( env6_weight * mpdata_.mem_env_log6()( aa, layer1, static_cast< int >( fcen6 ) ) );
		Real score10_layer1( env10_weight * mpdata_.mem_env_log10()( aa, layer1, static_cast< int >( fcen10 ) ) );

		// Increment Score
		menv_score = f * ( score6_layer2 + score10_layer2 ) + ( 1 - f ) * ( score6_layer1 + score10_layer1 );

		if ( z_position <= -18.0 || z_position >= 18.0 ) {
			layer = 2;
		} else {
			layer = 3;
		}

		// Pure Interface Phase
	} else if ( ( z_position > -17.0 && z_position < -13.0 ) || ( z_position > 13.0 && z_position < 17.0 ) ) {

		layer = 2; //interface layer

		Real score6 ( env6_weight * mpdata_.mem_env_log6()( aa, layer, static_cast< int >( fcen6 ) ) );
		Real score10 ( env10_weight * mpdata_.mem_env_log10()( aa, layer, static_cast< int >( fcen10 ) ) );

		menv_score = score6 + score10;

	} else if ( ( z_position >= -13.0 && z_position <= -11.0 ) || ( z_position >= 11.0 && z_position <= 13.0 ) ) {

		// interpolate between interface and hydrophobic phases

		layer1 = 1; //hydrophobic layer
		layer2 = 2; //interface layer

		if ( z_position <= -11.0 ) low = -11.0;
		else low = 11.0;

		z = 2*std::abs( ( z_position - low) ) / t3;
		zn = std::pow( z, s3 );
		f = zn/(1 + zn);

		// Layer Based Statistics
		Real score6_layer2( env6_weight * mpdata_.mem_env_log6()( aa, layer2, static_cast< int >( fcen6 ) ) );
		Real score10_layer2( env10_weight * mpdata_.mem_env_log10()( aa, layer2, static_cast< int >( fcen10 ) ) );
		Real score6_layer1( env6_weight * mpdata_.mem_env_log6()( aa, layer1, static_cast< int >( fcen6 ) ) );
		Real score10_layer1( env10_weight * mpdata_.mem_env_log10()( aa, layer1, static_cast< int >( fcen10 ) ) );

		// Update Score
		menv_score = f * ( score6_layer2  + score10_layer2 ) + ( 1 - f ) * ( score6_layer1 + score10_layer1 );

		// recompute layer
		if ( z_position <= -12.0 || z_position >= 12.0 ) layer = 2;
		else layer = 1;

		// Pure Hydrophobic Phase
	} else {

		layer = 1;

		Real score6( env6_weight * mpdata_.mem_env_log6()(  aa, layer, static_cast< int >( fcen6  ) ) );
		Real score10( env10_weight * mpdata_.mem_env_log10()( aa, layer, static_cast< int >( fcen10 ) ) );

		menv_score = score6 + score10;
	}

	menv_score *= 0.5; // pretty sure this weight should NOT be hard coded @ralford 3/29/14

	// Return score
	return menv_score;
}

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPEnvEnergy_cc
