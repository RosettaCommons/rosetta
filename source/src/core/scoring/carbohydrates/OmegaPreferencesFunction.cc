// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/carbohydrates/OmegaPreferencesFunction.cc
/// @brief   Method definitions for OmegaPreferencesFunction.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit Headers
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.hh>
#include <core/scoring/carbohydrates/database_io.hh>

// Project Header
#include <core/types.hh>

// Basic Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <cmath>


static basic::Tracer TR( "core.scoring.carbohydrates.OmegaPreferencesFunction" );


namespace core {
namespace scoring {
namespace carbohydrates {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
/// @details  This class is only intended to be instantiated by the ScoringManager.
OmegaPreferencesFunction::OmegaPreferencesFunction()
{
	init();
}

OmegaPreferencesFunction::~OmegaPreferencesFunction()
{}


// Other Public Methods ///////////////////////////////////////////////////////
/// @details  TBA
/// @param    <preference>: AXIAL or GAUCHE_EFFECT
/// @param    <x>: an angle, in degrees, between 0 and 360
Energy
OmegaPreferencesFunction::operator()( OmegaPreferenceType preference, core::Angle x ) const {
	return evaluate_function( preference, x );
}

/// @details  TBD
/// @param    <preference>: AXIAL or GAUCHE_EFFECT
/// @param    <x>: an angle, in degrees, between 0 and 360
Real
OmegaPreferencesFunction::evaluate_derivative( OmegaPreferenceType preference, core::Angle x ) const
{
	if ( preference == PREFERENCE_NA ) { return 0.0; }
	set_parameters( preference, x );
	return 2 * k_ * (x - theta_);
}


// Private methods ////////////////////////////////////////////////////////////
void
OmegaPreferencesFunction::init()
{
	// TODO: Fill in, if we ever decide to make this function for complicated.
}

void
OmegaPreferencesFunction::set_parameters( OmegaPreferenceType preference, Angle x ) const
{
	// Note: I am only using private variables here for ease of conversion to another function in the future.  ~Labonte

	switch ( preference ) {
	case ANTI :
		if ( x <= 120 ) {
			theta_ = 60; b_ = 0.0;
		} else if ( x > 240 ) {
			theta_ = 300; b_ = 1.0;
		} else /* 120 < x <= 240 */ {
			theta_ = 180; b_ = 0.3;
		}
		break;
	case GAUCHE_EFFECT :
		if ( x <= 120 ) {
			theta_ = 60; b_ = 0.21;
		} else if ( x > 240 ) {
			theta_ = 300; b_ = 0.0;
		} else /* 120 < x <= 240 */ {
			theta_ = 180; b_ = 1.39;
		}
		break;
	case PREFERENCE_NA :
		;
	}
}

Energy
OmegaPreferencesFunction::evaluate_function( OmegaPreferenceType preference, Angle x ) const
{
	if ( preference == PREFERENCE_NA ) { return 0.0; }
	set_parameters( preference, x );
	return k_ * pow( ( x - theta_ ), 2 ) + b_;
}



// Sampling Methods ///////////////////////////////////////////////////////////

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
void
OmegaPreferencesFunction::setup_for_sampling( core::Real step_size ){
	using utility::to_string;

	//Note:
	// Probability from energy: -ln(p)=E -> p = e^-E

	//Get phi/get psi linkages.
	TR << "Setting up omega preference sampling" << std::endl;

	// Write 2 for loops.  Difference will be -180, 180; 360.
	for ( core::Size i = 1; i <= core::Size(N_OMEGA_PREFERENCE_TYPES); ++i ) {
		//TR << "linkage: " << i << std::endl;
		//std::cout << "Angle,Energy,Probability"<<std::endl;

		OmegaPreferenceType preference_type = static_cast< OmegaPreferenceType >( i );

		OmegaPreferenceSamplingData sampling_data;

		sampling_data.preference_type = preference_type;
		sampling_data.step_size = step_size;

		for ( Angle dih = -180.0; dih <= 180.0; dih+=step_size ) {
			Energy e = evaluate_function(preference_type, dih);
			Probability prob = std::exp( -e );

			sampling_data.angles.push_back( dih );
			sampling_data.probabilities.push_back( prob );
			//std::cout <<dih << "," << e << "," << prob << std::endl;

		}

		dihedral_sampling_data_[ preference_type ] = sampling_data;
	}
}

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
OmegaPreferencesFunction::sampling_data_setup() const {
	if ( dihedral_sampling_data_.empty() ) {
		return false;
	} else {
		return true;
	}
}

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
OmegaPreferencesFunction::sampling_data_setup( OmegaPreferenceType const linkage_type ) const {
	std::map< OmegaPreferenceType, OmegaPreferenceSamplingData>::const_iterator const_iter;
	const_iter = dihedral_sampling_data_.find( linkage_type );

	if ( const_iter != dihedral_sampling_data_.end() ) {
		return true;
	} else {
		return false;
	}
}

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
OmegaPreferenceSamplingData const &
OmegaPreferencesFunction::get_sampling_data( OmegaPreferenceType const linkage_type ) const {

	std::map< OmegaPreferenceType, OmegaPreferenceSamplingData>::const_iterator const_iter;
	const_iter = dihedral_sampling_data_.find( linkage_type );


	if ( const_iter != dihedral_sampling_data_.end() ) {
		return const_iter->second;
	} else {
		std::string m = "Omega Preference Data not found for linkage type "+utility::to_string( linkage_type );
		utility_exit_with_message( m );
	}
}




}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core
