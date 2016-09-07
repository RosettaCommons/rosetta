// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#include <core/scoring/SecondaryStructureWeights.hh>

// numeric headers
#include <numeric/random/random.hh>

#include <utility/vector1.hh>

// utility headers
// ObjexxFCL headers
// C++ headers

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {


SecondaryStructureWeights::~SecondaryStructureWeights() = default;

/// @brief initialize to default values, also load proper data bins from external files
void
SecondaryStructureWeights::initialize()
{
	// weights
	setup_parallel_antiparallel_weights();

	// additional settings
	set_ss_lowstrand( 0 );
	set_ss_cutoff( 6 );
	// set_ss_lowstrand( 1 );
	// set_ss_cutoff( 11 );

	set_localstrandpair_penalty( 0.0 );
	set_seq_sep_scale( 20.0 );
	set_max_strand_dist_cutoff( 12.0 );
	set_strand_dist_cutoff( 6.5 );
	set_stretch_strand_dist_cutoff( false );
	set_handedness_score_flag( false );
}


/// @details if not randomizing weights, all weights set to 1.0
void
SecondaryStructureWeights::setup_parallel_antiparallel_weights(
	bool const & randomize_weights
) {
	if ( randomize_weights ) {
		// Choose whether to weight up parallel or antiparallel
		Real randomnumber = numeric::random::rg().uniform();
		Real randomweight = 10.0 * numeric::random::rg().uniform(); // Pretty drastic reweighting...
		if ( randomnumber < 0.5 ) {
			set_parallel_weight( randomweight );
			set_antiparallel_weight( 1.0/randomweight );
		} else {
			set_parallel_weight( 1.0/randomweight );
			set_antiparallel_weight( randomweight );
		}
	} else { // default weights
		set_parallel_weight( 1.0 );
		set_antiparallel_weight( 1.0 );
	}
}


}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::SecondaryStructureWeights::save( Archive & arc ) const {
	arc( CEREAL_NVP( parallel_weight_ ) ); // Real
	arc( CEREAL_NVP( antiparallel_weight_ ) ); // Real
	arc( CEREAL_NVP( ss_lowstrand_ ) ); // int
	arc( CEREAL_NVP( ss_cutoff_ ) ); // int
	arc( CEREAL_NVP( localstrandpair_penalty_ ) ); // Real
	arc( CEREAL_NVP( seq_sep_scale_ ) ); // Real
	arc( CEREAL_NVP( max_strand_dist_cutoff_ ) ); // Real
	arc( CEREAL_NVP( strand_dist_cutoff_ ) ); // Real
	arc( CEREAL_NVP( stretch_strand_dist_cutoff_ ) ); // _Bool
	arc( CEREAL_NVP( handedness_score_flag_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::SecondaryStructureWeights::load( Archive & arc ) {
	arc( parallel_weight_ ); // Real
	arc( antiparallel_weight_ ); // Real
	arc( ss_lowstrand_ ); // int
	arc( ss_cutoff_ ); // int
	arc( localstrandpair_penalty_ ); // Real
	arc( seq_sep_scale_ ); // Real
	arc( max_strand_dist_cutoff_ ); // Real
	arc( strand_dist_cutoff_ ); // Real
	arc( stretch_strand_dist_cutoff_ ); // _Bool
	arc( handedness_score_flag_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::SecondaryStructureWeights );
CEREAL_REGISTER_TYPE( core::scoring::SecondaryStructureWeights )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_SecondaryStructureWeights )
#endif // SERIALIZATION
