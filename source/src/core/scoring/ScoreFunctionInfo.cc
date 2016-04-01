// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunctionInfo.cc
/// @brief  Score function class descriptor
/// @author Phil Bradley

// Unit Headers
#include <core/scoring/ScoreFunctionInfo.hh>

// Package Headers
#include <core/scoring/ContextGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/ContextGraphFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

/// default constructor -- fill this in
ScoreFunctionInfo::ScoreFunctionInfo() : utility::pointer::ReferenceCount(),
	max_atomic_interaction_distance_( 0.0 ),
	max_context_neighbor_cutoff_( 0.0 ),
	context_graphs_required_( num_context_graph_types, false ),
	energy_method_options_( /* 0 */ )
{}

/// @details Make sure to make a copy of the EMOpts, instead of just copying the pointer;
/// we don't want SFIs to share pointers to the same objects or the object could change
/// underneath the SFI.
ScoreFunctionInfo::ScoreFunctionInfo( ScoreFunctionInfo const & src ) : utility::pointer::ReferenceCount(),
	max_atomic_interaction_distance_( src.max_atomic_interaction_distance_ ),
	max_context_neighbor_cutoff_( src.max_context_neighbor_cutoff_ ),
	context_graphs_required_( src.context_graphs_required_ ),
	scores_present_( src.scores_present_ ),
	energy_method_options_( src.energy_method_options_ ? new methods::EnergyMethodOptions( * src.energy_method_options_ ) : 0 )
{}

ScoreFunctionInfo::~ScoreFunctionInfo() {}

ScoreFunctionInfo::ScoreFunctionInfo( ScoreFunction const & scorefxn )
:
	context_graphs_required_( num_context_graph_types, false )
{
	initialize_from( scorefxn );
}

/// @brief initializes three peices of data that describe the score function,
/// the atomic interaction distance, the context neighbor distance, and the
/// context graphs required by the scoring function to be properly evaluated
///
/// now also including the EnergyMethodOptions object, which holds eg the etable
/// name and (currently) the reference energy aa-weights
///
void
ScoreFunctionInfo::initialize_from( ScoreFunction const & scorefxn ) {

	max_atomic_interaction_distance_ = scorefxn.max_atomic_interaction_cutoff();

	std::fill( context_graphs_required_.begin(), context_graphs_required_.end(), false );
	scorefxn.indicate_required_context_graphs( context_graphs_required_ );

	max_context_neighbor_cutoff_ = 0;
	for ( uint ii = 1; ii <= num_context_graph_types; ++ii ) {
		if ( context_graphs_required_[ ii ] ) {
			ContextGraphOP example_graph = ContextGraphFactory::create_context_graph( ContextGraphType(ii) );
			if ( example_graph->neighbor_cutoff() > max_context_neighbor_cutoff_ ) {
				max_context_neighbor_cutoff_ = example_graph->neighbor_cutoff();
			}
		}
	}
	for ( uint ii = 1; ii <= n_score_types; ++ii ) {
		scores_present_[ ScoreType ( ii ) ] = scorefxn.has_zero_weight( ScoreType( ii ) ) ? 0 : 1;
	}

	/// Use the copy constructor to create a new object -- do not simply take a pointer to the existing object, which might change.
	if ( ! energy_method_options_ ) {
		energy_method_options_ = methods::EnergyMethodOptionsOP( new methods::EnergyMethodOptions( scorefxn.energy_method_options() ) );
	} else {
		*energy_method_options_ = scorefxn.energy_method_options();
	}
}


bool
operator==( ScoreFunctionInfo const & a, ScoreFunctionInfo const & b ) /* PHIL */ {

	if ( a.max_atomic_interaction_distance_ != b.max_atomic_interaction_distance_ ) return false;
	if ( a.max_context_neighbor_cutoff_ != b.max_context_neighbor_cutoff_ ) return false;
	for ( uint ii = 1; ii <= num_context_graph_types; ++ii ) {
		if ( a.context_graphs_required_[ ii ] != b.context_graphs_required_[ ii ] ) return false;
	}
	for ( uint ii = 1; ii <= n_score_types; ++ii ) {
		if ( a.scores_present_[ ScoreType( ii ) ] != b.scores_present_[ ScoreType ( ii ) ] ) return false;
	}
	// SFIs are not the same if either of them have been default constructed without originating from
	// an existing ScoreFunction.

	if ( ! a.energy_method_options_ || ! b.energy_method_options_ ) return false;
	if ( (*a.energy_method_options_) != (*b.energy_method_options_) ) return false;

	return true;
}

bool
ScoreFunctionInfo::requires_context_graph( ContextGraphType cgt ) const {
	return context_graphs_required_[ cgt ];
}


}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::ScoreFunctionInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( max_atomic_interaction_distance_ ) ); // Distance
	arc( CEREAL_NVP( max_context_neighbor_cutoff_ ) ); // Distance
	arc( CEREAL_NVP( context_graphs_required_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( scores_present_ ) ); // EnergyMap
	arc( CEREAL_NVP( energy_method_options_ ) ); // methods::EnergyMethodOptionsOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::ScoreFunctionInfo::load( Archive & arc ) {
	arc( max_atomic_interaction_distance_ ); // Distance
	arc( max_context_neighbor_cutoff_ ); // Distance
	arc( context_graphs_required_ ); // utility::vector1<_Bool>
	arc( scores_present_ ); // EnergyMap
	arc( energy_method_options_ ); // methods::EnergyMethodOptionsOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::ScoreFunctionInfo );
CEREAL_REGISTER_TYPE( core::scoring::ScoreFunctionInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_ScoreFunctionInfo )
#endif // SERIALIZATION
