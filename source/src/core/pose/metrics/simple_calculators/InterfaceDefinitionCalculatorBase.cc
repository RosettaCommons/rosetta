// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/simple_metrics/InterfaceDefinitionCalculatorBase.cc
/// @brief  InterfaceDefinitionCalculatorBase class
/// @author John Karanicolas

// Unit headers
#include <core/pose/metrics/simple_calculators/InterfaceDefinitionCalculatorBase.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace metrics {
namespace simple_calculators {


InterfaceDefinitionCalculator::InterfaceDefinitionCalculator( Size const chain1_number, Size const chain2_number ) :
	StructureDependentCalculator(),
	chain1_number_(chain1_number),
	chain2_number_(chain2_number),
	got_chain_numbers_(true)
{
}

InterfaceDefinitionCalculator::InterfaceDefinitionCalculator( char const chain1_letter, char const chain2_letter ) :
	StructureDependentCalculator(),
	chain1_letter_(chain1_letter),
	chain2_letter_(chain2_letter),
	got_chain_numbers_(false)
{
}


void InterfaceDefinitionCalculator::verify_chain_setup( pose::Pose const & pose ) {
	if ( ! got_chain_numbers_ ) {
		chain1_number_ = chain_letter_to_number( pose, chain1_letter_ );
		chain2_number_ = chain_letter_to_number( pose, chain2_letter_ );
		got_chain_numbers_ = true;
	}
	fill_in_chain_terminii( pose );
	return;
}


// helper method to convert a single chain letter to a single chain number.
// returns 0 if chain letter not found
Size InterfaceDefinitionCalculator::chain_letter_to_number( pose::Pose const & pose, char const chain_id ) {
	char temp_letter_ = chain_id;
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.pdb_info()->chain( i ) == temp_letter_ ) {
			return pose.chain( i );
		} else continue;
	}
	return 0; // if this happens, we haven't found our chain letter in the pose
}


// figures out chain terminii from chain numbers
void InterfaceDefinitionCalculator::fill_in_chain_terminii( pose::Pose const & pose ) {
	ch1_begin_num_ = pose.conformation().chain_begin( chain1_number_ );
	ch1_end_num_ = pose.conformation().chain_end( chain1_number_ );
	ch2_begin_num_ = pose.conformation().chain_begin( chain2_number_ );
	ch2_end_num_ = pose.conformation().chain_end( chain2_number_ );
	return;
}

} // simple_calculators
} // metrics
} // pose
} // core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator::InterfaceDefinitionCalculator() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( ch1_begin_num_ ) ); // core::Size
	arc( CEREAL_NVP( ch1_end_num_ ) ); // core::Size
	arc( CEREAL_NVP( ch2_begin_num_ ) ); // core::Size
	arc( CEREAL_NVP( ch2_end_num_ ) ); // core::Size
	arc( CEREAL_NVP( chain1_number_ ) ); // core::Size
	arc( CEREAL_NVP( chain2_number_ ) ); // core::Size
	arc( CEREAL_NVP( chain1_letter_ ) ); // char
	arc( CEREAL_NVP( chain2_letter_ ) ); // char
	arc( CEREAL_NVP( got_chain_numbers_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( ch1_begin_num_ ); // core::Size
	arc( ch1_end_num_ ); // core::Size
	arc( ch2_begin_num_ ); // core::Size
	arc( ch2_end_num_ ); // core::Size
	arc( chain1_number_ ); // core::Size
	arc( chain2_number_ ); // core::Size
	arc( chain1_letter_ ); // char
	arc( chain2_letter_ ); // char
	arc( got_chain_numbers_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator );
CEREAL_REGISTER_TYPE( core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_metrics_simple_calculators_InterfaceDefinitionCalculatorBase )
#endif // SERIALIZATION
