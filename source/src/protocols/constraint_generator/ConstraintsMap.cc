// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ConstraintsMap.cc
/// @brief Cacheable data map to cache constraint pointers in the pose.
/// @author Tom Linsky (tlinsky@uw.edu)

#include <protocols/constraint_generator/ConstraintsMap.hh>

#include <basic/Tracer.hh>

#ifdef    SERIALIZATION
// Required for serializing Constraints
#include <core/scoring/constraints/Constraint.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh> // 13
#include <utility/vector1.srlz.hh> // 11

// Numeric serialization headers
//#include <numeric/xyz.serialization.hh> // 7

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION
static basic::Tracer TR( "protocols.constraint_generator.ConstraintsMap" );


namespace protocols {
namespace constraint_generator {

ConstraintsMap::ConstraintsMap():
	basic::datacache::CacheableData(),
	cst_map_()
{}

ConstraintsMap::~ConstraintsMap() = default;

basic::datacache::CacheableDataOP
ConstraintsMap::clone() const
{
	return basic::datacache::CacheableDataOP( new ConstraintsMap( *this ) );
}

/// @brief Insert csts into the ConstraintsMap under the name given.
/// @param[in] name Map key name under which constraints will be stored
/// @param[in] csts Constraints to store
/// @returns ConstraintsMap::iterator to new map item.
ConstraintsMap::iterator
ConstraintsMap::insert( std::string const & name, ConstraintCOPs const & csts )
{
	std::pair< std::string, ConstraintCOPs > const pair( name, csts );
	return cst_map_.insert( pair ).first;
}

void
ConstraintsMap::erase( iterator const & erase_me )
{
	cst_map_.erase( erase_me );
}

ConstraintsMap::iterator
ConstraintsMap::find( std::string const & name )
{
	return cst_map_.find( name );
}

ConstraintsMap::const_iterator
ConstraintsMap::find( std::string const & name ) const
{
	return cst_map_.find( name );
}

ConstraintsMap::iterator
ConstraintsMap::begin()
{
	return cst_map_.begin();
}

ConstraintsMap::const_iterator
ConstraintsMap::begin() const
{
	return cst_map_.begin();
}

ConstraintsMap::iterator
ConstraintsMap::end()
{
	return cst_map_.end();
}

ConstraintsMap::const_iterator
ConstraintsMap::end() const
{
	return cst_map_.end();
}

std::string
ConstraintsMap::valid_names_string() const
{
	std::stringstream stream;
	utility::vector1< std::string > const names = valid_names();
	for ( auto const & name : names ) {
		if ( !stream.str().empty() ) stream << ", ";
		stream << name;
	}
	return stream.str();
}

utility::vector1< std::string >
ConstraintsMap::valid_names() const
{
	utility::vector1< std::string > names;
	for ( auto const & pair : cst_map_ ) {
		names.push_back( pair.first );
	}
	return names;
}

} //protocols
} //constraint_generator


#ifdef    SERIALIZATION // 1

typedef protocols::constraint_generator::ConstraintsMap ConstraintsMap;

template< class Archive >
void
ConstraintsMap::save( Archive & arc ) const
{
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) ); // 3
	arc( CEREAL_NVP( cst_map_ ) ); // 6
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
ConstraintsMap::load( Archive & arc )
{
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) ); // 5

	typedef std::map< std::string, core::scoring::constraints::ConstraintOPs > NameToConstraintOPs;
	NameToConstraintOPs local_cst_map;
	arc( local_cst_map );
	for ( NameToConstraintOPs::const_iterator pair=local_cst_map.begin(); pair!=local_cst_map.end(); ++pair ) {
		cst_map_[ pair->first ] = pair->second;
	}
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::constraint_generator::ConstraintsMap ); // 12
CEREAL_REGISTER_TYPE( ConstraintsMap ) // 14

CEREAL_REGISTER_DYNAMIC_INIT( protocols_constraint_generator_ConstraintsMap )
#endif // SERIALIZATION
