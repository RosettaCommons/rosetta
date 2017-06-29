// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueTypeSetCache.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

// Unit headers
#include <core/chemical/ResidueTypeSetCache.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/Patch.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tools/make_vector1.hh>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.chemical.ResidueTypeSetCache" );

namespace core {
namespace chemical {

//Constructor
ResidueTypeSetCache::ResidueTypeSetCache( ResidueTypeSet const & rsd_type_set ):
	rsd_type_set_( rsd_type_set )
{
#ifdef MULTI_THREADED
	read_write_mutex_.reset( new utility::thread::ReadWriteMutex );
#endif
}

//Destructor
ResidueTypeSetCache::~ResidueTypeSetCache() = default;

// @details A shallow copy of the cache -
ResidueTypeSetCacheOP
ResidueTypeSetCache::clone( ResidueTypeSet const & rsd_type_set ) const {
	ResidueTypeSetCacheOP cloned( new ResidueTypeSetCache(rsd_type_set) );

	// This is a shallow copy - (immutable) ResidueTypes can be shared by multiple ResidueTypeSets
	cloned->name_map_ = name_map_;
	cloned->prohibited_types_ = prohibited_types_;

	cloned->clear_cached_maps(); // clear the caches, and allow them to regenerate.

	return cloned;
}

/// @details Main accessor function into ResidueTypeSetCache.  Must be read locked
/// by calling function.
ResidueTypeCOP
ResidueTypeSetCache::name_map( std::string const & name_in ) const
{
	auto it = name_map_.find( name_in );
	runtime_assert( it != name_map_.end() );
	return it->second;
}

/// @details Must be write-locked by calling function
void
ResidueTypeSetCache::add_residue_type( ResidueTypeCOP residue_type )
{
	if ( name_map_.find( residue_type->name() ) != name_map_.end() ) {
		TR.Error << "Residue type " << residue_type->name() << " is already in thge ResidueTypeSetCache!" << std::endl;
		utility_exit_with_message( "Error in core::chemical::ResidueTypeSetCache::add_residue_type(): Attempting to add a new residue type, but it already exists in the cache.  (Did you load a .params file with the -extra_res_fa commandline option that was already listed in residue_types.txt, perhaps?)" );
	}
	name_map_[ residue_type->name() ] = residue_type;
	//  clear_cached_maps(); // no can't do this
}

/// @details Danger Danger Danger.  Unsafe in a multi-threaded environment.
void
ResidueTypeSetCache::remove_residue_type( std::string const & name )
{
	auto it = name_map_.find( name );
	runtime_assert( it != name_map_.end() );
	name_map_.erase( it );
	clear_cached_maps();
}

/// @details Danger Danger Danger. Unsafe in a multi-threaded environment
void
ResidueTypeSetCache::update_residue_type( ResidueTypeCOP residue_type_original, ResidueTypeCOP residue_type_new )
{
	std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.find( residue_type_original->name() );
	runtime_assert( it != name_map_.end() );
	runtime_assert( residue_type_original->name() == residue_type_new->name() );
	name_map_[ residue_type_original->name() ] = residue_type_new;
	clear_cached_maps();
}

/// @details Must be read-locked by calling function
bool
ResidueTypeSetCache::has_generated_residue_type( ResidueTypeCOP residue_type ) const {
	auto it = name_map_.find( residue_type->name() );
	if ( it == name_map_.end() ) return false;
	// The provided residue type could come from a different ResidueTypeSet, so it could
	// be a different object, even if the names match.
	return it->second == residue_type;
}

/// @details Must be read-locked by calling function
bool
ResidueTypeSetCache::has_generated_residue_type( std::string const & rsd_name ) const {
	return ( name_map_.find( rsd_name ) != name_map_.end() );
}


///// @details Must be read-locked by calling function
//ResidueTypeCOPs
//ResidueTypeSetCache::generated_residue_types() {
// ResidueTypeCOPs residue_types;
// for ( std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.begin();
//   it != name_map_.end(); ++it ) {
//  residue_types.push_back( it->second );
// }
// return residue_types;
//}

/// @details Needs to be write-locked by calling function
void
ResidueTypeSetCache::add_prohibited( std::string const & rsd_name )
{
	prohibited_types_.insert( rsd_name );
}

/// @details Needs to be read-locked by calling function
bool
ResidueTypeSetCache::is_prohibited( std::string const & rsd_name ) const
{
	return prohibited_types_.count( rsd_name );
}

bool
ResidueTypeSetCache::all_types_with_variants_aa_already_cached(
	AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions
) const
{
	AA_VariantsExceptions query( std::make_pair( aa, std::make_pair( variants, exceptions ) ) );
	return cached_aa_variants_map_.find( query ) != cached_aa_variants_map_.end();
}

void
ResidueTypeSetCache::cache_all_types_with_variants_aa(
	AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions,
	ResidueTypeCOPs cached_types
)
{
	AA_VariantsExceptions query( std::make_pair( aa, std::make_pair( variants, exceptions ) ) );
	cached_aa_variants_map_[ query ] = cached_types;
}

ResidueTypeCOPs
ResidueTypeSetCache::retrieve_all_types_with_variants_aa(
	AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions
) const
{
	AA_VariantsExceptions query( std::make_pair( aa, std::make_pair( variants, exceptions ) ) );
	return cached_aa_variants_map_.find( query )->second;
}

/// @details Danger danger danger.  Unsafe in a multithreaded environment.
void
ResidueTypeSetCache::clear_cached_maps()
{
	cache_up_to_date_ = false;
	//aa_map_.clear();
	//name3_map_.clear();
	cached_aa_variants_map_.clear();
	name3_generated_by_base_residue_name_.clear();
	interchangeability_group_generated_by_base_residue_name_.clear();
}

/// @brief information on residue types whose name3's can be changed by patches.
std::map< std::string, std::set< std::string > > const &
ResidueTypeSetCache::name3_generated_by_base_residue_name() {
	if ( ! cache_up_to_date_ ) { regenerate_cached_maps(); }
	return name3_generated_by_base_residue_name_;
}

/// @brief interchangeability groups that appear upon patch application.
std::map< std::string, std::set< std::string > > const &
ResidueTypeSetCache::interchangeability_group_generated_by_base_residue_name() {
	if ( ! cache_up_to_date_ ) { regenerate_cached_maps(); }
	return interchangeability_group_generated_by_base_residue_name_;
}

/// @details following assumes that all new name3 and interchangeability groups for residue types
///    can be discovered by applying patches to base residue types -- i.e. on the 'first patch'.
///    Probably should set up a runtime_assert in ResidueTypeFinder to check this assumption.
void
ResidueTypeSetCache::regenerate_cached_maps() {

	name3_generated_by_base_residue_name_.clear();
	interchangeability_group_generated_by_base_residue_name_.clear();

	for ( PatchCOP p : rsd_type_set_.patches() ) {
		for ( ResidueTypeCOP rsd_type : rsd_type_set_.base_residue_types() ) {
			if ( p->applies_to( *rsd_type ) ) {

				// Check if any patches change name3 of residue_types.
				// Such patches must be applicable to base residue types -- probably should
				// check this somewhere (e.g. ResidueTypeFinder).
				std::string const & new_name3 = p->generates_new_name3( *rsd_type );
				if ( new_name3.size() > 0 && rsd_type->name3() != new_name3 ) {
					name3_generated_by_base_residue_name_[ rsd_type->name() ].insert( new_name3 );
				}

				// similarly, check if any patches create interchangeability_groups from this base residue type.
				// Used by ResidueTypeFinder.
				std::string const & interchangeability_group = p->generates_interchangeability_group( *rsd_type );
				if ( interchangeability_group.size() > 0 ) {
					interchangeability_group_generated_by_base_residue_name_[ rsd_type->name() ].insert( interchangeability_group );
				}

			}
		}
	}

	cache_up_to_date_ = true;
}

#ifdef MULTI_THREADED
//std::recursive_mutex & ResidueTypeSetCache::cache_mutex() { return cache_mutex_; }

utility::thread::ReadWriteMutex &
ResidueTypeSetCache::read_write_mutex() { return *read_write_mutex_; }

#endif

} //chemical
} //core

#ifdef    SERIALIZATION

template< class Archive >
void
core::chemical::ResidueTypeSetCache::save( Archive & arc ) const {
	arc( ::cereal::make_nvp( "rsd_type_set_", rsd_type_set_.get_self_ptr() ) );
	arc( CEREAL_NVP( name_map_ ) );
	arc( CEREAL_NVP( prohibited_types_ ) );

	// All the cached data is exempt, as it can be regenerated later.
	// EXEMPT cached_aa_variants_map_
	// EXEMPT cache_up_to_date_
	// ... yes, even cache_up_to_date_, because we're forcing it to false on reconstruction
	// EXEMPT name3_generated_by_base_residue_name_ interchangeability_group_generated_by_base_residue_name_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueTypeSetCache::load_and_construct( Archive & arc, cereal::construct< core::chemical::ResidueTypeSetCache > & construct ) {
	using namespace core::chemical;

	ResidueTypeSetCOP parent;
	arc( parent ); // EXEMPT rsd_type_set_

	construct( *parent );
	arc( construct->name_map_ );
	arc( construct->prohibited_types_ );

	construct->cache_up_to_date_ = false;

	// All the cached data is exempt, as it can be regenerated later.
	// EXEMPT cached_aa_variants_map_
	// EXEMPT cache_up_to_date_
	// EXEMPT name3_generated_by_base_residue_name_ interchangeability_group_generated_by_base_residue_name_
}

SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE( core::chemical::ResidueTypeSetCache );
CEREAL_REGISTER_TYPE( core::chemical::ResidueTypeSetCache )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ResidueTypeSetCache )
#endif // SERIALIZATION
