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

#include <core/chemical/ResidueTypeSetCache.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.chemical.ResidueTypeSetCache" );

namespace core {
namespace chemical {

//Constructor
ResidueTypeSetCache::ResidueTypeSetCache( ResidueTypeSet const & rsd_type_set ):
	rsd_type_set_( rsd_type_set )
{}

//Destructor
ResidueTypeSetCache::~ResidueTypeSetCache()
{}

/// @details Main accessor function into ResidueTypeSetCache
ResidueTypeCOP
ResidueTypeSetCache::name_map( std::string const & name_in ) const
{
	std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.find( name_in );
	runtime_assert( it != name_map_.end() );
	return it->second;
}

void
ResidueTypeSetCache::clear_cached_maps()
{
	aa_map_.clear();
	name3_map_.clear();
	cached_aa_variants_map_.clear();
}

void
ResidueTypeSetCache::add_residue_type( ResidueTypeCOP residue_type )
{
	runtime_assert_string_msg( name_map_.find( residue_type->name() ) == name_map_.end(), "Error in core::chemical::ResidueTypeSetCache::add_residue_type(): Attempting to add a new residue type, but it already exists in the cache.  (Did you load a .params file with the -extra_res_fa commandline option that was already listed in residue_types.txt, perhaps?)" );
	name_map_[ residue_type->name() ] = residue_type;
	//  clear_cached_maps(); // no can't do this
}

void
ResidueTypeSetCache::remove_residue_type( std::string const & name )
{
	std::map< std::string, ResidueTypeCOP >::iterator it = name_map_.find( name );
	runtime_assert( it != name_map_.end() );
	name_map_.erase( it );
	clear_cached_maps();
}

void
ResidueTypeSetCache::update_residue_type( ResidueTypeCOP residue_type_original, ResidueTypeCOP residue_type_new )
{
	std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.find( residue_type_original->name() );
	runtime_assert( it != name_map_.end() );
	runtime_assert( residue_type_original->name() == residue_type_new->name() );
	name_map_[ residue_type_original->name() ] = residue_type_new;
	clear_cached_maps();
}

bool
ResidueTypeSetCache::has_generated_residue_type( ResidueTypeCOP residue_type ) const {
	std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.find( residue_type->name() );
	if ( it == name_map_.end() ) return false;
	runtime_assert( it->second == residue_type );
	return true;
}

bool
ResidueTypeSetCache::has_generated_residue_type( std::string const & rsd_name ) const {
	return ( name_map_.find( rsd_name ) != name_map_.end() );
}

ResidueTypeCOPs
ResidueTypeSetCache::generated_residue_types() {
	ResidueTypeCOPs residue_types;
	for ( std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.begin();
			it != name_map_.end(); ++it ) {
		residue_types.push_back( it->second );
	}
	return residue_types;
}

void
ResidueTypeSetCache::add_prohibited( std::string const & rsd_name )
{
	prohibited_types_.insert( rsd_name );
}

bool
ResidueTypeSetCache::is_prohibited( std::string const & rsd_name ) const
{
	return prohibited_types_.count( rsd_name );
}

ResidueTypeCOPs
ResidueTypeSetCache::get_all_types_with_variants_aa( AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions )
{
	AA_VariantsExceptions query( std::make_pair( aa, std::make_pair( variants, exceptions ) ) );
	if ( cached_aa_variants_map_.find( query ) == cached_aa_variants_map_.end() ) {
		cached_aa_variants_map_[ query ] = ResidueTypeFinder( rsd_type_set_ ).aa( aa ).variants( variants ).variant_exceptions( exceptions ).get_all_possible_residue_types();
	}
	return cached_aa_variants_map_[ query ];
}


} //chemical
} //core
