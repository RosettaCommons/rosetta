// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibraryFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary SingleResidueRotamerLibrarys
///         from a string --> SingleResidueRotamerLibraryCreator map
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

// Package headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleBasicRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryCreator.hh>

// Program Headers

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/BasicRotamerLibrarySpecification.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

namespace core {
namespace pack {
namespace rotamers {

static THREAD_LOCAL basic::Tracer TR("core.pack.rotamers.SingleResidueRotamerLibraryFactory");

#ifdef MULTI_THREADED
std::mutex SingleResidueRotamerLibraryFactory::cache_mutex_{};
#endif

void
SingleResidueRotamerLibraryFactory::factory_register( SingleResidueRotamerLibraryCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string err_msg = "Factory Name Conflict: Two or more SingleResidueRotamerLibraryCreators registered with the name " + creator->keyname();
		utility_exit_with_message( err_msg );
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool
SingleResidueRotamerLibraryFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

/// @details Will also check to make sure that the Creator is in the map and appropriately registered.
std::string
SingleResidueRotamerLibraryFactory::type_for_residuetype( core::chemical::ResidueType const & restype ) const {
	core::chemical::rotamers::RotamerLibrarySpecificationCOP const & rotlibspec( restype.rotamer_library_specification() );
	if ( rotlibspec ) {
		std::string type( rotlibspec->keyname() );

		// Check to make sure that it's appropriately registered.
		if ( ! has_type(type) ) {
			TR << "For residue " << restype.name() << ", can't find SingleResidueRotamerLibraryCreator for type '" << type << "'. Known types are: " << std::endl;
			TR << "\t\t";
			for ( CreatorMap::const_iterator iter( creator_map_.begin() ), iter_end( creator_map_.end() ); iter != iter_end; ++iter ) {
				TR << iter->first << ", ";
			}
			TR << std::endl;
			utility_exit_with_message( "Unable to find creator for SingleResidueRotamerLibrary type " + type );
		}

		return type;

	} else {
		return "";
	}
}

std::string
SingleResidueRotamerLibraryFactory::get_cachetag( core::chemical::ResidueType const & restype ) const {
	core::chemical::rotamers::RotamerLibrarySpecificationCOP const & rotlibspec( restype.rotamer_library_specification() );
	if ( rotlibspec ) {
		return rotlibspec->cache_tag( restype );
	} else {
		return "";
	}
}

/// @details For thread safety in accessing the cache, we use a "Software Transactional Memory"-type approach.
/// We lock the cache to see if we have something cached. If not, we release the lock and build the new SRRL.
/// This allows other threads to potentially access the shared cache for other residue types in the meantime.
/// We then re-lock the cache, check that the generated library is still needed, and then put it in the cache.
/// If another thread got in before us, we discard the work we did and use the previously made one.
/// (The cache is write-once.) This should work so long as cacheable SRRLs made from SRRLSpecifications with the
/// same type tag and cache string are functionally identical.
core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueRotamerLibraryFactory::get( core::chemical::ResidueType const & restype, bool forcebasic /* false */ ) const {
	std::string type( type_for_residuetype( restype ) );
	if ( type == "" ) {
		// A null pointer means that there isn't a rotamer library
		if ( forcebasic ) {
			return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
		} else {
			return SingleResidueRotamerLibraryCOP(0);
		}
	}
	std::string cachetag( get_cachetag(restype) );
	std::pair< std::string, std::string > cachepair(type,cachetag);
	if ( cachetag.size() ) {
#ifdef MULTI_THREADED
		std::lock_guard<std::mutex> lock(cache_mutex_); // mutex will release when this scope exits
#endif
		if ( cache_.count( cachepair ) ) {
			if ( ! cache_[ cachepair ] && forcebasic ) {
				return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
			} else {
				return cache_[ cachepair ];
			}
		}
	}

	CreatorMap::const_iterator entry(creator_map_.find( type ) );
	assert( entry != creator_map_.end() );
	core::pack::rotamers::SingleResidueRotamerLibraryCOP library( entry->second->create( restype ) );

	// Some creators (e.g. Dunbrack for Gly/Ala) will return empty residue library pointers
	// We can cache a null pointer (don't cache a SingleBasicRotamerLibrary, as not all calls will forcebasic)

	if ( cachetag.size() ) {
#ifdef MULTI_THREADED
		std::lock_guard<std::mutex> lock(cache_mutex_); // mutex will release when this scope exits
#endif
		if ( ! cache_.count( cachepair ) ) {
			cache_[ cachepair ] = library;
		}
		if ( ! cache_[ cachepair ] && forcebasic ) {
			return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
		} else {
			return cache_[ cachepair ]; // Catch case where another thread already assigned cachepair
		}
	}
	if ( ! library && forcebasic ) {
		return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
	} else {
		return library;
	}
}

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueRotamerLibraryFactory::get( core::chemical::ResidueType const & restype, core::conformation::Residue const & residue ) const {
	std::string cachetag( get_cachetag(restype) );
	if ( cachetag.size() ) {
		// Cachable means that you can't alter the library creator based on the Residue being passed
		return get( restype );
	}
	std::string type( type_for_residuetype( restype ) );
	CreatorMap::const_iterator entry(creator_map_.find( type ) );
	assert( entry != creator_map_.end() );
	return entry->second->create( restype, residue );
}

SingleResidueRotamerLibraryFactory::SingleResidueRotamerLibraryFactory()
{}

} //namespace rotamers
} //namespace pack
} //namespace core
