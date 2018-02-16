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
#include <utility/thread/threadsafe_creation.hh>
#include <basic/Tracer.hh>

namespace core {
namespace pack {
namespace rotamers {

static basic::Tracer TR("core.pack.rotamers.SingleResidueRotamerLibraryFactory");

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
			for ( auto const & iter : creator_map_ ) {
				TR << iter.first << ", ";
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

/// @brief Get the SingleResidueRotamerLibrary coresponding to the given ResidueType.
/// If forcebasic is true, a SingleBasicRotamerLibrary will be returned instead of a null pointer.
/// @details If a cachetag exists for the ResidueType, store the information in the cache.
/// If not, always regenerate the RotamerLibrary
/// @author Rocco Moretti
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- implemented more efficient locking scheme.
core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueRotamerLibraryFactory::get( core::chemical::ResidueType const & restype, bool forcebasic /* false */ ) const {
	std::string type( type_for_residuetype( restype ) );
	if ( type == "" ) {
		// A null pointer means that there isn't a rotamer library
		if ( forcebasic ) {
			return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
		} else {
			return SingleResidueRotamerLibraryCOP(nullptr);
		}
	}
	std::string cachetag( get_cachetag(restype) );
	std::pair< std::string, std::string > cachepair(type,cachetag);

	if ( cachetag.size() ) {
		{ // Scope for read lock
#ifdef MULTI_THREADED
			utility::thread::ReadLockGuard readlock(cache_mutex_); // mutex will release when this scope exits
#endif
			if ( cache_.count( cachepair ) ) {
				if ( ! cache_[ cachepair ] && forcebasic ) {
					return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
				} else {
					return cache_[ cachepair ];
				}
			}
		} // End scope for read lock

		{ // Scope for write lock
			// We grab the write lock early, as we want to avoid memory spikes when
			// large numbers of threads all attempt to speculatively create the library at the same time.
#ifdef MULTI_THREADED
			utility::thread::WriteLockGuard writelock(cache_mutex_); // mutex will release when this scope exits
#endif
			//Check again -- creation might have happened between read guard release and write guard acquire:
			if ( cache_.count( cachepair ) ) {
				if ( ! cache_[ cachepair ] && forcebasic ) {
					return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
				} else {
					return cache_[ cachepair ];
				}
			}

			auto entry( creator_map_.find( type ) );
			debug_assert( entry != creator_map_.end() );
			core::pack::rotamers::SingleResidueRotamerLibraryCOP library( entry->second->create( restype ) );

			// Some creators (e.g. Dunbrack for Gly/Ala) will return empty residue library pointers
			// We can cache a null pointer (don't cache a SingleBasicRotamerLibrary, as not all calls will forcebasic)
			if ( ! cache_.count( cachepair ) ) {
				cache_[ cachepair ] = library;
			}

			if ( ! library && forcebasic ) {
				return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
			} else {
				return library;
			}
		} // End scope for write lock
	} else { // No cache tag - no locking needed
		auto entry( creator_map_.find( type ) );
		debug_assert( entry != creator_map_.end() );
		core::pack::rotamers::SingleResidueRotamerLibraryCOP library( entry->second->create( restype ) );

		if ( ! library && forcebasic ) {
			return SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
		} else {
			return library;
		}
	}
	return nullptr; // Should never get here.
}

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueRotamerLibraryFactory::get( core::chemical::ResidueType const & restype, core::conformation::Residue const & residue ) const {
	std::string cachetag( get_cachetag(restype) );
	if ( cachetag.size() ) {
		// Cachable means that you can't alter the library creator based on the Residue being passed
		return get( restype );
	}
	std::string type( type_for_residuetype( restype ) );
	auto entry(creator_map_.find( type ) );
	debug_assert( entry != creator_map_.end() );
	return entry->second->create( restype, residue );
}

SingleResidueRotamerLibraryFactory::SingleResidueRotamerLibraryFactory() :
	creator_map_(),
#ifdef MULTI_THREADED
	cache_mutex_(),
#endif
	cache_()
	//TODO: init vars here.
{}

} //namespace rotamers
} //namespace pack
} //namespace core
