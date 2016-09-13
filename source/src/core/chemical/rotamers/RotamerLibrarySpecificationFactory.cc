// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/RotamerLibrarySpecificationFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary RotamerLibrarySpecifications
///         from a string --> RotamerLibrarySpecificationCreator map
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/RotamerLibrarySpecificationFactory.hh>

// Package headers
#include <core/chemical/rotamers/RotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecificationCreator.hh>

// Program Headers

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

namespace core {
namespace chemical {
namespace rotamers {

static THREAD_LOCAL basic::Tracer TR("core.chemical.rotamers.RotamerLibrarySpecificationFactory");

RotamerLibrarySpecificationFactory *
RotamerLibrarySpecificationFactory::create_singleton_instance()
{
	return new RotamerLibrarySpecificationFactory;
}

void
RotamerLibrarySpecificationFactory::factory_register( RotamerLibrarySpecificationCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string err_msg = "Factory Name Conflict: Two or more RotamerLibrarySpecificationCreators registered with the name " + creator->keyname();
		utility_exit_with_message( err_msg );
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool
RotamerLibrarySpecificationFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

RotamerLibrarySpecificationCreatorCOP
RotamerLibrarySpecificationFactory::get_creator( std::string const & tag ) const {
	CreatorMap::const_iterator entry( creator_map_.find( tag ) );
	if ( entry == creator_map_.end() ) {
		TR << "Can't find RotamerLibrarySpecificationCreators of type " << tag << ". Known types are: " << std::endl;
		TR << "\t";
		for ( CreatorMap::const_iterator iter( creator_map_.begin() ), iter_end( creator_map_.end() ); iter != iter_end; ++iter ) {
			TR << iter->first << ", ";
		}
		TR << std::endl;
		utility_exit_with_message("Cannot find RotamerLibrarySpecification of type " + tag );
	}
	return entry->second;
}

core::chemical::rotamers::RotamerLibrarySpecificationOP
RotamerLibrarySpecificationFactory::get( std::string const & tag) const {
	return get_creator( tag )->create();
}

core::chemical::rotamers::RotamerLibrarySpecificationOP
RotamerLibrarySpecificationFactory::get( std::string const & tag, std::istream & parameters ) const {
	return get_creator( tag )->create( parameters );
}


RotamerLibrarySpecificationFactory::RotamerLibrarySpecificationFactory()
{}

} //namespace rotamers
} //namespace chemical
} //namespace core
