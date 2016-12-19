// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/MMAtomTypeSet.cc
/// @brief
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/chemical/MMAtomTypeSet.hh>

// Project headers
#include <basic/Tracer.hh>

#include <core/chemical/MMAtomType.hh>

// C++ headers
#include <iostream>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer tr( "core.chemical" );


MMAtomTypeSet::MMAtomTypeSet( std::string const & name ):
	name_( name ),
	atom_type_index_(),
	atoms_()
{
}


MMAtomTypeSet::~MMAtomTypeSet() = default;


/// @details Initialize an MMAtomTypeSet from an external file "filename",
/// and set parameters and properties for each MMAtomType.
/// Refer to minirosetta_database_stock/chemical/mm_atom_type_sets/fa_standard/mm_atom_properties.txt
/// for file format
///
void
MMAtomTypeSet::read_file( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() );

	if ( !data.good() ) utility_exit_with_message( "Unable to open MM atom type set file: "+filename );

	// parse the header line
	utility::vector1< std::string > tags;
	{ // scope
		std::string line, tag;
		getline( data, line );
		std::istringstream l( line );
		l >> tag;
		if ( tag != "NAME" ) {
			utility_exit_with_message("MMAtomTypeSet::read_file: bad first line: "+ line );
		}
		l >> tag;
		while ( !l.fail() ) {
			tags.push_back( tag );
			l >> tag;
		}
	}

	// now parse the rest of the file
	int const ntags( tags.size() );
	{
		using namespace basic;

		std::string line/*, tag*/, name_wo_whitespace;
		while ( getline( data,line ) ) {
			std::istringstream l( line );
			l >> name_wo_whitespace;
			if ( name_wo_whitespace.find("#",0) == 0 ) continue; // skip comment lines
			if ( l.fail() ) {
				utility_exit_with_message("bad line: "+line);
			}

			MMAtomTypeOP mm_atom_type_ptr( new MMAtomType( name_wo_whitespace ) );

			// now parse the parameters
			for ( int i=1; i<= ntags; ++i ) {
				Real setting;
				l >> setting;
				mm_atom_type_ptr->set_parameter( tags[i], setting );
			}
			if ( l.fail() ) {
				utility_exit_with_message("bad line: "+line);
			}

			// add this to the list
			atoms_.push_back( mm_atom_type_ptr );
			if ( atom_type_index_.count( name_wo_whitespace ) ) {
				utility_exit_with_message("MMAtomTypeSet:: duplicate atom name "+name_wo_whitespace);
			}
			atom_type_index_[ name_wo_whitespace ] = atoms_.size();
			tr.Debug << "New MM atom type: " << name_wo_whitespace << std::endl;
		}
	}
}

/// @details This function iterates over each element in the atom_type_index_ map and
/// prints both keys. It is only used for debugging.
void
MMAtomTypeSet::print_all_types()
{
	for ( std::map< std::string, int >::const_iterator i = atom_type_index_.begin(), e = atom_type_index_.end(); i != e; ++i ) {
		std::cout << (*i).first << " " << (*i).second << std::endl;
	}
}

} // chemical
} // core

#ifdef SERIALIZATION

#include <core/chemical/ChemicalManager.hh>

template < class Archive >
void core::chemical::serialize_mm_atom_type_set( Archive & arc, core::chemical::MMAtomTypeSetCOP ptr )
{
	if ( ! ptr ) {
		bool ptr_is_nonnull( false );
		arc( CEREAL_NVP( ptr_is_nonnull ) );
	} else {
		bool ptr_is_nonnull( true );
		arc( CEREAL_NVP( ptr_is_nonnull ) );
		std::string typeset_name( ptr->name() ); // Assumes that the name can be used to extract it from the ChemicalManager
		arc( CEREAL_NVP( typeset_name ) );
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, core::chemical::serialize_mm_atom_type_set, core::chemical::MMAtomTypeSetCOP );

template < class Archive >
void core::chemical::deserialize_mm_atom_type_set( Archive & arc, core::chemical::MMAtomTypeSetCOP & ptr )
{
	bool ptr_is_nonnull( true ); arc( ptr_is_nonnull );
	if ( ptr_is_nonnull ) {
		std::string typeset_name;
		arc( typeset_name );
		ptr = core::chemical::ChemicalManager::get_instance()->mm_atom_type_set( typeset_name );
	} else {
		ptr = nullptr;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, core::chemical::deserialize_mm_atom_type_set, core::chemical::MMAtomTypeSetCOP & );

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_MMAtomTypeSet )
#endif // SERIALIZATION
