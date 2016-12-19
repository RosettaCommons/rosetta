// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for reading in the orbital type properties
///
/// @details
/// This class reads in the orbital_properties.txt file which contains the "chemical" information for orbitals.
/// This does not contain the actual properties, but sets the properties through the OrbitalType class.
/// This class is called by the ChemicalManager. Modeled off of atomtypeset.
///
///
///
/// @author
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////


#include <fstream>

#include <basic/Tracer.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/orbitals/OrbitalType.hh>
#include <utility/io/izstream.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace orbitals {

/// @details Auto-generated virtual destructor
OrbitalTypeSet::~OrbitalTypeSet() {}

static THREAD_LOCAL basic::Tracer TR( "core.chemical.orbitals" );

OrbitalTypeSet::OrbitalTypeSet( std::string const & directory, std::string const & name)
: name_( name ), directory_( directory )
{
	read_file( directory + "/orbital_properties.txt" );

	//currently commented out. This is if someone down the line wants to add extras, like
	//in the atomtype extras.
	/* std::ifstream data( ( directory+"/extras.txt" ).c_str() );
	if ( data.good() ) { // add extra data
	std::string line;
	while( getline( data, line ) ) {
	if ( line.size() && line[0] == '#' ) continue;
	add_parameters_from_file( directory+"/"+line );
	}
	}
	data.close();*/

}


//function take from AtomTypeSet.cc. Modified slightly for OrbitalTypeSet.cc.
//basically reads the file orbital_properties.txt found in your directory.
//
void OrbitalTypeSet::read_file(std::string const & filename)
{
	utility::io::izstream data( filename.c_str() );

	if ( !data.good() ) utility_exit_with_message( "Unable to open atomset file: "+filename );

	// parse the header line currently should look like orbital_type atom_type hybridization distance
	utility::vector1< std::string > tags;
	{ // scope
		std::string line, tag, tag2;
		getline( data, line );
		std::istringstream l( line );
		l >> tag >> tag2;
		if ( tag != "orbital_type" || tag2 != "atom_type" ) {
			utility_exit_with_message("AtomTypeSet::read_file: bad first line: "+ line );
		}
		l >> tag;
		while ( !l.fail() ) {
			tags.push_back( tag );
			l >> tag;
		}
	}

	// now parse the rest of the file
	core::Size const ntags( tags.size() );
	{
		using namespace basic;

		std::string line, tag, name_wo_whitespace;
		while ( getline( data,line ) ) {
			std::istringstream l( line );
			l >> name_wo_whitespace;
			if ( l.fail() || name_wo_whitespace.find("#",0) == 0 ) continue; // skip comment,blank lines
			l >> tag;
			if ( l.fail() || tag.size() < 1 ) {
				utility_exit_with_message("bad line: "+line);
			}

			//   std::string const name( line.substr(0,4) );
			std::string atom_type_name( tag );
			OrbitalTypeOP orbital_type_ptr( new OrbitalType( name_wo_whitespace, atom_type_name ) );

			// now parse the parameters
			for ( Size i=1; i<= ntags; ++i ) {
				Real setting;
				l >> setting;
				orbital_type_ptr->set_parameter( tags[i], setting );
			}
			if ( l.fail() ) {
				utility_exit_with_message("bad line: "+line);
			}

			// now parse the properties
			l >> tag;
			while ( !l.fail() && tag.length() != 0 && tag[0] != '#' ) { //tag.find("#",0) != 0 ) {
				orbital_type_ptr->set_property( tag, true );
				l >> tag;
			}

			// add this to the list
			orbitals_.push_back( orbital_type_ptr );
			//  atom_type_index_[ name ] = atoms_.size();
			if ( orbital_type_index_.count( name_wo_whitespace ) ) {
				utility_exit_with_message("AtomTypeSet:: duplicate atom name "+name_wo_whitespace);
			}
			orbital_type_index_[ name_wo_whitespace ] = orbitals_.size();
			TR.Debug << "New atom type: " << name_wo_whitespace << ' ' << atom_type_name << std::endl; //std::endl;
		}
	} // scope


}

int
OrbitalTypeSet::orbital_type_index( std::string const & orbital_type_name ) const
{
	std::map< std::string, int >::const_iterator
		iter( orbital_type_index_.find( orbital_type_name ) );
	if ( iter == orbital_type_index_.end() ) {
		utility_exit_with_message("unrecognized orbital type name "+orbital_type_name);
	}
	return iter->second;
}


int
OrbitalTypeSet::orbital_type_index( std::string & orbital_type_name ) const
{
	std::map< std::string, int >::const_iterator
		iter( orbital_type_index_.find( orbital_type_name ) );
	if ( iter == orbital_type_index_.end() ) {
		utility_exit_with_message("unrecognized orbital type name "+orbital_type_name);
	}
	return iter->second;
}

}
}
}

#ifdef SERIALIZATION
#include <core/chemical/ChemicalManager.hh>

template < class Archive >
void core::chemical::orbitals::serialize_orbital_type_set( Archive & arc, core::chemical::orbitals::OrbitalTypeSetCOP ptr )
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
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, core::chemical::orbitals::serialize_orbital_type_set, core::chemical::orbitals::OrbitalTypeSetCOP );

template < class Archive >
void core::chemical::orbitals::deserialize_orbital_type_set( Archive & arc, core::chemical::orbitals::OrbitalTypeSetCOP & ptr )
{
	bool ptr_is_nonnull( true ); arc( ptr_is_nonnull );
	if ( ptr_is_nonnull ) {
		std::string typeset_name;
		arc( typeset_name );
		ptr = core::chemical::ChemicalManager::get_instance()->orbital_type_set( typeset_name );
	} else {
		ptr = nullptr;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, core::chemical::orbitals::deserialize_orbital_type_set, core::chemical::orbitals::OrbitalTypeSetCOP & );

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_orbitals_OrbitalTypeSet )
#endif // SERIALIZATION
