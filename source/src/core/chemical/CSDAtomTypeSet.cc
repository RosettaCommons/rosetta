// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/CSDAtomTypeSet.cc
/// @brief
/// @author Ian W. Davis

// Unit headers
#include <core/chemical/CSDAtomTypeSet.hh>

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>
#include <iostream>

namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical");


CSDAtomTypeSet::CSDAtomTypeSet():
	atom_type_index_(),
	atoms_()
{}


CSDAtomTypeSet::~CSDAtomTypeSet() {}

Size
CSDAtomTypetSet::n_atomtypes() const {
	return atoms_.size();
}

bool
CSDAtomTypeSetcontains_atom_type( std::string const & atom_type_name ) const {
	std::map< std::string, int >::const_iterator
		iter( atom_type_index_.find( atom_type_name ) );
	return iter != atom_type_index_.end();
}


/// @details Initialize an CSDAtomTypeSet from an external file "filename",
/// and set parameters and properties for each CSDAtomType.
/// Refer to minirosetta_database_stock/chemical/mm_atom_type_sets/fa_standard/mm_atom_properties.txt
/// for file format
void
CSDAtomTypeSet::read_file( std::string const & filename ) {
	using namespace basic;

	utility::io::izstream data( filename.c_str() );

	if ( !data.good() ) utility_exit_with_message( "Unable to open CSDAtomTypeSet file: "+filename );

	std::string line, tag;
	while ( getline( data,line ) ) {
		std::istringstream l( line );
		l >> tag;
		if ( l.fail() ) {
			utility_exit_with_message("bad line: "+line);
		}

		std::string  name(tag); //( line.substr(0,4) );
		CSDAtomType* csd_atom_type_ptr( new CSDAtomType( name ) );


		// add this to the list
		atoms_.push_back( csd_atom_type_ptr );
		atom_type_index_[ name ] = atoms_.size();
		tr.Debug << "New CSD atom type: " << name << std::endl;
	}

	tr.Debug << "CSD atoms types added " << atoms_.size() << std::endl;
}

/// @details This function iterates over each element in the atom_type_index_ map and
/// prints both keys. It is only used for debugging.
void
CSDAtomTypeSet::print_all_types() {
	for( std::map< std::string, int >::const_iterator i = atom_type_index_.begin(), e = atom_type_index_.end(); i != e; ++i )
		{
			std::cout << (*i).first << " " << (*i).second << std::endl;
		}
}

int
CSDAtomTypeSet::atom_type_index( std::string const &
atom_type_name ) const {
	std::map< std::string, int >::const_iterator
		iter( atom_type_index_.find( atom_type_name ) );
	if ( iter == atom_type_index_.end() ) {
		utility_exit_with_message("unrecognized csd_atom_type_name "+atom_type_name);
	}
	return iter->second;
}

CSDAtomType const &
CSDAtomTypeSet::operator[] ( Size const index ) const
{
	return *( atoms_[ index ] );
}

} // chemical
} // core
