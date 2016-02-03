// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.cc
/// @brief Each line of a PDB file is a Record which is divided into Fields
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


#include <core/io/pdb/Field.hh>

// Utility headers
#include <utility/tools/make_map.hh>


namespace core {
namespace io {
namespace pdb {

using std::string;
using std::ostream;

// Various constructors, for convenience.
Field::Field() :
	value( "" ),
	start( 0 ),
	end( 0 )
{}

Field::Field( core::uint start_in, core::uint end_in, PDBDataType data_type_in ) :
	value( "" ),
	data_type( data_type_in ),
	start( start_in ),
	end( end_in )
{}


// Thank you, Frank diMaio! ~Labonte
void
Field::set_value_from_string( string source ) {
	if ( end != 0 ) {
		value = string( source.begin() + start - 1, source.begin() + end );
	} else {
		value = string( source.begin() + start - 1, source.end() );
	}
}


// Get the PDBDataType value from the corresponding string.
PDBDataType
get_pdb_data_type_from_string( std::string const & type )
{
	if ( type == "AChar" ) {
		return AChar;
	} else if ( type == "Atom" ) {
		return Atom;
	} else if ( type == "Character" ) {
		return Character;
	} else if ( type == "Continuation" ) {
		return Continuation;
	} else if ( type == "Date" ) {
		return Date;
	} else if ( type == "IDcode" ) {
		return IDcode;
	} else if ( type == "Integer" ) {
		return Integer;
	} else if ( type == "List" ) {
		return List;
	} else if ( type == "LString" ) {
		return LString;
	} else if ( type == "Real_type" ) {
		return Real_type;
	} else if ( type == "Record_name" ) {
		return Record_name;
	} else if ( type == "Residue_name" ) {
		return Residue_name;
	} else if ( type == "SList" ) {
		return SList;
	} else if ( type == "Specification" ) {
		return Specification;
	} else if ( type == "Specification_List" ) {
		return Specification_List;
	} else if ( type == "SymOP" ) {
		return SymOP;
	} else {  // String_type is assumed if not recognized, because that is safest.
		return String_type;
	}
}

// Get the string from the corresponding PDBDataType value.
std::string
get_pdb_data_type_from_string( PDBDataType type )
{
	switch ( type ) {
	case AChar : return "AChar";
	case Atom : return "Atom";
	case Character : return "Character";
	case Continuation : return "Continuation";
	case Date : return "Date";
	case IDcode : return "IDcode";
	case Integer : return "Integer";
	case List : return "List";
	case LString : return "LString";
	case Real_type : return "Real_type";
	case Record_name : return "Record_name";
	case SList : return "SList";
	case Specification : return "Specification";
	case Specification_List : return "Specification_List";
	case SymOP : return "SymOP";
	default : return "String_type";  // String_type is assumed if not recognized, because that is safest.
	}
}


std::ostream &
operator<<( std::ostream & os, Field const & F )
{
	os << "[" << F.start << ", " <<  F.end << "]=" << F.value << "";
	return os;
}

std::ostream &
operator<<( std::ostream & os, Record const & R )
{
	for ( Record::const_iterator p = R.begin(), end = R.end(); p != end; ++p ) {
		os << "<Record>{" << p->first << ":" << p->second << "}" << std::endl;
	}
	return os;
}

}  // namespace pdb
}  // namespace io
}  // namespace core
