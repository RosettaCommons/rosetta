// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.cc
/// @brief  Each line of a PDB file is a Record which is divided into Fields
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/pdb/Field.hh>


namespace core {
namespace io {
namespace pdb {

// Various constructors, for convenience.
Field::Field() :
		value( "" ),
		start( 0 ),
		end( 0 )
{}

Field::Field( core::uint start_in, core::uint end_in  ) :
		value( "" ),
		start( start_in ),
		end( end_in )
{}


// Read field value from given .pdb line and set.
// Thank you, Frank diMaio for idea to have zero indicate no upper bound! ~Labonte
void
Field::set_value_from_pdb_line( std::string source ) {
	if ( end != 0 ) {
		value = std::string( source.begin() + start - 1, source.begin() + end );
	} else {
		value = std::string( source.begin() + start - 1, source.end() );
	}
}


// Helper Function ////////////////////////////////////////////////////////////
std::ostream &
operator<<( std::ostream & os, Field const & field )
{
	os << "[" << field.start << ", " <<  field.end << "]=" << field.value << "";
	return os;
}

}  // namespace pdb
}  // namespace io
}  // namespace core
