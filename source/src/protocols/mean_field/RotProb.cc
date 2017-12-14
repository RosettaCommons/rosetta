// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RotProb.cc

/// @brief   Method definitions for RotProb.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/RotProb.hh>

// Package headers


// Project headers
#include <core/conformation/Residue.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers


// C++ headers
#include <iostream>
#include <iomanip>

// Boost headers
#include <boost/format.hpp>
#include <utility>

// Construct tracer.
static basic::Tracer TR("protocols.mean_field.RotProb");


namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details initialized to nonsense values
RotProb::RotProb() :
	probability_( Real( 0.0 ) ),
	rot_ind_( Size( 0 ) ),
	pos_( Size ( 0 ) ),
	res_( /* 0 */ )
{}

/// @details main constructor
/// @param [in] prob - probability of Rotamer occurring at the site
/// @param [in] rot_ind - used to uniquely identify rotamer across backbones
/// @param [in] pos - position of the site in pose numbering
/// @param [in] res - COP to current residue at site (useful for AA identity)
RotProb::RotProb( Real prob, Size rot_ind, Size pos, core::conformation::ResidueCOP res ) :
	probability_( prob ),
	rot_ind_( rot_ind ),
	pos_( pos ),
	res_(std::move( res ))
{}

/// @details Copy constructor uses private member function
RotProb::RotProb( RotProb const & object_to_copy )
{
	copy_data( *this, object_to_copy );
}

/// @details Assignment operator uses private member function
RotProb &
RotProb::operator=( RotProb const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data( *this, object_to_copy );
	return *this;
}

/// @brief Destructor
RotProb::~RotProb() = default;


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

/// @details prints RotProb name, pos, and probability. ex: LEU132 0.5
void
RotProb::show( std::ostream & output ) const
{
	std::string name = core::chemical::name_from_aa( aa_ind( ) );
	output << boost::format( "%1%%2%\t%3%" ) % name % pos_ % probability_;
}


// Accessors/Mutators
core::conformation::ResidueCOP
RotProb::res() const
{
	return res_;
}

void
RotProb::res( core::conformation::ResidueCOP r )
{
	res_ = r;
}

/// @details uses ResidueCOP to access AA identity of RotProb
core::chemical::AA
RotProb::aa_ind() const
{
	return res()->aa();
}



// Private methods /////////////////////////////////////////////////////////////
// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RotProb::copy_data(
	RotProb & object_to_copy_to, RotProb const & object_to_copy_from )
{
	object_to_copy_to.probability( object_to_copy_from.probability() );
	object_to_copy_to.rot_ind( object_to_copy_from.rot_ind() );
	object_to_copy_to.pos( object_to_copy_from.pos() );
	object_to_copy_to.res( object_to_copy_from.res() );

}


// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that AAProb can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, RotProb const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace mean_field
}  // namespace protocols
