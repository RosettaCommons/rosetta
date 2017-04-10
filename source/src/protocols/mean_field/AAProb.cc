// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    AAProb.cc

/// @brief   Method definitions for AAProb.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/AAProb.hh>

// Package headers
#include <protocols/mean_field/RotProb.hh>

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

// Construct tracer.
static THREAD_LOCAL basic::Tracer TR("protocols.mean_field.AAProb");


namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details initialized to nonsense values
AAProb::AAProb() :
	probability_( Real( 0.0 ) ),
	aa_ind_( core::chemical::num_aa_types ),
	pos_( Size( 0 ) ),
	nrot_( Size( 0 ) )
{}

/// @param [in] prob - probability of Amino Acid occurring at the site
/// @param [in] aa_ind - AA enum of identity
/// @param [in] pos - position of the site in pose numbering
/// @param [in] nrot - number of rotamers corresponding to this amino acid at that site
AAProb::AAProb( Real prob, core::chemical::AA aa_ind, Size pos, Size nrot ) :
	probability_( prob ),
	aa_ind_( aa_ind ),
	pos_( pos ),
	nrot_( nrot )
{}

/// @details constructs AAProb from RotProb
AAProb::AAProb( RotProb const & rp ) :
	probability_( rp.probability() ),
	aa_ind_( rp.res()->aa() ),
	pos_( rp.pos() ),
	nrot_( Size( 1 ) )
{}

/// @details Copy constructor uses private member function
AAProb::AAProb( AAProb const & object_to_copy )
{
	copy_data( *this, object_to_copy );
}

/// @details Assignment constructor uses private member function
AAProb &
AAProb::operator=( AAProb const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data( *this, object_to_copy );
	return *this;
}


/// @brief Destructor
AAProb::~AAProb() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

/// @details prints AAProb name, pos, and probability. ex: LEU132 0.5
void
AAProb::show( std::ostream & output ) const
{
	std::string name = core::chemical::name_from_aa( aa_ind() );
	output << boost::format("%1%%2%\t%3%") % name % pos_ % probability_;
}



// Private methods /////////////////////////////////////////////////////////////

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
AAProb::copy_data(
	AAProb & object_to_copy_to,
	AAProb const & object_to_copy_from )
{
	object_to_copy_to.probability( object_to_copy_from.probability() );
	object_to_copy_to.aa_ind( object_to_copy_from.aa_ind() );
	object_to_copy_to.pos( object_to_copy_from.pos() );
	object_to_copy_to.nrot( object_to_copy_from.nrot() );
}


// Friend methods //////////////////////////////////////////////////////////////

// Insertion operator (overloaded so that AAProb can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, AAProb const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace mean_field
}  // namespace protocols
