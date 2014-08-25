// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Atom.cc
/// @brief  Method definitions for conformation::Atom
/// @note   not to be confused with chemical::Atom
/// @author Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/conformation/Atom.hh>


namespace core {
namespace conformation {

void
Atom::show( std::ostream & output ) const
{
	output << xyz_.x() << ", " << xyz_.y() << ", " << xyz_.z();
}


std::ostream &
operator << ( std::ostream & out, Atom const & atom )
{
	atom.show( out );
	return out;
}

}  // conformation
}  // core
