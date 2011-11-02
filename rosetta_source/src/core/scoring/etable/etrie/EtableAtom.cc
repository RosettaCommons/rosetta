// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/EtableAtom.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/EtableAtom.hh>

// Project Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {


EtableAtom::EtableAtom() : parent(), is_hydrogen_( false ) {}

EtableAtom::EtableAtom( conformation::Residue const & res, Size atom_index )
:
	parent( res.atom( atom_index ) ),
	is_hydrogen_( false )
{}

EtableAtom::~EtableAtom() {}

/// @brief send a description of the atom to standard out
void
EtableAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
EtableAtom::print( std::ostream & os ) const
{
	os << "atom type" << type() << " ";
	os << "(" << xyz().x();
	os << ", " << xyz().y();
	os << ", " << xyz().z() << ")" << std::endl;
}

std::ostream & operator << ( std::ostream & os, EtableAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

