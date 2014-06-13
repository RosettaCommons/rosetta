// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/mm/mmtrie/MMEnergyTableAtom.cc
/// @brief  Implimentation for the MMEnergyTableAtom. Heavily coppied from the EtableAtom.cc
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit Headers
#include <core/scoring/mm/mmtrie/MMEnergyTableAtom.hh>

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
namespace mm {
namespace mmtrie {


MMEnergyTableAtom::MMEnergyTableAtom() : parent(), is_hydrogen_( false ) {}

MMEnergyTableAtom::MMEnergyTableAtom( conformation::Residue const & res, Size atom_index )
:
	parent( res.atom( atom_index ) ),
	is_hydrogen_( false )
{}

MMEnergyTableAtom::~MMEnergyTableAtom() {}

/// @brief send a description of the atom to standard out
void
MMEnergyTableAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
MMEnergyTableAtom::print( std::ostream & os ) const
{
	os << "mm atom type" << mm_type() << " ";
	os << "(" << xyz().x();
	os << ", " << xyz().y();
	os << ", " << xyz().z() << ")" << std::endl;
}

std::ostream & operator << ( std::ostream & os, MMEnergyTableAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace mmtrie
} // namespace mm
} // namespace scoring
} // namespace core

