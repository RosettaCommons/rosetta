// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/hbtrie/LKBAtom.hh
/// @brief
/// @author

// Unit Headers
#include <core/scoring/lkball/lkbtrie/LKBAtom.hh>

// Project Headers
#include <core/types.hh>

// STL Headers
#include <iostream>

// Numceric Headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace lkball {
namespace lkbtrie {


LKBAtom::LKBAtom() {}

LKBAtom::~LKBAtom() {}

/// @brief send a description of the atom to standard out
void
LKBAtom::print() const { print( std::cout ); }

/// @brief send a description of the atom to an output stream
void
LKBAtom::print( std::ostream & os ) const
{
	os << "LKBAtom" <<  " ";
	os << "(" << base_.xyz().x();
	os << ", " << base_.xyz().y();
	os << ", " << base_.xyz().z() << ")";
}

std::ostream & operator << ( std::ostream & os, LKBAtom const & atom )
{
	atom.print( os );
	return os;
}

} // namespace lkbtrie
} // namespace lkball
} // namespace scoring
} // namespace core

