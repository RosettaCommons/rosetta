// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/serialization/ubyte.srlz.cc
/// @brief  Serlialization routines for ubytes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/serialization/ObjexxFCL/ubyte.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace ObjexxFCL {

template < class Archive >
void save( Archive & arc, ubyte const & ub )
{
	unsigned short b( ub );
	arc( b );
	// The b_ data member cannot be directly accessed
	// EXEMPT b_
}

template < class Archive >
void load( Archive & arc, ubyte & ir )
{
	unsigned short b;
	arc( b );
	ir = b;
	// The b_ data member cannot be directly accessed
	// EXEMPT b_
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( ubyte );

}


#endif
