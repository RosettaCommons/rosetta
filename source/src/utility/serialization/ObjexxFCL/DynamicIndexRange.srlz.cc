// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/serialization/DynamicIndexRange.srlz.cc
/// @brief  Serlialization routines for DynamicIndexRanges
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/serialization/ObjexxFCL/DynamicIndexRange.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace ObjexxFCL {

template < class Archive >
void save( Archive & arc, DynamicIndexRange const & ir )
{
	// EXEMPT l_dim_p_ u_dim_p_
	arc( ir.l(), ir.u() );
}

template < class Archive >
void load( Archive & arc, DynamicIndexRange & ir )
{
	// EXEMPT l_dim_p_ u_dim_p_
	int l, u;
	arc( l, u );
	ir.assign( l, u );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( DynamicIndexRange );

}


#endif
