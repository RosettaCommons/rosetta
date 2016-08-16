// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathVector.srlz.cc
/// @brief  Serlialization routines for MathVectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <numeric/MathVector.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace numeric {

template < class Archive > void save( Archive & archive, numeric::MathVector< float  > const & vect )
{
  numeric::save_math_vector( archive, vect );
}

template < class Archive > void save( Archive & archive, numeric::MathVector< double > const & vect )
{
  numeric::save_math_vector( archive, vect );
}


template < class Archive > void load( Archive & archive, numeric::MathVector< float  > & vect )
{
  numeric::load_math_vector( archive, vect );
}

template < class Archive > void load( Archive & archive, numeric::MathVector< double > & vect )
{
  numeric::load_math_vector( archive, vect );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::MathVector< float  > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::MathVector< double > );

}


#endif
