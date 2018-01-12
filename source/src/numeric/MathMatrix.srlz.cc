// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathMatrix.srlz.cc
/// @brief  Serlialization routines for MathMatrix
/// @author Rebecca Alford (ralford3@jhu.edu)

#ifdef SERIALIZATION

// Unit headers
#include <numeric/MathMatrix.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace numeric {

template < class Archive >
void
save( Archive & archive, numeric::MathMatrix< float  > const & matrix )
{
	numeric::save_math_matrix( archive, matrix );
}

template < class Archive >
void
save( Archive & archive, numeric::MathMatrix< double > const & matrix )
{
	numeric::save_math_matrix( archive, matrix );
}


template < class Archive >
void
load( Archive & archive, numeric::MathMatrix< float  > & matrix )
{
	numeric::load_math_matrix( archive, matrix );
}

template < class Archive >
void
load( Archive & archive, numeric::MathMatrix< double > & matrix )
{
	numeric::load_math_matrix( archive, matrix );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::MathMatrix< float  > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( numeric::MathMatrix< double > );

}


#endif
