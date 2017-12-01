// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathVector.srlz.hh
/// @brief  Serlialization routines for MathVector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_numeric_MatVector_srlz_HH
#define INCLUDED_numeric_MatVector_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <numeric/MathVector.hh>

namespace numeric {

template < class Archive, class T >
void save_math_vector( Archive & archive, numeric::MathVector< T > const & vect )
{
	archive( vect.size() );
	for ( platform::Size ii = 0; ii < vect.size(); ++ii ) {
		archive( vect(ii) );
	}
}

template < class Archive, class T >
void load_math_vector( Archive & archive, numeric::MathVector< T > & vect )
{
	platform::Size n; archive( n );
	numeric::MathVector< T > second_vect( n );
	for ( platform::Size ii = 0; ii < second_vect.size(); ++ii ) {
		archive( second_vect(ii) );
	}
	vect = second_vect;
}

/// @brief Serialization function for arbitary MathVector classes.
///
/// @details You may serialize MathVector's holding arbitrary data using this serialization
/// function, or you may define a specialization of this template for a particular
/// datatype that's commonly serialized and put that specialization in a .cc file (to
/// cut down on the compilers effort to instantiate the serialization templates, and on
/// disk space to represent the object files that result from repeated instantiation of
/// the serialization templates.  To be clear: the template specialization will look
/// almost exactly like this function below -- it will invoke "save_math_vector" -- the
/// difference will simply be the replacement of "T" with your type.  After you declare
/// the specialization in your header file (and this specialization must live in namespace
/// utility, even if the thing your vector holds lives in another namespace), you will
/// implement the function your .cc file, and below that implementation, you'll use the
/// EXTERNAL_SAVE_AND_LOAD( numeric::MathVector< yourclass > ); macro to force the template
/// instantiation for the set of active archives.  Look for MathVector.srlz.cc for an example
/// of how these template specialization functions and macros should be written.
template < class Archive, class T >
void save( Archive & archive, typename numeric::MathVector< T > const & vect )
{
	save_math_vector( archive, vect );
}

/// @brief Deserialization function for MathVector's holding arbitrary data
template < class Archive, class T >
void load( Archive & archive, numeric::MathVector< T > & vect )
{
	load_math_vector( archive, vect );
}


// template specialization for a handful of commonly used MathVectors; the implementation for these
// methods will live in MathVector.srlz.cc so that they don't have to be regenerated in as many
// object files.

/// @brief partial template specialization for MathVectors of integers; the purpose of this specialization
/// is solely to reduce compile time and disk space, since MathVectors of integers are fairly common
template < class Archive > void save( Archive & archive, numeric::MathVector< float          > const & );
template < class Archive > void save( Archive & archive, numeric::MathVector< double         > const & );

template < class Archive > void load( Archive & archive, numeric::MathVector< float          > & );
template < class Archive > void load( Archive & archive, numeric::MathVector< double         > & );

}

#endif // SERIALIZATION

#endif // INCLUDED_numeric_MathVector_srlz_HH
