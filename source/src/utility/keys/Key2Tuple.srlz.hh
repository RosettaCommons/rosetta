// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key2Tuple.srlz.hh
/// @brief  Serlialization routines for Key2Tuples
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_keys_Key2Tuple_srlz_HH
#define INCLUDED_utility_keys_Key2Tuple_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <utility/keys/Key2Tuple.hh>

// Project headers
#include <platform/types.hh>

// C++ headers
#include <string>

namespace utility {
namespace keys {

template < class Archive, class S, class T >
void save_key2tuple( Archive & arc, Key2Tuple< S, T > const & key2tup )
{
	arc( key2tup.key1(), key2tup.key2() );
}

template < class Archive, class S, class T >
void load_key2tuple( Archive & arc, Key2Tuple< S, T > & key2tup )
{
	arc( key2tup.key1(), key2tup.key2() );
}

/// @brief Serialization function for arbitary Key2Tuple classes.
///
/// @details You may serialize Key2Tuple's holding arbitrary data using this serialization
/// function, or you may define a specialization of this template for a particular
/// datatype that's commonly serialized and put that specialization in a .cc file (to
/// cut down on the compilers effort to instantiate the serialization templates, and on
/// disk space to represent the object files that result from repeated instantiation of
/// the serialization templates.  To be clear: the template specialization will look
/// almost exactly like this function below -- it will invoke "save_vector" -- the
/// difference will simply be the replacement of "T" with your type.  After you declare
/// the specialization in your header file (and this specialization must live in namespace
/// utility, even if the thing your vector holds lives in another namespace), you will
/// implement the function your .cc file, and below that implementation, you'll use the
/// EXTERNAL_SAVE_AND_LOAD( Key2Tuple< yourclass > ); macro to force the template
/// instantiation for the set of active archives.  Look for Key2Tuple.srlz.cc for an example
/// of how these template specialization functions and macros should be written.
template < class Archive, class S, class T >
void save( Archive & arc, Key2Tuple< S, T > const & k2t )
{
  save_key2tuple( arc, k2t );
}

/// @brief Deserialization function for Key2Tuple's holding arbitrary data
template < class Archive, class S, class T >
void load( Archive & arc, Key2Tuple< S, T > & k2t )
{
  load_key2tuple( arc, k2t );
}


/// @brief partial template specialization for Key2Tuples of integers; the purpose of this specialization
/// is solely to reduce compile time and disk space
template < class Archive > void save( Archive & arc, Key2Tuple< platform::Size, platform::Size > const & );
template < class Archive > void load( Archive & arc, Key2Tuple< platform::Size, platform::Size > & );

template < class Archive > void save( Archive & arc, Key2Tuple< double, double > const & );
template < class Archive > void load( Archive & arc, Key2Tuple< double, double > & );

template < class Archive > void save( Archive & arc, Key2Tuple< float, float > const & );
template < class Archive > void load( Archive & arc, Key2Tuple< float, float > & );


}
}

#endif // SERIALIZATION

#endif // INCLUDED_utility_Key2Tuple_srlz_HH
