// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key3Tuple.srlz.hh
/// @brief  Serlialization routines for Key3Tuples
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_keys_Key3Tuple_srlz_HH
#define INCLUDED_utility_keys_Key3Tuple_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <utility/keys/Key3Tuple.hh>

// Project headers
#include <platform/types.hh>

// C++ headers
#include <string>

namespace utility {
namespace keys {

template < class Archive, class R, class S, class T >
void save_key3tuple( Archive & arc, Key3Tuple< R, S, T > const & key3tup )
{
	arc( key3tup.key1(), key3tup.key2(), key3tup.key3() );
}

template < class Archive, class R, class S, class T >
void load_key3tuple( Archive & arc, Key3Tuple< R, S, T > & key3tup )
{
	arc( key3tup.key1(), key3tup.key2(), key3tup.key3() );
}

/// @brief Serialization function for arbitary Key3Tuple classes.
///
/// @details You may serialize Key3Tuple's holding arbitrary data using this serialization
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
/// EXTERNAL_SAVE_AND_LOAD( Key3Tuple< yourclass > ); macro to force the template
/// instantiation for the set of active archives.  Look for Key3Tuple.srlz.cc for an example
/// of how these template specialization functions and macros should be written.
template < class Archive, class R, class S, class T >
void save( Archive & arc, Key3Tuple< R, S, T > const & k3t )
{
  save_key3tuple( arc, k3t );
}

/// @brief Deserialization function for Key3Tuple's holding arbitrary data
template < class Archive, class R, class S, class T >
void load( Archive & arc, Key3Tuple< R, S, T > & k3t )
{
  load_key3tuple( arc, k3t );
}


/// @brief partial template specialization for Key3Tuples of integers; the purpose of this specialization
/// is solely to reduce compile time and disk space
template < class Archive > void save( Archive & arc, Key3Tuple< platform::Size, platform::Size, platform::Size > const & );
template < class Archive > void load( Archive & arc, Key3Tuple< platform::Size, platform::Size, platform::Size > & );

}
}

#endif // SERIALIZATION

#endif // INCLUDED_utility_Key3Tuple_srlz_HH
