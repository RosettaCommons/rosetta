// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vector0.srlz.hh
/// @brief  Serlialization routines for vector0s
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_vector0_srlz_HH
#define INCLUDED_utility_vector0_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <utility/vector0.hh>

// C++ headers
#include <string>

namespace utility {

template < class Archive, class T >
void save_vector( Archive & archive, utility::vector0< T > const & vect )
{
  archive( vect.size() );
  for ( platform::Size ii = 0; ii < vect.size(); ++ii ) {
    archive( vect[ii] );
  }
}

template < class Archive, class T >
void load_vector( Archive & archive, utility::vector0< T > & vect )
{
	platform::Size n; archive( n );
	vect.resize( n );
  for ( platform::Size ii = 0; ii < vect.size(); ++ii ) {
    archive( vect[ii] );
  }
}

/// @brief Serialization function for arbitary vector0 classes.
///
/// @details You may serialize vector0's holding arbitrary data using this serialization
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
/// EXTERNAL_SAVE_AND_LOAD( utility::vector0< yourclass > ); macro to force the template
/// instantiation for the set of active archives.  Look for vector0.srlz.cc for an example
/// of how these template specialization functions and macros should be written.
template < class Archive, class T >
void save( Archive & archive, typename utility::vector0< T > const & vect )
{
  save_vector( archive, vect );
}

/// @brief Deserialization function for vector0's holding arbitrary data
template < class Archive, class T >
void load( Archive & archive, utility::vector0< T > & vect )
{
  load_vector( archive, vect );
}


// template specialization for a handful of commonly used vector0s; the implementation for these
// methods will live in vector0.srlz.cc so that they don't have to be regenerated in as many
// object files.

/// @brief partial template specialization for vector0s of integers; the purpose of this specialization
/// is solely to reduce compile time and disk space, since vector0s of integers are fairly common
template < class Archive > void save( Archive & archive, utility::vector0< std::string    > const & );
template < class Archive > void save( Archive & archive, utility::vector0< int            > const & );
template < class Archive > void save( Archive & archive, utility::vector0< platform::Size > const & );
template < class Archive > void save( Archive & archive, utility::vector0< float          > const & );
template < class Archive > void save( Archive & archive, utility::vector0< double         > const & );
template < class Archive > void save( Archive & archive, utility::vector0< bool           > const & );

template < class Archive > void load( Archive & archive, utility::vector0< std::string    > & );
template < class Archive > void load( Archive & archive, utility::vector0< int            > & );
template < class Archive > void load( Archive & archive, utility::vector0< platform::Size > & );
template < class Archive > void load( Archive & archive, utility::vector0< float          > & );
template < class Archive > void load( Archive & archive, utility::vector0< double         > & );
template < class Archive > void load( Archive & archive, utility::vector0< bool           > & );

}

#endif // SERIALIZATION

#endif // INCLUDED_utility_vector0_srlz_HH
