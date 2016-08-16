// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/fixedsizearray1.srlz.hh
/// @brief  Serlialization routines for fixedsizearray1s
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_fixedsizearray1_srlz_HH
#define INCLUDED_utility_fixedsizearray1_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <utility/fixedsizearray1.hh>

// C++ headers
#include <string>

namespace utility {

template < class Archive, class T, platform::Size S >
void save_fixedsizearray1( Archive & archive, utility::fixedsizearray1< T, S > const & vect )
{
  for ( platform::Size ii = 1; ii <= vect.size(); ++ii ) {
    archive( vect[ii] );
  }
}

template < class Archive, class T, platform::Size S >
void load_fixedsizearray1( Archive & archive, utility::fixedsizearray1< T, S > & vect )
{
  for ( platform::Size ii = 1; ii <= vect.size(); ++ii ) {
    archive( vect[ii] );
  }
}

/// @brief Serialization function for arbitary fixedsizearray1 classes.
///
/// @details You may serialize fixedsizearray1's holding arbitrary data using this serialization
/// function, or you may define a specialization of this template for a particular
/// datatype that's commonly serialized and put that specialization in a .cc file (to
/// cut down on the compilers effort to instantiate the serialization templates, and on
/// disk space to represent the object files that result from repeated instantiation of
/// the serialization templates.  To be clear: the template specialization will look
/// almost exactly like this function below -- it will invoke "save_fixedsizearray1" -- the
/// difference will simply be the replacement of "T" with your type.  After you declare
/// the specialization in your header file (and this specialization must live in namespace
/// utility, even if the thing your fixedsizearray1 holds lives in another namespace), you will
/// implement the function your .cc file, and below that implementation, you'll use the
/// EXTERNAL_SAVE_AND_LOAD( utility::fixedsizearray1< yourclass > ); macro to force the template
/// instantiation for the set of active archives.  Look for fixedsizearray1.srlz.cc for an example
/// of how these template specialization functions and macros should be written.
template < class Archive, class T, platform::Size S >
void save( Archive & archive, typename utility::fixedsizearray1< T, S > const & vect )
{
  save_fixedsizearray1( archive, vect );
}

/// @brief Deserialization function for fixedsizearray1's holding arbitrary data
template < class Archive, class T, platform::Size S >
void load( Archive & archive, utility::fixedsizearray1< T, S > & vect )
{
  load_fixedsizearray1( archive, vect );
}



}

#endif // SERIALIZATION

#endif // INCLUDED_utility_fixedsizearray1_srlz_HH
