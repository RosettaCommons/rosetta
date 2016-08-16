// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/serialization/ObjexxFCL/FArray3D.srlz.hh
/// @brief  Serlialization routines for ObjexxFCL::FArray3Ds; these functions have to live in namespace
///         ObjexxFCL, even though they are being written and compiled in the utility library.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_serialization_ObjexxFCL_FArray3D_srlz_HH
#define INCLUDED_utility_serialization_ObjexxFCL_FArray3D_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>

// Package headers
#include <utility/serialization/ObjexxFCL/DynamicIndexRange.srlz.hh>

// Project headers
#include <platform/types.hh>

// C++ headers
#include <string>

namespace ObjexxFCL {

/// @brief Serialization routine for an FArray3D< T >
///
/// @details The serialization support for FArrays is somewhat limited: if the FArray you are working with
/// uses the sophisticated "DimensionExpression" functionality, then this code will not successfully
/// serialize your FArray.  It instead uses the more primitive "IndexRange" (as opposed to DynamicIndexRange)
/// class, which represents the range the FArray covers with a lower and upper integer.
template < class Archive, class T >
void save_farray3d( Archive & archive, FArray3D< T > const & farray3d )
{
	DynamicIndexRange i1( farray3d.I1() ), i2( farray3d.I2() ), i3( farray3d.I3() );
  archive( i1, i2, i3 );
  for ( platform::Size ii = 0; ii < farray3d.size(); ++ii ) {
    archive( farray3d[ii] );
  }
}

/// @brief Deserialization routine for an FArray3D< T >
///
/// @details The serialization support for FArrays is somewhat limited: if the FArray you are working with
/// uses the sophisticated "DimensionExpression" functionality, then this code will not successfully
/// serialize your FArray.  It instead uses the more primitive "IndexRange" (as opposed to DynamicIndexRange)
/// class, which represents the range the FArray covers with a lower and upper integer.
template < class Archive, class T >
void load_farray3d( Archive & archive, FArray3D< T > & farray3d )
{
	DynamicIndexRange ir1, ir2, ir3; archive( ir1, ir2, ir3 );
	farray3d.dimension( ir1, ir2, ir3 );

  for ( platform::Size ii = 0; ii < farray3d.size(); ++ii ) {
    archive( farray3d[ii] );
  }
}

/// @brief Serialization function for arbitary FArray3D classes.
///
/// @details You may serialize FArray3D's holding arbitrary data using this serialization
/// function, or you may define a specialization of this template for a particular
/// datatype that's commonly serialized and put that specialization in a .cc file (to
/// cut down on the compilers effort to instantiate the serialization templates, and on
/// disk space to represent the object files that result from repeated instantiation of
/// the serialization templates.  To be clear: the template specialization will look
/// almost exactly like this function below -- it will invoke "save_farray3d" -- the
/// difference will simply be the replacement of "T" with your type.  After you declare
/// the specialization in your header file (and this specialization must live in namespace
/// ObjexxFCL, even if the thing your farray3d holds lives in another namespace), you will
/// implement the function your .cc file, and below that implementation, you'll use the
/// EXTERNAL_SAVE_AND_LOAD( ObjexxFCL::FArray3D< yourclass > ); macro to force the template
/// instantiation for the set of active archives.  Look for FArray3D.srlz.cc for an example
/// of how these template specialization functions and macros should be written.
template < class Archive, class T >
void save( Archive & archive, typename ObjexxFCL::FArray3D< T > const & farray3d )
{
  save_farray3d( archive, farray3d );
}

/// @brief Deserialization function for FArray3D's holding arbitrary data
template < class Archive, class T >
void load( Archive & archive, typename ObjexxFCL::FArray3D< T > & farray3d )
{
  load_farray3d( archive, farray3d );
}


// template specialization for a handful of commonly used FArray3Ds; the implementation for these
// methods will live in FArray3D.srlz.cc so that they don't have to be regenerated in as many
// object files.

/// @brief partial template specialization for FArray3Ds of integers; the purpose of this specialization
/// is solely to reduce compile time and disk space, since FArray3Ds of integers are fairly common
template < class Archive > void save( Archive & archive, FArray3D< std::string    > const & );
template < class Archive > void save( Archive & archive, FArray3D< int            > const & );
template < class Archive > void save( Archive & archive, FArray3D< platform::Size > const & );
template < class Archive > void save( Archive & archive, FArray3D< float          > const & );
template < class Archive > void save( Archive & archive, FArray3D< double         > const & );

template < class Archive > void load( Archive & archive, FArray3D< std::string    > & );
template < class Archive > void load( Archive & archive, FArray3D< int            > & );
template < class Archive > void load( Archive & archive, FArray3D< platform::Size > & );
template < class Archive > void load( Archive & archive, FArray3D< float          > & );
template < class Archive > void load( Archive & archive, FArray3D< double         > & );

}

#endif // SERIALIZATION

#endif // INCLUDED_utility_serialization_ObjexxFCL_FArray3D_srlz_HH
