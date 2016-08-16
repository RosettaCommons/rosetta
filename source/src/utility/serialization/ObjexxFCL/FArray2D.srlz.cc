// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/serialization/FArray2D.srlz.cc
/// @brief  Serlialization routines for FArray2Ds
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/serialization/ObjexxFCL/FArray2D.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace ObjexxFCL {

template < class Archive > void save( Archive & archive, FArray2D< std::string    > const & farray2d )
{
  save_farray2d( archive, farray2d );
}

template < class Archive > void save( Archive & archive, FArray2D< int            > const & farray2d )
{
  save_farray2d( archive, farray2d );
}

template < class Archive > void save( Archive & archive, FArray2D< platform::Size > const & farray2d )
{
  save_farray2d( archive, farray2d );
}

template < class Archive > void save( Archive & archive, FArray2D< float          > const & farray2d )
{
  save_farray2d( archive, farray2d );
}

template < class Archive > void save( Archive & archive, FArray2D< double         > const & farray2d )
{
  save_farray2d( archive, farray2d );
}

template < class Archive > void load( Archive & archive, FArray2D< std::string    > & farray2d )
{
  load_farray2d( archive, farray2d );
}

template < class Archive > void load( Archive & archive, FArray2D< int            > & farray2d )
{
  load_farray2d( archive, farray2d );
}

template < class Archive > void load( Archive & archive, FArray2D< platform::Size > & farray2d )
{
  load_farray2d( archive, farray2d );
}

template < class Archive > void load( Archive & archive, FArray2D< float          > & farray2d )
{
  load_farray2d( archive, farray2d );
}

template < class Archive > void load( Archive & archive, FArray2D< double         > & farray2d )
{
  load_farray2d( archive, farray2d );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray2D< std::string    > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray2D< int            > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray2D< platform::Size > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray2D< float          > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray2D< double         > );

}


#endif
