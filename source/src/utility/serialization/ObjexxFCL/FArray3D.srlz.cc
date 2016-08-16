// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/serialization/FArray3D.srlz.cc
/// @brief  Serlialization routines for FArray3Ds
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/serialization/ObjexxFCL/FArray3D.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace ObjexxFCL {

template < class Archive > void save( Archive & archive, FArray3D< std::string    > const & farray3d )
{
  save_farray3d( archive, farray3d );
}

template < class Archive > void save( Archive & archive, FArray3D< int            > const & farray3d )
{
  save_farray3d( archive, farray3d );
}

template < class Archive > void save( Archive & archive, FArray3D< platform::Size > const & farray3d )
{
  save_farray3d( archive, farray3d );
}

template < class Archive > void save( Archive & archive, FArray3D< float          > const & farray3d )
{
  save_farray3d( archive, farray3d );
}

template < class Archive > void save( Archive & archive, FArray3D< double         > const & farray3d )
{
  save_farray3d( archive, farray3d );
}

template < class Archive > void load( Archive & archive, FArray3D< std::string    > & farray3d )
{
  load_farray3d( archive, farray3d );
}

template < class Archive > void load( Archive & archive, FArray3D< int            > & farray3d )
{
  load_farray3d( archive, farray3d );
}

template < class Archive > void load( Archive & archive, FArray3D< platform::Size > & farray3d )
{
  load_farray3d( archive, farray3d );
}

template < class Archive > void load( Archive & archive, FArray3D< float          > & farray3d )
{
  load_farray3d( archive, farray3d );
}

template < class Archive > void load( Archive & archive, FArray3D< double         > & farray3d )
{
  load_farray3d( archive, farray3d );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray3D< std::string    > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray3D< int            > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray3D< platform::Size > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray3D< float          > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray3D< double         > );

}


#endif
