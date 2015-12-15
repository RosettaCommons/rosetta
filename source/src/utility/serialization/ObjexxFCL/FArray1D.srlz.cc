// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/serialization/ObjexxFCL/FArray1D.srlz.cc
/// @brief  Serlialization routines for FArray1Ds
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/serialization/ObjexxFCL/FArray1D.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace ObjexxFCL {

template < class Archive > void save( Archive & archive, FArray1D< std::string    > const & farray1d )
{
  save_farray1d( archive, farray1d );
}

template < class Archive > void save( Archive & archive, FArray1D< int            > const & farray1d )
{
  save_farray1d( archive, farray1d );
}

template < class Archive > void save( Archive & archive, FArray1D< platform::Size > const & farray1d )
{
  save_farray1d( archive, farray1d );
}

template < class Archive > void save( Archive & archive, FArray1D< float          > const & farray1d )
{
  save_farray1d( archive, farray1d );
}

template < class Archive > void save( Archive & archive, FArray1D< double         > const & farray1d )
{
  save_farray1d( archive, farray1d );
}

template < class Archive > void load( Archive & archive, FArray1D< std::string    > & farray1d )
{
  load_farray1d( archive, farray1d );
}

template < class Archive > void load( Archive & archive, FArray1D< int            > & farray1d )
{
  load_farray1d( archive, farray1d );
}

template < class Archive > void load( Archive & archive, FArray1D< platform::Size > & farray1d )
{
  load_farray1d( archive, farray1d );
}

template < class Archive > void load( Archive & archive, FArray1D< float          > & farray1d )
{
  load_farray1d( archive, farray1d );
}

template < class Archive > void load( Archive & archive, FArray1D< double         > & farray1d )
{
  load_farray1d( archive, farray1d );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray1D< std::string    > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray1D< int            > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray1D< platform::Size > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray1D< float          > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( FArray1D< double         > );

}


#endif
