// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector0.srlz.cc
/// @brief  Serlialization routines for vector0s
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/vector0.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>

namespace utility {

template < class Archive > void save( Archive & archive, utility::vector0< std::string    > const & vect )
{
  utility::save_vector( archive, vect );
}

template < class Archive > void save( Archive & archive, utility::vector0< int            > const & vect )
{
  utility::save_vector( archive, vect );
}

template < class Archive > void save( Archive & archive, utility::vector0< platform::Size > const & vect )
{
  utility::save_vector( archive, vect );
}

template < class Archive > void save( Archive & archive, utility::vector0< float          > const & vect )
{
  utility::save_vector( archive, vect );
}

template < class Archive > void save( Archive & archive, utility::vector0< double         > const & vect )
{
  utility::save_vector( archive, vect );
}

template < class Archive > void save( Archive & archive, utility::vector0< bool           > const & vect )
{
  archive( vect.size() );
  for ( platform::Size ii = 0; ii < vect.size(); ++ii ) {
		bool iival = vect[ii];
    archive( iival );
  }
}

template < class Archive > void load( Archive & archive, utility::vector0< std::string    > & vect )
{
  utility::load_vector( archive, vect );
}

template < class Archive > void load( Archive & archive, utility::vector0< int            > & vect )
{
  utility::load_vector( archive, vect );
}

template < class Archive > void load( Archive & archive, utility::vector0< platform::Size > & vect )
{
  utility::load_vector( archive, vect );
}

template < class Archive > void load( Archive & archive, utility::vector0< float          > & vect )
{
  utility::load_vector( archive, vect );
}

template < class Archive > void load( Archive & archive, utility::vector0< double         > & vect )
{
  utility::load_vector( archive, vect );
}

template < class Archive > void load( Archive & archive, utility::vector0< bool           > & vect )
{
	platform::Size n; archive( n );
	vect.resize( n );
  for ( platform::Size ii = 0; ii < vect.size(); ++ii ) {
		bool iival; archive( iival );
    vect[ii] = iival;
  }
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( utility::vector0< int            > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( utility::vector0< platform::Size > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( utility::vector0< float          > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( utility::vector0< double         > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( utility::vector0< bool           > );

}


#endif
