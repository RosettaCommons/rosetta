// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/Key2Tuple.srlz.cc
/// @brief  Serlialization routines for Key2Tuples
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/keys/Key2Tuple.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Boost headers
#include <boost/preprocessor/punctuation/comma.hpp>

namespace utility {
namespace keys {

template < class Archive > void save( Archive & arc, Key2Tuple< platform::Size, platform::Size > const & k2t )
{ save_key2tuple( arc, k2t ); }

template < class Archive > void load( Archive & arc, Key2Tuple< platform::Size, platform::Size > & k2t )
{ load_key2tuple( arc, k2t ); }

template < class Archive > void save( Archive & arc, Key2Tuple< double, double > const & k2t )
{ save_key2tuple( arc, k2t ); }
template < class Archive > void load( Archive & arc, Key2Tuple< double, double > & k2t )
{ load_key2tuple( arc, k2t ); }

template < class Archive > void save( Archive & arc, Key2Tuple< float, float > const & k2t )
{ save_key2tuple( arc, k2t ); }
template < class Archive > void load( Archive & arc, Key2Tuple< float, float > & k2t )
{ load_key2tuple( arc, k2t ); }


EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( Key2Tuple< platform::Size BOOST_PP_COMMA() platform::Size > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( Key2Tuple< float BOOST_PP_COMMA() float > );
EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( Key2Tuple< double BOOST_PP_COMMA() double > );

}
}

#endif
