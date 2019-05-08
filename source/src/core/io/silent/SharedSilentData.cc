// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/SharedSilentData.cc
/// @brief  Serialization routines for shared silent data.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/io/silent/SharedSilentData.hh>


#ifdef SERIALIZATION
#include <utility/serialization/serialization.hh>
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>


template< class Archive >
void core::io::silent::SharedSilentData::save( Archive & ) const {}

template< class Archive >
void core::io::silent::SharedSilentData::load( Archive & ) {}


template< class Archive >
void core::io::silent::SimpleSequenceData::save( Archive & arc ) const
{
	arc( cereal::base_class<core::io::silent::SharedSilentData>( this ) );
	arc( sequence_ );
}

template< class Archive >
void core::io::silent::SimpleSequenceData::load( Archive & arc )
{
	arc( cereal::base_class<core::io::silent::SharedSilentData>( this ) );
	arc( sequence_ );
}


SAVE_AND_LOAD_SERIALIZABLE( core::io::silent::SharedSilentData );
SAVE_AND_LOAD_SERIALIZABLE( core::io::silent::SimpleSequenceData );
CEREAL_REGISTER_TYPE( core::io::silent::SharedSilentData )
CEREAL_REGISTER_TYPE( core::io::silent::SimpleSequenceData )

CEREAL_REGISTER_DYNAMIC_INIT( core_io_silent_SharedSilentData )


#endif



