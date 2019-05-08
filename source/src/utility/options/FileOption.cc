// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   utility/options/FileOption.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <utility/options/FileOption.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>

/// @brief Automatically generated serialization method
template< class Archive >
void
utility::options::FileOption::save( Archive & arc ) const {
	arc( cereal::base_class< utility::options::ScalarOption_T_<class utility::options::FileOptionKey, class utility::file::FileName> >( this ) );
	arc( CEREAL_NVP( cached_value_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
utility::options::FileOption::load( Archive & arc ) {
	arc( cereal::base_class< utility::options::ScalarOption_T_<class utility::options::FileOptionKey, class utility::file::FileName> >( this ) );
	arc( cached_value_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( utility::options::FileOption );
CEREAL_REGISTER_TYPE( utility::options::FileOption )

CEREAL_REGISTER_DYNAMIC_INIT( utility_options_FileOption )

#endif // SERIALIZATION
