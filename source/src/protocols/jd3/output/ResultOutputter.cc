// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/output/ResultOutputter.cc
/// @brief  Method definitions for the %ResultOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/output/ResultOutputter.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace output {

ResultOutputter::ResultOutputter() = default;
ResultOutputter::~ResultOutputter() = default;

} // namespace output
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::output::ResultOutputter::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::output::ResultOutputter::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::output::ResultOutputter );
CEREAL_REGISTER_TYPE( protocols::jd3::output::ResultOutputter )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_output_ResultOutputter )
#endif // SERIALIZATION
