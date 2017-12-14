// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobSummary.cc
/// @brief  The implementation for class protocols::jd3::JobSummary's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/JobSummary.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

JobSummary::JobSummary() = default;
JobSummary::~JobSummary() = default;

} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::JobSummary::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::JobSummary::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::JobSummary );
CEREAL_REGISTER_TYPE( protocols::jd3::JobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_JobSummary )
#endif // SERIALIZATION
