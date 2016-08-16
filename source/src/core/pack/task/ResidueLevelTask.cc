// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResidueLevelTask.cc
/// @brief  Task class to describe packer's behavior
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/task/ResidueLevelTask.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {

ResidueLevelTask::~ResidueLevelTask() {}

}
}
}


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::ResidueLevelTask::save( Archive & ) const {
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::ResidueLevelTask::load( Archive & ) {
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::ResidueLevelTask );
CEREAL_REGISTER_TYPE( core::pack::task::ResidueLevelTask )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_ResidueLevelTask )
#endif // SERIALIZATION
