// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/STMStoredTask.cc
/// @brief  CacheableData wrapper for PackerTask storage
/// @author Neil King (neilking@uw.edu)

// Unit headers
#include <protocols/toolbox/task_operations/STMStoredTask.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/chemical/ResidueConnection.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace task_operations {

static basic::Tracer TR( "protocols.toolbox.task_operations.STMStoredTask" );

// @brief default constructor
STMStoredTask::STMStoredTask() {}

// @brief copy constructor
STMStoredTask::STMStoredTask(const STMStoredTask & rval) :
	basic::datacache::CacheableData( rval ),
	tasks_( rval.tasks_)
{ }

// @brief setter
void STMStoredTask::set_task( core::pack::task::PackerTaskOP task, std::string const & task_name ) {
	tasks_[ task_name ] = task;
}

// @brief getter
core::pack::task::PackerTaskOP STMStoredTask::get_task( const std::string & task_name ) const { return tasks_.find( task_name )->second; }

// @brief check to see if the task you're interested in is in the object
bool STMStoredTask::has_task( const std::string & task_name ) const { return ( tasks_.find( task_name ) != tasks_.end() ); }

basic::datacache::CacheableDataOP
STMStoredTask::clone() const {
	return basic::datacache::CacheableDataOP( new STMStoredTask ( *this ) );
}

basic::datacache::CacheableDataOP
STMStoredTask::fresh_instance() const {
	return basic::datacache::CacheableDataOP( new STMStoredTask );
}

} // task_operations
} // toolbox
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::task_operations::STMStoredTask::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( tasks_ ) ); // std::map<std::string, core::pack::task::PackerTaskOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::task_operations::STMStoredTask::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( tasks_ ); // std::map<std::string, core::pack::task::PackerTaskOP>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::task_operations::STMStoredTask );
CEREAL_REGISTER_TYPE( protocols::toolbox::task_operations::STMStoredTask )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_task_operations_STMStoredTask )
#endif // SERIALIZATION
