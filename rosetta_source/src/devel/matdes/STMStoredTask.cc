// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   devel/matdes/STMStoredTask.cc
/// @brief	CacheableData wrapper for PackerTask storage
/// @author Neil King (neilking@uw.edu)

// Unit headers
#include <devel/matdes/STMStoredTask.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

namespace devel {
namespace matdes {

static basic::Tracer TR("devel.matdes.STMStoredTask");

// @brief default constructor
STMStoredTask::STMStoredTask() {}

// @brief setter
void STMStoredTask::set_task( core::pack::task::PackerTaskOP task, std::string task_name ) { tasks_[ task_name ] = task; }

// @brief getter
core::pack::task::PackerTaskOP STMStoredTask::get_task( const std::string task_name ) const { return tasks_.find( task_name )->second; }

// @brief check to see if the task you're interested in is in the object
bool STMStoredTask::has_task( const std::string task_name ) const { return ( tasks_.find( task_name ) != tasks_.end() ); }

basic::datacache::CacheableDataOP
STMStoredTask::clone() const {
  return new STMStoredTask ( *this );
}

basic::datacache::CacheableDataOP
STMStoredTask::fresh_instance() const {
  return new STMStoredTask;
}

} // matdes
} // devel
