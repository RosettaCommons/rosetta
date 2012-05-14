// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   devel/matdes/STMStoredTask.hh
/// @brief
/// @author Neil King (neilking@uw.edu)


#ifndef INCLUDED_devel_matdes_STMStoredTask_hh
#define INCLUDED_devel_matdes_STMStoredTask_hh

#include <devel/matdes/STMStoredTask.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <map>

namespace devel {
namespace matdes {

class STMStoredTask: public basic::datacache::CacheableData {
public:

	// constructors
	STMStoredTask();
	virtual basic::datacache::CacheableDataOP clone() const;
	virtual basic::datacache::CacheableDataOP fresh_instance() const;

	// setter
	void set_task( core::pack::task::PackerTaskOP const task, std::string const task_name );

	// getters
	core::pack::task::PackerTaskOP get_task( std::string task_name ) const;
	bool has_task( std::string const task_name ) const;

private:
	std::map< std::string, core::pack::task::PackerTaskOP > tasks_;

};

} // matdes
} // devel

#endif
