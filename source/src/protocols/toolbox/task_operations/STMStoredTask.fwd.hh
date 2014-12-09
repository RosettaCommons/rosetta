// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   STMStoredTask.fwd.hh
/// @brief
/// @author Neil King (neilking@uw.edu)


#ifndef INCLUDED_protocols_toolbox_task_operations_STMStoredTask_fwd_hh
#define INCLUDED_protocols_toolbox_task_operations_STMStoredTask_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

class STMStoredTask;
typedef utility::pointer::shared_ptr< STMStoredTask > STMStoredTaskOP;
typedef utility::pointer::shared_ptr< STMStoredTask const > STMStoredTaskCOP;

} // task_operations
} // toolbox
} // protocols

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#endif

#endif
