// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RetrieveStoredTaskOperation.hh
/// @brief  Retrieves a stored task from the pose's CacheableData. Must be used in
///         conjunction with the StoredTaskMover. Allows storage/retrieval of a task
///         so that particular sets of residues can be stably addressed throughout
///         the entirety of a RosettaScripts protocol.
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_RetrieveStoredTaskOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RetrieveStoredTaskOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RetrieveStoredTaskOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RetrieveStoredTaskOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RetrieveStoredTaskOperation : public core::pack::task::operation::TaskOperation {
public:

	RetrieveStoredTaskOperation();
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	virtual ~RetrieveStoredTaskOperation();

	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	virtual void parse_tag( TagCOP, DataMap & );
	virtual void parse_def( utility::lua::LuaObject const & def);

private:
	std::string task_name_;

};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif // INCLUDED_protocols_toolbox_task_operations_RetrieveStoredTaskOperation_HH
