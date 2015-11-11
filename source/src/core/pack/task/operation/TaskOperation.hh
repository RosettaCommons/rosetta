// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperation.hh
/// @brief  An operation to perform on a packer task --
///         usually, by a PackerTaskFactory right after the task's construction
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// only generalized base classes go here. TaskOperations that actually do things do not belong here.

#ifndef INCLUDED_core_pack_task_operation_TaskOperation_hh
#define INCLUDED_core_pack_task_operation_TaskOperation_hh

// Unit Headers
#include <core/pack/task/operation/TaskOperation.fwd.hh>

// VC++ will not compile if these are forward declaration headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

class TaskOperation : public utility::pointer::ReferenceCount
{
public:
	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef pose::Pose Pose;

public:
	virtual ~TaskOperation();

	/// @brief Create another task operation of the type matching the most-derived
	/// version of the class.
	virtual TaskOperationOP clone() const = 0;

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	virtual void apply( Pose const &, PackerTask & ) const = 0;

	/// @brief Used to parse an xml-like tag to load parameters and properties.
	virtual void parse_tag( TagCOP, DataMap & );

};

// do not add any derived classes to this file, unless they are generalized abstract base classes and do not actually 'do any work'

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
