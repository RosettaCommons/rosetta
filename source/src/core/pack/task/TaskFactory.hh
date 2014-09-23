 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/PackerTaskFactory.hh
/// @brief  Factory class for the creation and initialization of PackerTask objects
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_TaskFactory_hh
#define INCLUDED_core_pack_task_TaskFactory_hh

// Unit Headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/types.hh>

// Package Headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>

// Project Headers
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <list>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {

/// @brief  Factory class for the creation and initialization of PackerTask objects
class TaskFactory : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;
	typedef operation::TaskOperation TaskOperation;
	typedef operation::TaskOperationOP TaskOperationOP;
	typedef operation::TaskOperationCOP TaskOperationCOP;
	typedef std::list< TaskOperationOP > OperationList;
	typedef OperationList::const_iterator const_iterator;

public:

	TaskFactory();
	TaskFactory( TaskFactory const & );
	virtual ~TaskFactory();

	TaskFactory const & operator = ( TaskFactory const & );

	virtual TaskFactoryOP clone() const;

	/// @brief  Non static version.
	PackerTaskOP
	create_task_and_apply_taskoperations( pose::Pose const & pose ) const;

	void modify_task( core::pose::Pose const & pose, PackerTaskOP task ) const;

	/// @brief clones the input task, and pushes it back into the list
	void
	push_back( TaskOperationCOP );

	/// @brief begin iterator of operations_
	const_iterator
	begin() const;

	/// @brief end iterator of operations_
	const_iterator
	end() const;

	void
	clear();

/// @brief return the size of the operations_ list
	core::Size size() const;

public:
	/// @brief Static construction of a task
	static
	PackerTaskOP
	create_packer_task( pose::Pose const & );

private:

	void copy_operations( TaskFactory const & src );

private:
	OperationList operations_;

};

} //namespace task
} //namespace pack
} //namespace core


#endif
