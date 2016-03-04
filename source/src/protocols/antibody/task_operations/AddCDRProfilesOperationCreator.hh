// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/task_operations/AddCDRProfileSetsOperation.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_task_operations_AddCDRProfilesOperationCreator_hh
#define INCLUDED_protocols_antibody_task_operations_AddCDRProfilesOperationCreator_hh

#include <core/pack/task/operation/TaskOperationCreator.hh>

namespace protocols {
namespace antibody {
namespace task_operations {



class AddCDRProfilesOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "ddCDRProfilesOperation"; }
	//core::pack::task::operation::TaskOperationOP clone() const;
};

} //task_operations
} //antibody
} //protocols

#endif //INCLUDED_AddCDRProfilesOperationCreator_hh


