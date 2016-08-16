// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperation.fwd.hh
/// @brief  Forward declaration of a class that performs an operation on ResidueLevelTask
///         usually by a TaskFactory right after the task's construction.
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperation_fwd_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

class ResLvlTaskOperation;
typedef utility::pointer::shared_ptr< ResLvlTaskOperation > ResLvlTaskOperationOP;
typedef utility::pointer::shared_ptr< ResLvlTaskOperation > ResLvlTaskOperationCOP;

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
