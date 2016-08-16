// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperations.fwd.hh
/// @brief  Forward declaration for core-level (very general) derived classes that wrap widely-used methods of the ResidueLevelTask interface. These are used by higher-level TaskOperations that allow the user to configure the behavior of PackerTasks that are created by TaskFactory.
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperations_fwd_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperations_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictToRepackingRLT;
class RestrictAbsentCanonicalAASRLT;
class DisallowIfNonnativeRLT;
class PreventRepackingRLT;
class AddBehaviorRLT;

typedef utility::pointer::shared_ptr< RestrictToRepackingRLT > RestrictToRepackingRLTOP;
typedef utility::pointer::shared_ptr< RestrictAbsentCanonicalAASRLT > RestrictAbsentCanonicalAASRLTOP;
typedef utility::pointer::shared_ptr< DisallowIfNonnativeRLT > DisallowIfNonnativeRLTOP;
typedef utility::pointer::shared_ptr< PreventRepackingRLT > PreventRepackingRLTOP;
typedef utility::pointer::shared_ptr< AddBehaviorRLT > AddBehaviorRLTOP;

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
