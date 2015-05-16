// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/task_operations/AddFrameworkProfilesOperation.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/task_operations/AddFrameworkProfilesOperation.hh>
#include <protocols/toolbox/task_operations/ResidueProbDesignOperation.hh>
#include <protocols/toolbox/task_operations/ConservativeDesignOperation.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static thread_local basic::Tracer TR("protocols.antibody.task_operations.AddFrameworkProfilesOperation");

namespace protocols {
namespace antibody {
namespace task_operations {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

AddFrameworkProfilesOperation::AddFrameworkProfilesOperation() {

}

AddFrameworkProfilesOperation::~AddFrameworkProfilesOperation(){}

//AddFrameworkProfilesOperation::AddFrameworkProfilesOperation(AddFrameworkProfilesOperation const & src) {

//}

AddFrameworkProfilesOperation::~AddFrameworkProfilesOperation() {
}


} //task_operations
} //antibody
} //protocols











