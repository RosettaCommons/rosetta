// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/PreventChainFromRepackingOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
// AUTO-REMOVED #include <set>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.PreventChainFromRepackingOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;

PreventChainFromRepackingOperation::PreventChainFromRepackingOperation() {}

PreventChainFromRepackingOperation::PreventChainFromRepackingOperation( core::Size const chain )
	: parent(), chain_( chain )
{
}

PreventChainFromRepackingOperation::~PreventChainFromRepackingOperation() {}

core::pack::task::operation::TaskOperationOP
PreventChainFromRepackingOperationCreator::create_task_operation() const
{
	return new PreventChainFromRepackingOperation;
}

core::pack::task::operation::TaskOperationOP PreventChainFromRepackingOperation::clone() const
{
	return new PreventChainFromRepackingOperation( *this );
}

void
PreventChainFromRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	runtime_assert(chain_);
  if (chain()>pose.conformation().num_chains()) {
		utility_exit_with_message("Number of chains in pose is smaller than the number defined in the xml under \"chain=\"! Aborting!");
  }
	core::Size const chain_begin( pose.conformation().chain_begin( chain_ ) );
	core::Size const chain_end( pose.conformation().chain_end( chain_ ) );

	core::pack::task::operation::PreventRepacking pp;
	for( core::Size i( chain_begin ); i<=chain_end; ++i )
		pp.include_residue( i );

	pp.apply( pose, task );
}

void
PreventChainFromRepackingOperation::chain( core::Size const chain )
{
	runtime_assert( chain );
	chain_ = chain;
}

core::Size
PreventChainFromRepackingOperation::chain() const
{
	return chain_;
}

void
PreventChainFromRepackingOperation::parse_tag( TagCOP tag , DataMap & )
{
  chain( tag->getOption< core::Size >( "chain", 1 ) );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
