// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictChainToRepackingOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// C++ Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictChainToRepackingOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;


RestrictChainToRepackingOperation::RestrictChainToRepackingOperation() {}

RestrictChainToRepackingOperation::RestrictChainToRepackingOperation( core::Size const chain )
: parent(), chain_( chain )
{
}

RestrictChainToRepackingOperation::~RestrictChainToRepackingOperation() {}

core::pack::task::operation::TaskOperationOP RestrictChainToRepackingOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictChainToRepackingOperation( *this ) );
}

void
RestrictChainToRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	runtime_assert(chain_);
	if ( chain()>pose.conformation().num_chains() ) {
		utility_exit_with_message("Number of chains in pose is smaller than the number defined in the xml under \"chain=\"! Aborting!");
	}
	core::Size const chain_begin( pose.conformation().chain_begin( chain_ ) );
	core::Size const chain_end( pose.conformation().chain_end( chain_ ) );

	core::pack::task::operation::RestrictResidueToRepacking rrtr;
	for ( core::Size i( chain_begin ); i<=chain_end; ++i ) {
		rrtr.include_residue( i );
	}

	rrtr.apply( pose, task );
}

void
RestrictChainToRepackingOperation::chain( core::Size const chain )
{
	runtime_assert( chain );
	chain_ = chain;
}

core::Size
RestrictChainToRepackingOperation::chain() const
{
	return chain_;
}

void
RestrictChainToRepackingOperation::parse_tag( TagCOP tag , DataMap & )
{
	chain( tag->getOption< core::Size >( "chain", 1 ) );
}

std::string RestrictChainToRepackingOperation::keyname() { return "RestrictChainToRepacking"; }

void RestrictChainToRepackingOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace::utility::tag;
	AttributeList attributes;

	attributes + XMLSchemaAttribute::attribute_w_default( "chain", xsct_non_negative_integer, "XRW TO DO", "1");

	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO" );

}

core::pack::task::operation::TaskOperationOP
RestrictChainToRepackingOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictChainToRepackingOperation );
}

std::string RestrictChainToRepackingOperationCreator::keyname() const { return RestrictChainToRepackingOperation::keyname(); }

void RestrictChainToRepackingOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictChainToRepackingOperation::provide_xml_schema( xsd );
}


} //namespace protocols
} //namespace toolbox
} //namespace task_operations
