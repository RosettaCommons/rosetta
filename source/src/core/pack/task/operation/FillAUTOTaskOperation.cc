// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/FillAUTOTaskOperation.cc
/// @brief fills the AUTO behavior for all residues in Task. Useful if a protocol expects AUTO-style resfile, but no resfile present.
/// @author Steven Lewis (smlewi@gmail.com)

#include <core/pack/task/operation/FillAUTOTaskOperation.hh>
#include <core/pack/task/operation/FillAUTOTaskOperationCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.pack.task.operation.FillAUTOTaskOperation" );

namespace core {
namespace pack {
namespace task {
namespace operation {

FillAUTOTaskOperation::FillAUTOTaskOperation():
	TaskOperation()
{}

FillAUTOTaskOperation::~FillAUTOTaskOperation() = default;

TaskOperationOP
FillAUTOTaskOperation::clone() const {
	return TaskOperationOP( new FillAUTOTaskOperation( *this ) );
}

FillAUTOTaskOperation::FillAUTOTaskOperation( FillAUTOTaskOperation const & /*src*/ ) = default;

//parse_tag is a no-op
void
FillAUTOTaskOperation::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap& ){}

void
FillAUTOTaskOperation::apply( core::pose::Pose const & /*pose*/, core::pack::task::PackerTask & task ) const {

	for ( core::Size i( 1 ); i <= task.total_residue(); ++i ) {
		task.nonconst_residue_task( i ).add_behavior("AUTO");
	}
}

std::string
FillAUTOTaskOperation::keyname() {
	return "FillAUTOTaskOperation";
}

void
FillAUTOTaskOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	task_op_schema_empty( xsd, keyname(), "Adds the AUTO behavior for all residues, as if the resfile had specified default AUTO.  Useful for some protocols that expect an AUTO-enabled resfile (in particular CoupledMoves" );
}

core::pack::task::operation::TaskOperationOP
FillAUTOTaskOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new FillAUTOTaskOperation );
}

std::string
FillAUTOTaskOperationCreator::keyname() const
{
	return FillAUTOTaskOperation::keyname();
}

void
FillAUTOTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FillAUTOTaskOperation::provide_xml_schema( xsd );
}


} //core
} //pack
} //task
} //operation
