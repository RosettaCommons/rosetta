// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictToCDRH3Loop.cc
/// @brief
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit headers
#include <protocols/toolbox/task_operations/RestrictToCDRH3Loop.hh>
#include <protocols/toolbox/task_operations/RestrictToCDRH3LoopCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers
#include <utility/exit.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;

using core::pose::Pose;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using namespace utility::tag;
using utility::tag::TagCOP;


RestrictToCDRH3Loop::RestrictToCDRH3Loop() : parent() {}

RestrictToCDRH3Loop::RestrictToCDRH3Loop( RestrictToCDRH3Loop const & ) : parent() {}

RestrictToCDRH3Loop::~RestrictToCDRH3Loop() {}

TaskOperationOP RestrictToCDRH3Loop::clone() const
{
	return TaskOperationOP( new RestrictToCDRH3Loop( *this ) );
}

void RestrictToCDRH3Loop::apply( Pose const & pose, PackerTask & task ) const
{
	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_on_packing;

	for ( Size residue_number = 1; residue_number <= pose.total_residue(); ++residue_number ) {
		if ( residue_is_in_h3_loop( pose, residue_number ) ) {
			turn_on_packing.include_residue( residue_number );
		} else {
			turn_off_packing.include_residue( residue_number );
		}
	}

	turn_off_packing.apply( pose, task );
	turn_on_packing.apply( pose, task );
}

bool RestrictToCDRH3Loop::residue_is_in_h3_loop( Pose const & pose, Size residue_number ) const
{
	Size const pose_numbered_h3_loop_start( pose.pdb_info()->pdb2pose( heavy_chain, pdb_numbered_h3_loop_start )  );
	Size const pose_numbered_h3_loop_end( pose.pdb_info()->pdb2pose( heavy_chain, pdb_numbered_h3_loop_end ) );

	return ( residue_number >= pose_numbered_h3_loop_start ) && ( residue_number <= pose_numbered_h3_loop_end );
}

// No parse_tag or other data?
void RestrictToCDRH3Loop::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}


TaskOperationOP RestrictToCDRH3LoopCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictToCDRH3Loop );
}

void RestrictToCDRH3LoopCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToCDRH3Loop::provide_xml_schema( xsd );
}

std::string RestrictToCDRH3LoopCreator::keyname() const
{
	return RestrictToCDRH3Loop::keyname();
}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
