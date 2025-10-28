// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictToCDRH3Loop.cc
/// @brief
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit headers
#include <protocols/task_operations/RestrictToCDRH3Loop.hh>
#include <protocols/task_operations/RestrictToCDRH3LoopCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers


namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;

using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using namespace utility::tag;
using utility::tag::TagCOP;
using core::Size;

std::string const RestrictToCDRH3Loop::heavy_chain = "H";

RestrictToCDRH3Loop::RestrictToCDRH3Loop() : parent() {}

RestrictToCDRH3Loop::RestrictToCDRH3Loop( RestrictToCDRH3Loop const & ) : parent() {}

RestrictToCDRH3Loop::~RestrictToCDRH3Loop() = default;

TaskOperationOP RestrictToCDRH3Loop::clone() const
{
	return utility::pointer::make_shared< RestrictToCDRH3Loop >( *this );
}

void RestrictToCDRH3Loop::apply( Pose const & pose, PackerTask & task ) const
{
	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_on_packing;

	for ( core::Size residue_number = 1; residue_number <= pose.size(); ++residue_number ) {
		if ( residue_is_in_h3_loop( pose, residue_number ) ) {
			turn_on_packing.include_residue( residue_number );
		} else {
			turn_off_packing.include_residue( residue_number );
		}
	}

	turn_off_packing.apply( pose, task );
	turn_on_packing.apply( pose, task );
}

bool RestrictToCDRH3Loop::residue_is_in_h3_loop( Pose const & pose, core::Size residue_number ) const
{
	core::Size const pose_numbered_h3_loop_start( pose.pdb_info()->pdb2pose( heavy_chain, pdb_numbered_h3_loop_start )  );
	core::Size const pose_numbered_h3_loop_end( pose.pdb_info()->pdb2pose( heavy_chain, pdb_numbered_h3_loop_end ) );

	return ( residue_number >= pose_numbered_h3_loop_start ) && ( residue_number <= pose_numbered_h3_loop_end );
}

// No parse_tag or other data?
void RestrictToCDRH3Loop::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname(), "XRW TO DO" );
}


TaskOperationOP RestrictToCDRH3LoopCreator::create_task_operation() const
{
	return utility::pointer::make_shared< RestrictToCDRH3Loop >();
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
} //namespace protocols
