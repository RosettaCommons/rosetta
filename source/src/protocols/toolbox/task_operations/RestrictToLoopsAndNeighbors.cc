// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.cc
/// @brief
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit headers
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighborsCreator.hh>

// Package headers
#include <protocols/toolbox/task_operations/RestrictToLoops.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/loops_main.hh>

// Project headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using utility::tag::TagCOP;

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
RestrictToLoopsAndNeighbors::RestrictToLoopsAndNeighbors() : parent()
{
	init();
}

/// @brief copy constructor
RestrictToLoopsAndNeighbors::RestrictToLoopsAndNeighbors( RestrictToLoopsAndNeighbors const & rhs ) : parent(rhs)
{
	init_for_copy( *this, rhs );
}

/// @brief assignment operator
RestrictToLoopsAndNeighbors & RestrictToLoopsAndNeighbors::operator=( RestrictToLoopsAndNeighbors const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	parent::operator=( rhs );
	init_for_copy( *this, rhs );
	return *this;
}

//destructor
RestrictToLoopsAndNeighbors::~RestrictToLoopsAndNeighbors() {}

//@brief clone operator, calls the copy constructor
TaskOperationOP
RestrictToLoopsAndNeighbors::clone() const
{
	return TaskOperationOP( new RestrictToLoopsAndNeighbors( *this ) );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


void RestrictToLoopsAndNeighbors::parse_tag(TagCOP tag, DataMap & data)
{
	parent::parse_tag(tag, data);

	set_include_neighbors(
		tag->getOption<bool>("include_neighbors", include_neighbors()));

	set_design_neighbors(
		tag->getOption<bool>("design_neighbors", design_neighbors()));

	set_cutoff_distance(
		tag->getOption<Real>("cutoff_dist", cutoff_distance()));


}

void RestrictToLoopsAndNeighbors::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	// From parent.
	RestrictToLoops::provide_attributes( attributes );

	attributes.push_back( XMLSchemaAttribute( "include_neighbors", xs_boolean ) );
	attributes.push_back( XMLSchemaAttribute( "design_neighbors", xs_boolean ) );
	attributes.push_back( XMLSchemaAttribute( "cutoff_dist", xs_decimal ) );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

void RestrictToLoopsAndNeighbors::apply( Pose const & pose, PackerTask & task ) const
{
	apply_helper(pose, task, include_neighbors(), cutoff_distance(), design_neighbors());
}

void RestrictToLoopsAndNeighbors::init()
{
	parent::init();
	set_include_neighbors( true );
	set_design_neighbors( false );
	set_cutoff_distance( 10.0 );
}

void RestrictToLoopsAndNeighbors::init_for_copy( RestrictToLoopsAndNeighbors & lhs, RestrictToLoopsAndNeighbors const & rhs)
{
	parent::copy(lhs, rhs);
	lhs.include_neighbors_ = rhs.include_neighbors_;
	lhs.cutoff_distance_ = rhs.cutoff_distance_;
	lhs.design_neighbors_ = rhs.design_neighbors_;
}

bool RestrictToLoopsAndNeighbors::include_neighbors() const
{
	return include_neighbors_;
}

void RestrictToLoopsAndNeighbors::set_include_neighbors( bool include_neighbors )
{
	include_neighbors_ = include_neighbors;
}

bool RestrictToLoopsAndNeighbors::design_neighbors() const
{
	return design_neighbors_;
}

void RestrictToLoopsAndNeighbors::set_design_neighbors(bool design_neighbors)
{
	design_neighbors_ = design_neighbors;
}

core::Real RestrictToLoopsAndNeighbors::cutoff_distance() const
{
	return cutoff_distance_;
}

void RestrictToLoopsAndNeighbors::set_cutoff_distance( core::Real cutoff_distance )
{
	if ( cutoff_distance >= 0.0 && cutoff_distance <= 10.0 ) {
		cutoff_distance_ = cutoff_distance;
	}
}

TaskOperationOP RestrictToLoopsAndNeighborsCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictToLoopsAndNeighbors );
}

void RestrictToLoopsAndNeighborsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToLoopsAndNeighbors::provide_xml_schema( xsd );
}

std::string RestrictToLoopsAndNeighborsCreator::keyname() const
{
	return RestrictToLoopsAndNeighbors::keyname();
}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
