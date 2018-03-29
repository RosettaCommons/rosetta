// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_task_operations/RestrictToLoops.cc
/// @brief
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit headers
#include <protocols/simple_task_operations/RestrictToLoops.hh>
#include <protocols/simple_task_operations/RestrictToLoopsCreator.hh>

// Package headers
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/loops_main.hh>

// Project headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// C++ headers
#include <string>

namespace protocols {
namespace simple_task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using namespace std;
using namespace core;
using core::pose::Pose;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using protocols::loops::Loops;
using protocols::loops::LoopsCOP;
using utility::tag::TagCOP;
using basic::datacache::DataMap;

TaskOperationOP RestrictToLoopsCreator::create_task_operation() const {
	return TaskOperationOP( new RestrictToLoops );
}

void RestrictToLoopsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToLoops::provide_xml_schema( xsd );
}

std::string RestrictToLoopsCreator::keyname() const
{
	return RestrictToLoops::keyname();
}

RestrictToLoops::RestrictToLoops() : parent() {
	init();
}

RestrictToLoops::RestrictToLoops( RestrictToLoops const & rhs ) : parent(rhs) {
	copy(*this, rhs);
}

RestrictToLoops::~RestrictToLoops() = default;

RestrictToLoops & RestrictToLoops::operator = ( RestrictToLoops const & rhs ) {
	if ( this == &rhs ) return *this;
	parent::operator=(rhs);
	copy(*this, rhs);
	return *this;
}

TaskOperationOP RestrictToLoops::clone() const {
	return TaskOperationOP( new RestrictToLoops(*this) );
}

void RestrictToLoops::parse_tag( TagCOP tag, DataMap & ) {

	// Parse the 'design' option.

	set_design_loop(
		tag->getOption<bool>( "design", design_loop() ));

	set_restrict_only_design_to_loops(
		tag->getOption<bool>( "restrict_only_design_to_loops", restrict_only_design_to_loops() ));

	// Parse the 'loops_file' option.

	if ( tag->hasOption( "loops_file" ) ) {
		set_loops_from_file(
			tag->getOption<string>( "loops_file" ));
	}
}

void RestrictToLoops::provide_xml_schema( XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	provide_attributes( attributes );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO" );
}

void RestrictToLoops::provide_attributes( AttributeList & attributes ) {
	attributes
		+ XMLSchemaAttribute( "design", xsct_rosetta_bool , "XRW TO DO" )
		+ XMLSchemaAttribute( "restrict_only_design_to_loops", xsct_rosetta_bool , "XRW TO DO" )
		+ XMLSchemaAttribute( "loops_file", xs_string , "XRW TO DO" );
}

void RestrictToLoops::init() {
	set_design_loop( false );
	set_restrict_only_design_to_loops( false );
	set_loops( nullptr );
}

void RestrictToLoops::copy( RestrictToLoops & lhs, RestrictToLoops const & rhs ) {
	lhs.design_loops_ = rhs.design_loops_;
	lhs.restrict_only_design_ = rhs.restrict_only_design_;
	lhs.loops_ = loops::LoopsCOP( loops::LoopsOP( new Loops( *rhs.loops_ ) ) );
}

void RestrictToLoops::apply( Pose const & pose, PackerTask & task ) const {
	apply_helper( pose, task, false, 0, false );
}

void RestrictToLoops::apply_helper(
	Pose const & pose,
	PackerTask & task,
	bool include_neighbors,
	Real cutoff_distance,
	bool design_neighbors ) const {

	if ( ! loops() ) { return; }

	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	// Decide which residues should be packed.

	utility::vector1<bool> is_packable( pose.size(), false );
	select_loop_residues( pose, *loops(), include_neighbors, is_packable, cutoff_distance );
	core::pose::symmetry::make_residue_mask_symmetric( pose, is_packable );

	// Decide which residues should be designed.
	utility::vector1<bool> is_designable;
	if ( design_neighbors ) {
		is_designable = is_packable;
	} else if ( restrict_only_design_ || design_loop() ) {
		is_designable.resize( pose.size(), false );
		loops()->transfer_to_residue_vector( is_designable, true );
	} else {
		is_designable.resize( pose.size(), false );
	}

	// Update the packer task.
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( ! is_designable[ i ] ) {
			turn_off_design.include_residue( i );
		}
		if ( ! is_packable[ i ] ) {
			turn_off_packing.include_residue( i );
		}
	}

	turn_off_design.apply( pose, task );

	if ( ! restrict_only_design_ ) {
		turn_off_packing.apply( pose, task );
	}

}

bool RestrictToLoops::design_loop() const {
	return design_loops_;
}

void RestrictToLoops::set_design_loop( bool design_loop ) {
	design_loops_ = design_loop;
}

LoopsCOP RestrictToLoops::loops() const {
	return loops_;
}

bool RestrictToLoops::restrict_only_design_to_loops() const {
	return restrict_only_design_;
}

void RestrictToLoops::set_restrict_only_design_to_loops( bool restrict_only_design ) {
	restrict_only_design_ = restrict_only_design;

}

void RestrictToLoops::set_loops( LoopsCOP loops ) {
	loops_ = loops;
}


void RestrictToLoops::set_loops_from_file( string loops_file ) {
	loops_ = loops::LoopsCOP( loops::LoopsOP( new Loops( loops_file ) ) );
}


} //namespace simple_task_operations
} //namespace protocols
