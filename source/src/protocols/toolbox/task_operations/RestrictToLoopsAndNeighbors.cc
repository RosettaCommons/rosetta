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



namespace protocols {
namespace toolbox {
namespace task_operations {

using core::pose::Pose;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using utility::tag::TagCOP;

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

///@brief default constructor
RestrictToLoopsAndNeighbors::RestrictToLoopsAndNeighbors() : parent()
{
	init();
}

///@brief copy constructor
RestrictToLoopsAndNeighbors::RestrictToLoopsAndNeighbors( RestrictToLoopsAndNeighbors const & rhs ) : parent(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

///@brief assignment operator
RestrictToLoopsAndNeighbors & RestrictToLoopsAndNeighbors::operator=( RestrictToLoopsAndNeighbors const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	parent::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
RestrictToLoopsAndNeighbors::~RestrictToLoopsAndNeighbors() {}

//@brief clone operator, calls the copy constructor
TaskOperationOP
RestrictToLoopsAndNeighbors::clone() const
{
	return new RestrictToLoopsAndNeighbors( *this );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


void RestrictToLoopsAndNeighbors::apply( Pose const & pose, PackerTask & task ) const
{

	if ( ! loops() ) return;
	
	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	// Get an up-to-date list of the residues that can be packed
	utility::vector1<bool> is_packable( pose.total_residue(), false );
	select_loop_residues( pose, *loops(), include_neighbors(), is_packable, cutoff_distance() );

	core::pose::symmetry::make_residue_mask_symmetric( pose, is_packable ); // does nothing if pose is not symm

	// If we're designing, allow design at loop positions
	utility::vector1< bool > is_designable( pose.total_residue(), false );
	if ( design_loop() ) {
		loops()->transfer_to_residue_vector( is_designable, true );
	}

	for ( Size residue_number = 1; residue_number <= pose.total_residue(); ++residue_number )
	{
		if ( is_packable[ residue_number ] &&  ! is_designable[ residue_number ] )
		{
			turn_off_design.include_residue( residue_number );
		}
		else if ( ! is_packable[ residue_number ] )
		{
			turn_off_packing.include_residue( residue_number );
		}
	}

	turn_off_design.apply( pose, task );
	turn_off_packing.apply( pose, task );
}

void RestrictToLoopsAndNeighbors::init()
{
	set_design_loop( false );
	set_include_neighbors( true );
	set_cutoff_distance( 10.0 );
	set_loops( NULL );
}

void RestrictToLoopsAndNeighbors::init_for_equal_operator_and_copy_constructor( RestrictToLoopsAndNeighbors & lhs, RestrictToLoopsAndNeighbors const & rhs)
{
	// copy all data members from rhs to lhs
	lhs.design_loop_ = rhs.design_loop_;
	lhs.include_neighbors_ = rhs.include_neighbors_;
	lhs.cutoff_distance_ = rhs.cutoff_distance_;
	lhs.loops_ = rhs.loops_;
}

bool RestrictToLoopsAndNeighbors::design_loop() const
{
	return design_loop_;
}

void RestrictToLoopsAndNeighbors::set_design_loop( bool design_loop )
{
	design_loop_ = design_loop;
}

bool RestrictToLoopsAndNeighbors::include_neighbors() const
{
	return include_neighbors_;
}

void RestrictToLoopsAndNeighbors::set_include_neighbors( bool include_neighbors )
{
	include_neighbors_ = include_neighbors;
}

core::Real RestrictToLoopsAndNeighbors::cutoff_distance() const
{
	return cutoff_distance_;
}

void RestrictToLoopsAndNeighbors::set_cutoff_distance( core::Real cutoff_distance )
{
	if ( cutoff_distance >= 0.0 && cutoff_distance <= 10.0 )
	{
		cutoff_distance_ = cutoff_distance;
	}
}

loops::LoopsCOP RestrictToLoopsAndNeighbors::loops() const
{
	return loops_;
}

void RestrictToLoopsAndNeighbors::set_loops( loops::LoopsCOP loops )
{
	loops_ = loops;
}

TaskOperationOP RestrictToLoopsAndNeighborsCreator::create_task_operation() const
{
	return new RestrictToLoopsAndNeighbors;
}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
