// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/RestrictInterGroupVectorOperation.cc
/// @brief
/// @author  Ben Stranges stranges@unc.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictInterGroupVectorOperation.hh>
#include <protocols/toolbox/task_operations/RestrictInterGroupVectorOperationCreator.hh>
//#include <protocols/toolbox/task_operations/InterfaceTaskOperation.hh> //shouldn't need this
// Project Headers
#include <core/pose/Pose.hh>
#include <core/select/util/interface_vector_calculate.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/selection.hh>
// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// C++ Headers

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictInterGroupVectorOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
RestrictInterGroupVectorOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictInterGroupVectorOperation );
}

void RestrictInterGroupVectorOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictInterGroupVectorOperation::provide_xml_schema( xsd );
}

std::string RestrictInterGroupVectorOperationCreator::keyname() const
{
	return RestrictInterGroupVectorOperation::keyname();
}

/// @brief default constructor
RestrictInterGroupVectorOperation::RestrictInterGroupVectorOperation():
	parent(),
	//need something for groups
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
{ }

/// @brief full constructor
RestrictInterGroupVectorOperation::RestrictInterGroupVectorOperation(
	group_vector const & group,
	core::Real CB_dist_cutoff,
	core::Real nearby_atom_cutoff,
	core::Real vector_angle_cutoff,
	core::Real vector_dist_cutoff
):
	parent(),
	pair_vector_(group),
	CB_dist_cutoff_( CB_dist_cutoff ),
	nearby_atom_cutoff_( nearby_atom_cutoff ),
	vector_angle_cutoff_( vector_angle_cutoff ),
	vector_dist_cutoff_( vector_dist_cutoff )
{}

//@brief convienience contstuctor for one pair
RestrictInterGroupVectorOperation:: RestrictInterGroupVectorOperation(
	group_pair const & one_group,
	core::Real CB_dist_cutoff,
	core::Real nearby_atom_cutoff,
	core::Real vector_angle_cutoff,
	core::Real vector_dist_cutoff
):
	parent(),
	CB_dist_cutoff_( CB_dist_cutoff ),
	nearby_atom_cutoff_( nearby_atom_cutoff ),
	vector_angle_cutoff_( vector_angle_cutoff ),
	vector_dist_cutoff_( vector_dist_cutoff )
{
	pair_vector_.clear();
	pair_vector_.push_back(one_group);
}

/// @brief destructor
RestrictInterGroupVectorOperation::~RestrictInterGroupVectorOperation() {}

/// @details be warned if you use clone that you'll not get a new interface calculator
core::pack::task::operation::TaskOperationOP RestrictInterGroupVectorOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictInterGroupVectorOperation( *this ) );
}

/// @details setters
void RestrictInterGroupVectorOperation::insert_pair( group_pair pair){
	pair_vector_.push_back(pair);
}
void RestrictInterGroupVectorOperation::CB_dist_cutoff( core::Real CB_dist_cutoff){
	CB_dist_cutoff_ = CB_dist_cutoff;
}
void RestrictInterGroupVectorOperation::nearby_atom_cutoff(core::Real nearby_atom_cutoff){
	nearby_atom_cutoff_ = nearby_atom_cutoff;
}
void RestrictInterGroupVectorOperation::vector_angle_cutoff(core::Real vector_angle_cutoff){
	vector_angle_cutoff_ = vector_angle_cutoff;
}
void RestrictInterGroupVectorOperation::vector_dist_cutoff(core::Real vector_dist_cutoff){
	vector_dist_cutoff_ = vector_dist_cutoff;
}

/// @details apply function, uses inherited functionality
void
RestrictInterGroupVectorOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	//First check to make sure we can do something
	if ( pair_vector_.size() < 1 ) {
		utility_exit_with_message( "No residues defined to look at interface between." );
	}
	//make a list of residues allowed
	utility::vector1_bool repack_full(pose.total_residue(), false);
	//itterate over each set of pairs, make a union
	for ( core::Size jj=1 ; jj<= pair_vector_.size(); ++jj ) {//{utility::vector1<group_pair>::const_iterator jj = pair_vector_.begin() ; jj != pair_vector_.end() ; ++jj ){
		group_pair this_pair = pair_vector_[jj];
		utility::vector1_bool repack =
			core::select::util::calc_interacting_vector(
			pose,
			this_pair.first,
			this_pair.second,
			CB_dist_cutoff_,
			nearby_atom_cutoff_,
			vector_angle_cutoff_,
			vector_dist_cutoff_
		);
		//add repack true setting to repack_full
		for ( core::Size ii = 1; ii <= repack.size(); ++ii ) {
			if ( repack[ii] ) {
				repack_full[ii] = true;
			}
		}//end add repack to full_repack
	}//itterate over vector of pairs
	task.restrict_to_residues(repack_full);

}//apply

// /// @details parse_tag
// void
// RestrictInterGroupVectorOperation::parse_tag(utility::tag::TagCOP tag)
// {
// }//parse_tag

// parse tag not defined
void RestrictInterGroupVectorOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

}//task_operations
}//toolbox
}//protocols
