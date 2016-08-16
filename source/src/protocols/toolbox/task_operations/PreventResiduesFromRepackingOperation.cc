// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.cc
/// @brief
/// @author Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperationCreator.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>

// Project Headers

#include <core/pack/task/operation/TaskOperations.hh>

#include <boost/foreach.hpp>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// C++ Headers

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.PreventResiduesFromRepackingOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
PreventResiduesFromRepackingOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new PreventResiduesFromRepackingOperation );
}

void PreventResiduesFromRepackingOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PreventResiduesFromRepackingOperation::provide_xml_schema( xsd );
}

std::string PreventResiduesFromRepackingOperationCreator::keyname() const
{
	return PreventResiduesFromRepackingOperation::keyname();
}

PreventResiduesFromRepackingOperation::PreventResiduesFromRepackingOperation() {}

PreventResiduesFromRepackingOperation::PreventResiduesFromRepackingOperation( utility::vector1 < core::Size > residues )
: parent(), residues_( residues ), reference_pdb_id_( "" )
{
}

PreventResiduesFromRepackingOperation::~PreventResiduesFromRepackingOperation() {}

core::pack::task::operation::TaskOperationOP PreventResiduesFromRepackingOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new PreventResiduesFromRepackingOperation( *this ) );
}

void
PreventResiduesFromRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	if ( residues_.size() == 0 ) { // do nothing
		return;
	}

	core::pack::task::operation::PreventRepacking pp;
	for ( core::Size i = 1; i<=residues_.size(); ++i ) {
		pp.include_residue( residues_[i] );
	}

	pp.apply( pose, task );
}


utility::vector1< core::Size >
PreventResiduesFromRepackingOperation::get_residues() const
{
	return residues_;
}

void
PreventResiduesFromRepackingOperation::set_residues( utility::vector1  < core::Size > residues_vec )
{
	runtime_assert( residues_vec.size() != 0 );
	residues_.clear();
	BOOST_FOREACH ( core::Size const item, residues_vec ) {
		residues_.push_back( item );
	}
}


void
PreventResiduesFromRepackingOperation::parse_tag( TagCOP tag , DataMap & )
{
	reference_pdb_id_ = tag->getOption< std::string >( "reference_pdb_id", "" );
	if ( tag->getOption< std::string >( "residues" ) == "-1" ) {
		residues_.clear();
		TR<<"No residues specified to prevent repacking. I'm doing nothing"<<std::endl;
		return;
	}
	unparsed_residues_ = tag->getOption< std::string >( "residues", "" ) ;
	if ( unparsed_residues_ != "" ) {

		core::pose::Pose reference_pose;
		if ( reference_pdb_id_ != "" ) {
			core::import_pose::pose_from_file( reference_pose, reference_pdb_id_ , core::import_pose::PDB_file);
		}

		utility::vector1< std::string > const res_keys( utility::string_split( unparsed_residues_ , ',' ) );
		utility::vector1< core::Size > residues_vector;
		residues_.clear();

		BOOST_FOREACH ( std::string const key, res_keys ) {
			Size res( utility::string2int( key ) );
			if ( reference_pdb_id_ != "" ) {
				res = core::pose::parse_resnum( key, reference_pose );
			}
			TR.Debug<<"Prevent repack on residue: "<< key  <<" which is translated to residue "<<res<<std::endl;
			residues_.push_back( res );
		}
	}
}

void PreventResiduesFromRepackingOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "reference_pdb_id", xs_string, "" )
		+ XMLSchemaAttribute::attribute_w_default(  "residues", xs_string, "" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}


} //namespace protocols
} //namespace toolbox
} //namespace task_operations
