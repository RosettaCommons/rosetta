// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictToTerminiOperation.fwd.hh
/// @brief  Restrict to packing only the residues at either or both termini
/// @author Arpit Tandon
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToTerminiOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToTerminiOperationCreator.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project Headers
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers

#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/tag/Tag.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictToTerminiOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using utility::vector1;
using core::Size;
using core::pack::task::operation::TaskOperationOP;
using core::pose::Pose;
using core::pack::task::PackerTask;
using core::pack::task::operation::RestrictToRepackingRLT;

////////////////// Creator ////////////////
TaskOperationOP
RestrictToTerminiOperationCreator::create_task_operation() const {
	return TaskOperationOP( new RestrictToTerminiOperation );
}

void RestrictToTerminiOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToTerminiOperation::provide_xml_schema( xsd );
}

std::string RestrictToTerminiOperationCreator::keyname() const
{
	return RestrictToTerminiOperation::keyname();
}

///////////////// End Creator ////////////


RestrictToTerminiOperation::RestrictToTerminiOperation() :
	chain_(1),
	repack_n_terminus_(true),
	repack_c_terminus_(true)
{}

RestrictToTerminiOperation::RestrictToTerminiOperation(
	Size const chain,
	bool const repack_n_terminus,
	bool const repack_c_terminus) :
	chain_(chain),
	repack_n_terminus_(repack_n_terminus),
	repack_c_terminus_(repack_c_terminus)
{}

RestrictToTerminiOperation::RestrictToTerminiOperation(RestrictToTerminiOperation const & /*src*/) = default;

RestrictToTerminiOperation::~RestrictToTerminiOperation() = default;

TaskOperationOP
RestrictToTerminiOperation::clone() const {
	return TaskOperationOP( new RestrictToTerminiOperation( *this ) );
}

/// @brief restrict to pack only the N and/or C-termini
void
RestrictToTerminiOperation::apply(
	Pose const & pose,
	PackerTask & task
) const {

	if ( chain_ > pose.conformation().num_chains() ) {
		utility_exit_with_message(
			"The pose does not contain the chain you have specified.");
	}

	vector1<bool> repack_residues(pose.size(), false);

	// N-terminus
	if ( repack_n_terminus_ ) {
		Size const n_terminus(pose.conformation().chain_begin(chain_));
		repack_residues[n_terminus] = true;
	}

	// C-terminus
	if ( repack_c_terminus_ ) {
		Size const c_terminus(pose.conformation().chain_end(chain_));
		repack_residues[c_terminus] = true;
	}

	for ( Size i = 1; i <= repack_residues.size(); ++i ) {
		if ( repack_residues[i] ) {
			task.nonconst_residue_task(i).restrict_to_repacking();
		} else {
			task.nonconst_residue_task(i).prevent_repacking();
		}
	}
}


void
RestrictToTerminiOperation::parse_tag( TagCOP tag , DataMap & )
{
	chain_ = tag->getOption<Size>("chain", 1);
	repack_n_terminus_ = tag->getOption<bool>("repack_n_terminus", true);
	repack_c_terminus_ = tag->getOption<bool>("repack_c_terminus", true);
}

void RestrictToTerminiOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "chain", xsct_non_negative_integer, "XRW TO DO",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "repack_n_terminus", xsct_rosetta_bool, "XRW TO DO",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "repack_c_terminus", xsct_rosetta_bool, "XRW TO DO",  "true"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO" );
}


} //namespace
} //namespace
} //namespace
