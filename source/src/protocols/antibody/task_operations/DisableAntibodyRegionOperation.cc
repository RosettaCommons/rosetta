// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/task_operations/DisableAntibodyRegionOperation.cc
/// @brief Task Operation to disable a region of an antibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/task_operations/DisableAntibodyRegionOperation.hh>
#include <protocols/antibody/task_operations/DisableAntibodyRegionOperationCreator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.task_operations.DisableAntibodyRegionOperation");

namespace protocols {
namespace antibody {
namespace task_operations {
using namespace core::pack::task::operation;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility::tag;

DisableAntibodyRegionOperation::DisableAntibodyRegionOperation():
	TaskOperation(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

DisableAntibodyRegionOperation::DisableAntibodyRegionOperation(AntibodyInfoCOP ab_info):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
}

DisableAntibodyRegionOperation::DisableAntibodyRegionOperation(AntibodyInfoCOP ab_info, AntibodyRegionEnum region):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	region_ = region;
}

DisableAntibodyRegionOperation::DisableAntibodyRegionOperation(AntibodyInfoCOP ab_info, AntibodyRegionEnum region, bool disable_packing_and_design):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	region_ = region;
	disable_packing_and_design_ = disable_packing_and_design;
}

DisableAntibodyRegionOperation::~DisableAntibodyRegionOperation() {}

TaskOperationOP
DisableAntibodyRegionOperation::clone() const {
	return TaskOperationOP( new DisableAntibodyRegionOperation( *this));
}

DisableAntibodyRegionOperation::DisableAntibodyRegionOperation(DisableAntibodyRegionOperation const & src):
	TaskOperation(src),
	ab_info_(src.ab_info_),
	region_(src.region_),
	disable_packing_and_design_(src.disable_packing_and_design_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_)
{

}

void
DisableAntibodyRegionOperation::parse_tag(utility::tag::TagCOP tag, basic::datacache::DataMap&){
	AntibodyEnumManager manager = AntibodyEnumManager();

	if ( tag->hasOption("region") ) {
		region_ = manager.antibody_region_string_to_enum(tag->getOption< std::string >( "region" ));
	}
	disable_packing_and_design_ = tag->getOption< bool >( "disable_packing_and_design", disable_packing_and_design_);

	if ( tag->hasOption("cdr_definition") && tag->hasOption("numbering_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("numbering_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;

	}

}

void
DisableAntibodyRegionOperation::set_region(AntibodyRegionEnum region){
	region_ = region;
}

void
DisableAntibodyRegionOperation::set_disable_packing_and_design(bool disable_packing_and_design){
	disable_packing_and_design_ = disable_packing_and_design;
}

void
DisableAntibodyRegionOperation::set_defaults() {
	region_ = cdr_region;
	disable_packing_and_design_ = true;

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::numbering_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);
}


void
DisableAntibodyRegionOperation::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, cdr_definition_));
	} else {
		local_ab_info = ab_info_->clone();
	}

	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( local_ab_info->get_region_of_residue(pose, i) == region_ ) {
			turn_off_packing.include_residue(i);
			turn_off_design.include_residue(i);
		}
	}

	if ( disable_packing_and_design_ ) {
		turn_off_packing.apply( pose, task );
	} else {
		turn_off_design.apply( pose, task );
	}

}

std::string DisableAntibodyRegionOperation::keyname() { return "DisableAntibodyRegionOperation"; }

void DisableAntibodyRegionOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute( "region", xs_string )
		+ XMLSchemaAttribute::attribute_w_default(  "disable_packing_and_design", xs_boolean, "true" )
		+ XMLSchemaAttribute( "cdr_definition", xs_string )
		+ XMLSchemaAttribute( "numbering_scheme", xs_string );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}


core::pack::task::operation::TaskOperationOP
DisableAntibodyRegionOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DisableAntibodyRegionOperation );
}

void DisableAntibodyRegionOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisableAntibodyRegionOperation::provide_xml_schema( xsd );
}

std::string DisableAntibodyRegionOperationCreator::keyname() const
{
	return DisableAntibodyRegionOperation::keyname();
}

} //task_operations
} //antibody
} //protocols









