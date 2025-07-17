// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/task_operations/DisableCDRsOperation.cc
/// @brief Task Operation to disable specific CDRs
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/task_operations/DisableCDRsOperation.hh>
#include <protocols/antibody/task_operations/DisableCDRsOperationCreator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <basic/Tracer.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static basic::Tracer TR("protocols.antibody.task_operations.DisableCDRsOperation");

namespace protocols {
namespace antibody {
namespace task_operations {

using namespace core::pack::task::operation;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility::tag;

DisableCDRsOperation::DisableCDRsOperation():
	TaskOperation(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

DisableCDRsOperation::DisableCDRsOperation(AntibodyInfoCOP ab_info):
	TaskOperation(),
	ab_info_(std::move(ab_info))
{
	set_defaults();
}

DisableCDRsOperation::DisableCDRsOperation(AntibodyInfoCOP ab_info, const utility::vector1<bool>& cdrs):
	TaskOperation(),
	ab_info_(std::move(ab_info))
{
	set_defaults();
	set_cdrs( cdrs );

}

DisableCDRsOperation::DisableCDRsOperation(
	AntibodyInfoCOP ab_info,
	utility::vector1<bool> const & cdrs,
	bool disable_packing_and_design):

	TaskOperation(),
	ab_info_(std::move(ab_info))
{
	set_defaults();
	set_cdrs( cdrs );
	disable_packing_and_design_ = disable_packing_and_design;
}

void
DisableCDRsOperation::set_defaults() {
	disable_packing_and_design_ = true;
	cdrs_.clear();
	cdrs_.resize(8, true);
	cdrs_[ l4 ] = false;
	cdrs_[ h4 ] = false;

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::input_ab_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);

}

DisableCDRsOperation::~DisableCDRsOperation() = default;

DisableCDRsOperation::DisableCDRsOperation(DisableCDRsOperation const & src):
	core::pack::task::operation::TaskOperation( src ),
	cdrs_(src.cdrs_),
	disable_packing_and_design_(src.disable_packing_and_design_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_)
{
	if ( src.ab_info_ ) ab_info_ = utility::pointer::make_shared< AntibodyInfo >( *src.ab_info_ );
}



core::pack::task::operation::TaskOperationOP
DisableCDRsOperation::clone() const {
	return utility::pointer::make_shared< DisableCDRsOperation >(*this);
}

void
DisableCDRsOperation::set_cdrs(const utility::vector1<bool>& cdrs){
	cdrs_ = cdrs;
	if ( cdrs_.size() < CDRNameEnum_proto_total ) {
		for ( core::Size i = cdrs_.size() +1; i <= CDRNameEnum_proto_total; ++i ) {
			cdrs_.push_back( false );
		}
	}
	debug_assert( cdrs_.size() == 8);

}

void
DisableCDRsOperation::set_cdr_only(CDRNameEnum cdr){
	cdrs_.clear();
	cdrs_.resize(CDRNameEnum_proto_total, false);
	cdrs_[ cdr ] = true;
}

void
DisableCDRsOperation::set_disable_packing_and_design(bool disable_packing_and_design){
	disable_packing_and_design_ = disable_packing_and_design;
}

void
DisableCDRsOperation::parse_tag( TagCOP tag, DataMap & ) {
	AntibodyEnumManager manager = AntibodyEnumManager();

	if ( tag->hasOption("cdrs") ) {
		TR << "Setting CDRs from settings" << std::endl;
		cdrs_ = get_cdr_bool_from_tag(tag, "cdrs", true /* include_cdr4*/);
	}

	disable_packing_and_design_ = tag->getOption< bool >("disable_packing_and_design", disable_packing_and_design_);

	if ( tag->hasOption("cdr_definition") && tag->hasOption("input_ab_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("input_ab_scheme"));

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("input_ab_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;

	}
}


void
DisableCDRsOperation::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = utility::pointer::make_shared< AntibodyInfo >(pose, numbering_scheme_, cdr_definition_);
	} else {
		local_ab_info = ab_info_->clone();
	}

	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	for ( core::Size i = 1; i <= CDRNameEnum_proto_total; ++i ) {

		auto cdr = static_cast<CDRNameEnum>( i );


		if ( ! cdrs_[ i ] ) continue;
		if ( local_ab_info->is_camelid() && local_ab_info->get_CDR_chain( cdr ) == "L" ) continue;



		core::Size start = local_ab_info->get_CDR_start( cdr, pose );
		core::Size end = local_ab_info->get_CDR_end( cdr, pose );

		for ( core::Size resnum = start; resnum <= end; ++resnum ) {
			turn_off_packing.include_residue( resnum );
			turn_off_design.include_residue( resnum );
		}
	}

	if ( disable_packing_and_design_ ) {
		turn_off_packing.apply( pose, task );
	} else {
		turn_off_design.apply( pose, task);
	}
}

std::string DisableCDRsOperation::keyname() { return "DisableCDRsOperation"; }

void DisableCDRsOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes_for_get_cdr_bool_from_tag(attributes, "cdrs");

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "disable_packing_and_design", xsct_rosetta_bool, "Disable packing AND design of these CDRs",  "true"  )
		+ XMLSchemaAttribute( "cdr_definition", xs_string ,
		"The specific CDR definition" )
		+ XMLSchemaAttribute( "input_ab_scheme", xs_string ,
		"The numbering scheme of the antibody" );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n Disable specific CDRs");
}

core::pack::task::operation::TaskOperationOP
DisableCDRsOperationCreator::create_task_operation() const
{
	return utility::pointer::make_shared< DisableCDRsOperation >();
}

void DisableCDRsOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisableCDRsOperation::provide_xml_schema( xsd );
}

std::string DisableCDRsOperationCreator::keyname() const
{
	return DisableCDRsOperation::keyname();
}

} //task_operations
} //antibody
} //protocols



