// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/task_operations/RestrictToCDRsAndNeighbors.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/task_operations/RestrictToCDRsAndNeighbors.hh>
#include <protocols/antibody/task_operations/RestrictToCDRsAndNeighborsCreator.hh>


#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.task_operationsRestrictToCDRsAndNeighbors");

namespace protocols {
namespace antibody {
namespace task_operations {
using namespace core::pack::task::operation;
using namespace protocols::loops;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility::tag;

RestrictToCDRsAndNeighbors::RestrictToCDRsAndNeighbors():
	TaskOperation(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

RestrictToCDRsAndNeighbors::RestrictToCDRsAndNeighbors(AntibodyInfoCOP ab_info):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
}

RestrictToCDRsAndNeighbors::RestrictToCDRsAndNeighbors(AntibodyInfoCOP ab_info, utility::vector1<bool> const & cdrs):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	cdrs_ = cdrs;
}

RestrictToCDRsAndNeighbors::RestrictToCDRsAndNeighbors(AntibodyInfoCOP ab_info, utility::vector1<bool> const & cdrs, bool allow_cdr_design):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	set_cdrs( cdrs );
	design_cdrs_ = allow_cdr_design;
}

RestrictToCDRsAndNeighbors::RestrictToCDRsAndNeighbors(
	AntibodyInfoCOP ab_info,
	utility::vector1<bool> const & cdrs,
	bool allow_cdr_design,
	bool allow_neighbor_framework_design,
	bool allow_neighbor_antigen_design):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	set_cdrs( cdrs );
	design_cdrs_ = allow_cdr_design;
	design_framework_ = allow_neighbor_framework_design;
	design_antigen_ = allow_neighbor_antigen_design;
}

RestrictToCDRsAndNeighbors::RestrictToCDRsAndNeighbors(RestrictToCDRsAndNeighbors const & src):
	TaskOperation(src),
	cdrs_(src.cdrs_),
	neighbor_dis_(src.neighbor_dis_),
	design_cdrs_(src.design_cdrs_),
	design_antigen_(src.design_antigen_),
	design_framework_(src.design_framework_),
	stem_size_(src.stem_size_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_)
{
	if ( src.ab_info_ ) ab_info_ = AntibodyInfoOP( new AntibodyInfo( *src.ab_info_ ));
}

TaskOperationOP
RestrictToCDRsAndNeighbors::clone() const {
	return TaskOperationOP(new RestrictToCDRsAndNeighbors( *this ));
}

RestrictToCDRsAndNeighbors::~RestrictToCDRsAndNeighbors() {}

void
RestrictToCDRsAndNeighbors::set_defaults() {
	cdrs_.clear();
	cdrs_.resize(8, true);
	cdrs_[ l4 ] = false;
	cdrs_[ h4 ] = false;

	neighbor_dis_ = 6.0;
	design_cdrs_ = false;
	design_antigen_ = false;
	design_framework_ = false;
	stem_size_ = 0;

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::input_ab_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);


}

void
RestrictToCDRsAndNeighbors::parse_tag(utility::tag::TagCOP tag, basic::datacache::DataMap&){
	if ( tag->hasOption("cdrs") ) {
		cdrs_ = get_cdr_bool_from_tag(tag, "cdrs", true /* include cdr4 */);
	}

	neighbor_dis_ = tag->getOption< core::Real >("neighbor_dis", neighbor_dis_);
	design_cdrs_ = tag->getOption< bool >("design_cdrs", design_cdrs_);
	design_antigen_ = tag->getOption< bool >("design_antigen", design_antigen_);
	design_framework_ = tag->getOption< bool >("design_framework", design_framework_);
	stem_size_ = tag->getOption< core::Size >("stem_size", stem_size_);

	if ( tag->hasOption("cdr_definition") && tag->hasOption("input_ab_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("input_ab_scheme"));

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("input_ab_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme. These can also be set via cmd line options of the same name." << std::endl;

	}

}

std::string RestrictToCDRsAndNeighbors::keyname() { return "RestrictToCDRsAndNeighbors"; }

void RestrictToCDRsAndNeighbors::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes_for_get_cdr_bool_from_tag(attributes, "cdrs");

	attributes
		+ XMLSchemaAttribute::attribute_w_default( "neighbor_dis", xsct_real, "XRW TO DO", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "design_cdrs", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "design_antigen", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "design_framework", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "stem_size", xsct_non_negative_integer, "XRW TO DO", "0" )
		+ XMLSchemaAttribute( "cdr_definition", xs_string ,
		"cdr_definition requires input_ab_scheme to be set" )
		+ XMLSchemaAttribute( "input_ab_scheme", xs_string ,
		"input_ab_scheme requires cdr_definition to be set" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

void
RestrictToCDRsAndNeighbors::set_cdrs(const utility::vector1<bool>& cdrs) {
	cdrs_ = cdrs;
	if ( cdrs.size() < CDRNameEnum_proto_total ) {
		for ( core::Size i = cdrs.size() +1; i <= CDRNameEnum_proto_total; ++i ) {
			cdrs_.push_back( false );
		}
	}

	assert(cdrs_.size() == 8);
}

void
RestrictToCDRsAndNeighbors::set_cdr_only(CDRNameEnum cdr) {
	cdrs_.clear();
	cdrs_.resize(CDRNameEnum_proto_total, false);
	cdrs_[ cdr ] = true;
}

void
RestrictToCDRsAndNeighbors::set_allow_design_cdr(bool allow_cdr_design) {
	design_cdrs_ = allow_cdr_design;
}

void
RestrictToCDRsAndNeighbors::set_allow_design_neighbor_framework(bool allow_framework_design){
	design_framework_ = allow_framework_design;
}

void
RestrictToCDRsAndNeighbors::set_allow_design_neighbor_antigen(bool allow_antigen_design){
	design_antigen_ = allow_antigen_design;
}

void
RestrictToCDRsAndNeighbors::set_neighbor_distance(core::Real neighbor_dis) {
	neighbor_dis_ = neighbor_dis;
}

void
RestrictToCDRsAndNeighbors::set_stem_size(core::Size stem_size){
	stem_size_ = stem_size;
}

void
RestrictToCDRsAndNeighbors::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, cdr_definition_));
	} else {
		local_ab_info = ab_info_->clone();
	}

	LoopsOP cdr_loops_and_stems = get_cdr_loops(local_ab_info, pose, cdrs_, stem_size_);

	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	utility::vector1<bool> loops_and_neighbors( pose.size(), false );
	select_loop_residues( pose, *cdr_loops_and_stems, true , loops_and_neighbors, neighbor_dis_ );

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! loops_and_neighbors[ i ] ) {
			turn_off_packing.include_residue( i );
		} else {
			AntibodyRegionEnum region = local_ab_info->get_region_of_residue( pose, i );

			if ( region == cdr_region && !design_cdrs_ ) {
				turn_off_design.include_residue( i );
			} else if ( region == framework_region && ! design_framework_ ) {
				turn_off_design.include_residue( i );
			} else if ( region == antigen_region && ! design_antigen_ ) {
				turn_off_design.include_residue( i );
			}
		}

	}

	turn_off_design.apply(pose, task);
	turn_off_packing.apply(pose, task);
}

core::pack::task::operation::TaskOperationOP
RestrictToCDRsAndNeighborsCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictToCDRsAndNeighbors );
}

void RestrictToCDRsAndNeighborsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToCDRsAndNeighbors::provide_xml_schema( xsd );
}

std::string RestrictToCDRsAndNeighborsCreator::keyname() const
{
	return RestrictToCDRsAndNeighbors::keyname();
}

} //task_operations
} //antibody
} //protocols









