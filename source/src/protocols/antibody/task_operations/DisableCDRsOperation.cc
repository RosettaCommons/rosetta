// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.task_operations.DisableCDRsOperation");

namespace protocols {
namespace antibody {
namespace task_operations {

using namespace core::pack::task::operation;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;

DisableCDRsOperation::DisableCDRsOperation():
	TaskOperation(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

DisableCDRsOperation::DisableCDRsOperation(AntibodyInfoCOP ab_info):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
}

DisableCDRsOperation::DisableCDRsOperation(AntibodyInfoCOP ab_info, const utility::vector1<bool>& cdrs):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	set_cdrs( cdrs );

}

DisableCDRsOperation::DisableCDRsOperation(
	AntibodyInfoCOP ab_info,
	utility::vector1<bool> const & cdrs,
	bool disable_packing_and_design):

	TaskOperation(),
	ab_info_(ab_info)
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
	std::string numbering_scheme = option [OptionKeys::antibody::numbering_scheme]();
	std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);
	cdr_definition_ = manager.cdr_definition_string_to_enum(cdr_definition);

}

DisableCDRsOperation::~DisableCDRsOperation() {}

DisableCDRsOperation::DisableCDRsOperation(DisableCDRsOperation const & src):
	core::pack::task::operation::TaskOperation( src ),
	ab_info_(src.ab_info_),
	cdrs_(src.cdrs_),
	disable_packing_and_design_(src.disable_packing_and_design_),
	numbering_scheme_(src.numbering_scheme_),
	cdr_definition_(src.cdr_definition_)
{

}



core::pack::task::operation::TaskOperationOP
DisableCDRsOperation::clone() const {
	return TaskOperationOP(new DisableCDRsOperation(*this));
}

void
DisableCDRsOperation::set_cdrs(const utility::vector1<bool>& cdrs){
	cdrs_ = cdrs;
	if ( cdrs.size() < CDRNameEnum_proto_total ) {
		for ( core::Size i = cdrs.size() +1; i <= CDRNameEnum_proto_total; ++i ) {
			cdrs_.push_back( false );
		}
	}

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

	if ( tag->hasOption("cdr_definition") && tag->hasOption("numbering_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		cdr_definition_ = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("numbering_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;

	}

}

void
DisableCDRsOperation::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, cdr_definition_));
	} else {
		local_ab_info = ab_info_->clone();
	}

	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	for ( core::Size i = 1; i <= CDRNameEnum_proto_total; ++i ) {

		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );


		if ( ! cdrs_[ i ] ) continue;
		if ( local_ab_info->is_camelid() && local_ab_info->get_CDR_chain( cdr ) == 'L' ) continue;



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

core::pack::task::operation::TaskOperationOP
DisableCDRsOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DisableCDRsOperation );
}

} //task_operations
} //antibody
} //protocols



