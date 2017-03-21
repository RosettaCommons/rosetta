// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/task_operations/AddCDRProfileSetsOperation.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/task_operations/AddCDRProfileSetsOperation.hh>
#include <protocols/antibody/task_operations/AddCDRProfileSetsOperationCreator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/toolbox/task_operations/MutationSetDesignOperation.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/AA.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.task_operations.AddCDRProfileSetsOperation");

namespace protocols {
namespace antibody {
namespace task_operations {
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace protocols::toolbox::task_operations;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::antibody::design;
using namespace utility::tag;

AddCDRProfileSetsOperation::AddCDRProfileSetsOperation():
	TaskOperation(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

AddCDRProfileSetsOperation::AddCDRProfileSetsOperation(AntibodyInfoCOP ab_info):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
}

AddCDRProfileSetsOperation::AddCDRProfileSetsOperation(AntibodyInfoCOP ab_info, const utility::vector1<bool>& cdrs):
	TaskOperation(),
	ab_info_(ab_info)

{
	set_defaults();
	cdrs_ = cdrs;
}

AddCDRProfileSetsOperation::AddCDRProfileSetsOperation(AntibodyInfoCOP ab_info, const utility::vector1<bool>& cdrs, bool limit_only_to_length):
	TaskOperation(),
	ab_info_(ab_info)
{
	set_defaults();
	cdrs_ = cdrs;
	limit_only_to_length_ = limit_only_to_length;
}

AddCDRProfileSetsOperation::AddCDRProfileSetsOperation(AddCDRProfileSetsOperation const & src):
	TaskOperation(src),
	cdrs_(src.cdrs_),
	limit_only_to_length_(src.limit_only_to_length_),
	force_north_paper_db_(src.force_north_paper_db_),
	use_outliers_(src.use_outliers_),
	cutoff_(src.cutoff_),
	picking_rounds_(src.picking_rounds_),
	keep_task_allowed_aas_(src.keep_task_allowed_aas_),
	include_native_restype_(src.include_native_restype_),
	sequences_(src.sequences_),
	pre_loaded_data_(src.pre_loaded_data_),
	numbering_scheme_(src.numbering_scheme_),
	ignore_light_chain_( src.ignore_light_chain_ )
{
	if ( src.ab_info_ ) ab_info_ = AntibodyInfoOP( new AntibodyInfo( *src.ab_info_));
}

TaskOperationOP
AddCDRProfileSetsOperation::clone() const {
	return TaskOperationOP( new AddCDRProfileSetsOperation( *this));
}

AddCDRProfileSetsOperation::~AddCDRProfileSetsOperation(){}

void
AddCDRProfileSetsOperation::set_defaults(){
	cdrs_.clear();
	cdrs_.resize(8, true);
	limit_only_to_length_ = false;
	force_north_paper_db_ = false;
	use_outliers_ = false;
	cutoff_ = 10;

	//Profile Options
	picking_rounds_ = 1;
	keep_task_allowed_aas_ = false;
	include_native_restype_ = true;
	pre_loaded_data_ = false;
	sequences_.clear();

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::input_ab_scheme]();
	//std::string cdr_definition = option [OptionKeys::antibody::cdr_definition]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);

}

void
AddCDRProfileSetsOperation::parse_tag(utility::tag::TagCOP tag, basic::datacache::DataMap&){
	if ( tag->hasOption("cdrs") ) {
		TR << "Setting CDRs from settings" << std::endl;
		cdrs_ = get_cdr_bool_from_tag(tag, "cdrs", true /* include cdr4 */);
	}

	limit_only_to_length_ = tag->getOption< bool >("limit_only_to_length", limit_only_to_length_);
	force_north_paper_db_ = tag->getOption< bool >("force_north_paper_db", force_north_paper_db_);
	use_outliers_ = tag->getOption< bool >("use_outliers", use_outliers_);

	keep_task_allowed_aas_ = tag->getOption< bool>("add_to_current", keep_task_allowed_aas_);
	include_native_restype_ = tag->getOption< bool>("include_native_restype", include_native_restype_);
	picking_rounds_ = tag->getOption< core::Size >("picking_rounds", picking_rounds_);
	cutoff_ = tag->getOption< core::Size >("cutoff", cutoff_);

	if ( tag->hasOption("input_ab_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("input_ab_scheme"));

	}
	if ( tag->hasOption("cdr_definition") && tag->getOption<std::string>("cdr_definition") != "North" ) {
		TR <<"This operation only works with the North CDR definition." <<std::endl;
	}

}

void AddCDRProfileSetsOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute( "cdrs", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "limit_only_to_length", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "force_north_paper_db", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_outliers", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "add_to_current", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "include_native_restype", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "picking_rounds", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "cutoff", xsct_non_negative_integer, "XRW TO DO", "10" )
		+ XMLSchemaAttribute( "input_ab_scheme", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute( "cdr_definition", xs_string , "XRW TO DO" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}



void
AddCDRProfileSetsOperation::set_cdr_only(CDRNameEnum cdr) {
	cdrs_.clear();
	cdrs_.resize(8, false);
	cdrs_[ cdr ] = true;
}

void
AddCDRProfileSetsOperation::set_cdrs(const utility::vector1<bool>& cdrs) {
	cdrs_ = cdrs;
	if ( cdrs_.size() < CDRNameEnum_proto_total ) {
		for ( core::Size i = cdrs_.size() +1; i <= CDRNameEnum_proto_total; ++i ) {
			cdrs_.push_back( false );
		}
	}
	debug_assert( cdrs_.size() == 8);

}

void
AddCDRProfileSetsOperation::set_limit_only_to_length(bool limit_only_to_length){
	limit_only_to_length_ = limit_only_to_length;
}

void
AddCDRProfileSetsOperation::set_force_north_paper_db(bool force_north_db){
	force_north_paper_db_ = force_north_db;
}


void
AddCDRProfileSetsOperation::set_use_outliers(bool use_outliers) {
	use_outliers_ = use_outliers;
}

void
AddCDRProfileSetsOperation::set_picking_rounds(core::Size rounds){
	picking_rounds_ = rounds;
}

void
AddCDRProfileSetsOperation::set_include_native_type(bool use_native){
	include_native_restype_ = use_native;
}

void
AddCDRProfileSetsOperation::set_add_to_current(bool add_to_current) {
	keep_task_allowed_aas_ = add_to_current;
}

void
AddCDRProfileSetsOperation::set_cutoff(core::Size cutoff){
	cutoff_ = cutoff;
}

utility::vector1<bool>
AddCDRProfileSetsOperation::pre_load_data(const core::pose::Pose& pose){
	if ( ! ab_info_ ) {
		ab_info_ = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, North));
	}

	AntibodyDatabaseManager manager = AntibodyDatabaseManager(ab_info_, force_north_paper_db_);
	manager.set_outlier_use(use_outliers_);
	manager.ignore_light_chain( ignore_light_chain_ );
	
	sequences_ = manager.load_cdr_sequences(cdrs_, pose, limit_only_to_length_);

	pre_loaded_data_ = true;

	utility::vector1<bool> no_data_cdrs(8, false);

	for ( core::Size i = 1; i <= 8; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( ! cdrs_[ i ] ) continue;

		if ( cdr == l4 || cdr == h4 ) {
			TR << "Skipping CDR4 profile data!" << std::endl;
			no_data_cdrs[ i ] = true;
			continue;
		}

		if (  sequences_.count(cdr) == 0 || (  sequences_.count(cdr) != 0  &&  sequences_[ cdr ].size() <= cutoff_ ) ) {
			no_data_cdrs[ i ] = true;
		}
	}
	return no_data_cdrs;
}

void
AddCDRProfileSetsOperation::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{
	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, North));
	} else {
		local_ab_info = ab_info_->clone();
	}

	CDRDBSequenceSet sequences;

	if ( pre_loaded_data_ ) {
		sequences = sequences_;

	} else {
		AntibodyDatabaseManager manager = AntibodyDatabaseManager(local_ab_info, force_north_paper_db_);
		manager.set_outlier_use(use_outliers_);
		manager.ignore_light_chain( ignore_light_chain_ );
		
		sequences = manager.load_cdr_sequences(cdrs_, pose, limit_only_to_length_);
	}

	MutationSetDesignOperation mut_set_op = MutationSetDesignOperation();

	for ( core::Size i_cdr = 1; i_cdr <= core::Size(local_ab_info->get_total_num_CDRs( true /* inlude cdr4 */)); ++i_cdr ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i_cdr );
		if ( sequences.count(cdr) == 0 || sequences[ cdr ].size() <= cutoff_ ) continue;
		TR << "applying profile set operation for " << local_ab_info->get_CDR_name(cdr) << std::endl;
		utility::vector1< std::map< core::Size, core::chemical::AA > > mutation_sets;
		for ( core::Size i = 1; i <= sequences[ cdr ].size(); ++i ) {
			std::string sequence = sequences[ cdr ][ i ].sequence;
			//TR << sequence << std::endl;
			std::map< core::Size, core::chemical::AA> mutation_set = design::transform_sequence_to_mutation_set(local_ab_info, pose, cdr, sequence);
			mutation_sets.push_back(mutation_set);
		}
		if ( mutation_sets.size() == 0 ) continue;

		//Apply the task for each CDR
		mut_set_op.set_mutation_sets(mutation_sets);
		mut_set_op.add_to_allowed_aas(keep_task_allowed_aas_);
		mut_set_op.include_native_aa(include_native_restype_);
		mut_set_op.set_picking_rounds(picking_rounds_);
		mut_set_op.apply( pose, task );
	}
}

core::pack::task::operation::TaskOperationOP
AddCDRProfileSetsOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new AddCDRProfileSetsOperation );
}

void AddCDRProfileSetsOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddCDRProfileSetsOperation::provide_xml_schema( xsd );
}

std::string AddCDRProfileSetsOperationCreator::keyname() const {
	return AddCDRProfileSetsOperation::keyname();
}

} //task_operations
} //antibody
} //protocols










