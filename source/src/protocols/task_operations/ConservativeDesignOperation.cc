// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/ConservativeDesignOperation.cc
/// @brief TaskOperation to allow only conservative mutations at designable residues
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/task_operations/ConservativeDesignOperation.hh>
#include <protocols/task_operations/ConservativeDesignOperationCreator.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/database/sql_utils.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/task_operations.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <boost/algorithm/string.hpp>

#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <cppdb/frontend.h>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

#include <algorithm>

static basic::Tracer TR( "protocols.task_operations.ConservativeDesignOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using utility::vector1;

ConservativeDesignOperationCreator::ConservativeDesignOperationCreator(){}

ConservativeDesignOperationCreator::~ConservativeDesignOperationCreator(){}

void ConservativeDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConservativeDesignOperation::provide_xml_schema( xsd );
}

core::pack::task::operation::TaskOperationOP
ConservativeDesignOperationCreator::create_task_operation() const
{
	return utility::pointer::make_shared< ConservativeDesignOperation >();
}

std::string ConservativeDesignOperationCreator::keyname() const
{
	return ConservativeDesignOperation::keyname();
}


ConservativeDesignOperation::ConservativeDesignOperation() : core::pack::task::operation::TaskOperation(){
	set_defaults();
	set_data_source( basic::options::option[ basic::options::OptionKeys::task_operations::cons_design_data_source ]() );
}

ConservativeDesignOperation::ConservativeDesignOperation(std::string data_source) : core::pack::task::operation::TaskOperation(){
	set_defaults();
	set_data_source(data_source);
}

ConservativeDesignOperation::~ConservativeDesignOperation() = default;

void
ConservativeDesignOperation::set_defaults(){
	conserved_mutations_.resize(20);
	pose_sequence_.clear();
	add_to_allowed_aas_ = false;
	include_native_aa_ = true;
	data_source_ = "NA";
}

//void
//ConservativeDesignOperation::read_cmd_line_options() {

//}

void
ConservativeDesignOperation::set_data_source( std::string const & data_source ) {

	std::string original_data_source = data_source_;
	if ( data_source == "chothia_1976" ) {
		data_source_ = data_source;
	} else {
		std::string data_source_copy = data_source;
		boost::to_upper( data_source_copy );
		data_source_ = data_source_copy;
	}

	if ( original_data_source != data_source_ ) {
		load_data_from_db();
	}
}


core::pack::task::operation::TaskOperationOP
ConservativeDesignOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ConservativeDesignOperation( *this ) );
}

ConservativeDesignOperation::ConservativeDesignOperation(ConservativeDesignOperation const & rhs): core::pack::task::operation::TaskOperation(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

void
ConservativeDesignOperation::init_for_equal_operator_and_copy_constructor(ConservativeDesignOperation& lhs, ConservativeDesignOperation const & rhs){
	lhs.conserved_mutations_ = rhs.conserved_mutations_;
	lhs.add_to_allowed_aas_ = rhs.add_to_allowed_aas_;
	lhs.include_native_aa_ = rhs.include_native_aa_;
	lhs.positions_ = rhs.positions_;
	lhs.pose_sequence_ = rhs.pose_sequence_;
	lhs.data_source_ = rhs.data_source_;
	std::sort( lhs.positions_.begin(), lhs.positions_.end() );
}

ConservativeDesignOperation & ConservativeDesignOperation::operator =(ConservativeDesignOperation const & rhs){
	if ( this == &rhs ) return *this;
	core::pack::task::operation::TaskOperation::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

void
ConservativeDesignOperation::limit_to_positions( vector1< Size > const & positions){
	positions_ = positions;
	std::sort( positions_.begin(), positions_.end() );
}

void
ConservativeDesignOperation::include_residue( core::Size resnum ) {
	//keep positions_ sorted
	auto const iter = std::lower_bound( positions_.begin(), positions_.end(), resnum );
	positions_.insert( iter, resnum );
}

void
ConservativeDesignOperation::clear_positions() {
	positions_.clear();
}

void
ConservativeDesignOperation::use_pose_sequence_as_native( core::pose::Pose const & pose ){
	set_native_sequence( pose.sequence() );
}

void
ConservativeDesignOperation::load_data_from_db() {

	conserved_mutations_.resize(20);
	TR << "Loading conservative mutational data " <<  std::endl;
	std::string default_path = basic::database::full_name("/sequence/resinfo2.db");
	utility::sql_database::sessionOP session = basic::database::get_db_session(default_path);

	std::string base_statement =
		"SELECT \n"
		"\tid,\n"
		"\taa\n"
		"FROM \n"
		"\tconservation\n"
		"WHERE\n"
		"\tdata_source = ?;";

	cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, session));
	select_statement.bind(1, data_source_);

	cppdb::result amino_result(basic::database::safely_read_from_database(select_statement));

	while ( amino_result.next() ) {
		core::Size aa_index;
		std::string conserved_aminos;
		amino_result >> aa_index >> conserved_aminos;

		vector1<char> conserved_char(conserved_aminos.begin(), conserved_aminos.end());
		vector1<bool> allowed_aminos(20, false);

		for ( core::Size i = 1; i<= conserved_char.size(); ++i ) {
			core::chemical::AA  amino = core::chemical::aa_from_oneletter_code(conserved_char[ i ]);
			//TR <<"enabling " << i << " " << conserved_char[ i ] <<std::endl;
			allowed_aminos[amino] = true;
		}
		conserved_mutations_[aa_index] = allowed_aminos;
	}
}

void
ConservativeDesignOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const {

#ifndef NDEBUG
	for ( core::Size i = 1; i < positions_.size(); ++i ) {
		if ( positions_[ i ] > positions_[ i + 1 ] ) {
			utility_exit_with_message( "positions_ is not sorted" );
		}
	}
#endif
	//std::sort( positions_.begin(), positions_.end() );

	std::string const seq = ( pose_sequence_.empty() ? pose.sequence() : pose_sequence_ );

	utility::vector1< bool > const residues_allowed_by_selector =
		( residue_selector_ ? residue_selector_->apply( pose ) : utility::vector1< bool >( pose.size(), true ) );

	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		//If its not a design position or the residue doesn't exist in set positions, keep going
		if ( skip_resid( i, task, seq ) ) continue;

		if ( ! residues_allowed_by_selector[ i ] ) continue;

		core::chemical::AA native_amino = core::chemical::aa_from_oneletter_code(seq[i-1]);

		vector1< bool > allowed_aminos = conserved_mutations_[ native_amino ];

		//TR<< utility::to_string(allowed_aminos) << std::endl;
		if ( include_native_aa_ ) {
			allowed_aminos[ native_amino ] = true;
		}

		//Add the residues to the allowed list in task, or replace it.
		if ( add_to_allowed_aas_ ) {
			for ( core::Size aa_num = 1; aa_num <= 20; ++aa_num ) {
				auto amino = static_cast<core::chemical::AA>(aa_num);
				if ( allowed_aminos[aa_num] ) {
					task.nonconst_residue_task(i).allow_aa(amino);
				}
			}
		} else {
			task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aminos);
		}
		//task.show_residue_task(std::cout, i);
	}
}


void
ConservativeDesignOperation::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ){

	if ( tag->hasOption( "residue_selector" ) ) {
		std::string const & selector_name = tag->getOption< std::string >( "residue_selector" );
		try {
			residue_selector_ = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from ConservativeDesignOperation::parse_tag\n" + e.msg();
			throw CREATE_EXCEPTION( utility::excn::Exception,  error_message );
		}
	}

	if ( tag->hasOption( "add_to_allowed_aas" ) ) {
		add_to_allowed_aas( tag->getOption< bool >( "add_to_allowed_aas" ) );
	}

	if ( tag->hasOption( "include_native_aa" ) ) {
		include_native_aa( tag->getOption< bool >( "include_native_aa" ) );
	}

	if ( tag->hasOption( "data_source" ) ) {
		set_data_source( tag->getOption< std::string >( "data_source" ) );
	}

}


void ConservativeDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"residue_selector", xs_string,
		"Residue selector that indicates to which residues the operation will be applied.",
		"")
		+ XMLSchemaAttribute::attribute_w_default(
		"data_source", xs_string,
		"Set the source of the data used to define what is conservative. Options are: chothia_76 and the Blosum matrices from 30 to 100; designated as blosum30, 62, etc. Default is blosum62.  The higher the number, the more conservative the set of mutations (numbers are sequence identity cutoffs).",
		"blosum62")
		+ XMLSchemaAttribute::attribute_w_default(
		"add_to_allowed_aas", xsct_rosetta_bool,
		"Add to the allowed amino acids list instead of replacing it",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"include_native_aa", xsct_rosetta_bool,
		"Include native amino acid in the allowed_aas list",
		"true");

	task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"TaskOperation to allow only conservative mutations at designable residues");
}


bool ConservativeDesignOperation::skip_resid( core::Size const resid, core::pack::task::PackerTask const & task, std::string const & seq ) const {
	if ( ! task.being_designed( resid ) ) return true;
	if ( ! positions_.empty() ) { //if no positions are assigned, are we using all positions?
		if ( ! std::binary_search( positions_.begin(), positions_.end(), resid ) ) {
#ifndef NDEBUG
			for ( core::Size pos : positions_ ) {
				if ( pos == resid ) {
					TR << "positions_:";
					for ( core::Size pos2 : positions_ ) {
						TR << " " << pos2;
					}
					TR << std::endl;
					utility_exit_with_message( "binary_search could not find element: " + std::to_string( resid ) );
				}
			}
#endif
			return true;
		}
	}

	//Skip any unrecognized residues
	if ( ! core::chemical::oneletter_code_specifies_aa( seq[ resid - 1 ] ) ) {
		return true;
	}

	return false;
}

} //design
} //antibody
