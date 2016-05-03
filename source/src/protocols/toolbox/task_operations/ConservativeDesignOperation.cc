// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/ConservativeDesignOperation.cc
/// @brief TaskOperation to allow only conservative mutations at designable residues
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/toolbox/task_operations/ConservativeDesignOperation.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/sql_utils.hh>
#include <cppdb/frontend.h>
#include <basic/database/open.hh>

#include <basic/Tracer.hh>
#include <boost/algorithm/string.hpp>
#include <utility/string_util.hh>

#include <basic/options/keys/task_operations.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.ConservativeDesignOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using utility::vector1;

/* AMW: Creator not defined
void ConservativeDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
ConservativeDesignOperation::provide_xml_schema( xsd );
}

void ConservativeDesignOperationCreator::keyname() const
{
return ConservativeDesignOperation::keyname();
}
*/

ConservativeDesignOperation::ConservativeDesignOperation() : core::pack::task::operation::TaskOperation(){
	set_defaults();
	set_data_source(basic::options::option [basic::options::OptionKeys::task_operations::cons_design_data_source]());
}

ConservativeDesignOperation::ConservativeDesignOperation(std::string data_source) : core::pack::task::operation::TaskOperation(){
	set_defaults();
	set_data_source(data_source);

}

ConservativeDesignOperation::~ConservativeDesignOperation(){}

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
ConservativeDesignOperation::set_data_source(std::string const data_source) {

	std::string original_data_source = data_source_;
	if ( data_source == "chothia_1976" ) {
		data_source_ = data_source;
	} else {
		std::string data_source_copy = data_source;
		boost::to_upper(data_source_copy);
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
}

ConservativeDesignOperation & ConservativeDesignOperation::operator =(ConservativeDesignOperation const & rhs){
	if ( this == &rhs ) return *this;
	core::pack::task::operation::TaskOperation::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

void
ConservativeDesignOperation::limit_to_positions(vector1<Size> const positions){
	positions_ = positions;
}

void
ConservativeDesignOperation::include_residue(core::Size resnum) {
	positions_.push_back(resnum);
}

void
ConservativeDesignOperation::clear_positions() {
	positions_.clear();
}

void
ConservativeDesignOperation::include_native_aa(bool const & setting) {
	include_native_aa_ = setting;
}

void
ConservativeDesignOperation::add_to_allowed_aas(bool const & setting) {
	add_to_allowed_aas_ = setting;
}

void
ConservativeDesignOperation::use_pose_sequence_as_native(core::pose::Pose const & pose){
	pose_sequence_ = pose.sequence();
}

void
ConservativeDesignOperation::set_native_sequence(std::string seq) {
	pose_sequence_ = seq;
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
ConservativeDesignOperation::apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const {

	vector1< bool > design_positions = task.designing_residues();
	for ( core::Size i = 1; i <=pose.total_residue(); ++i ) {


		//If its not a design position or the residue doesn't exist in set positions, keep going
		std::string seq = pose_sequence_;
		if ( pose_sequence_.empty() ) {
			seq = pose.sequence();
		}

		if ( ! design_positions[i] || (! positions_.empty() && std::find(positions_.begin(), positions_.end(), i) == positions_.end() ) ) { continue;}

		//Skip any unrecognized residues
		if ( ! core::chemical::oneletter_code_specifies_aa(seq[i-1]) ) { continue; };

		//task.nonconst_residue_task(i).get_original_residue()
		//TR << "Using conservative mutations for position: " << i << std::endl;
		core::chemical::AA native_amino = core::chemical::aa_from_oneletter_code(seq[i-1]);

		vector1<bool> allowed_aminos = conserved_mutations_[native_amino];

		//TR<< utility::to_string(allowed_aminos) << std::endl;
		if ( include_native_aa_ ) {
			allowed_aminos[native_amino] = true;
		}

		//Add the residues to the allowed list in task, or replace it.
		if ( add_to_allowed_aas_ ) {
			for ( core::Size aa_num = 1; aa_num <= 20; ++aa_num ) {
				core::chemical::AA amino = static_cast<core::chemical::AA>(aa_num);
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

// AMW: There is no parse_tag!
void ConservativeDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}


} //design
} //antibody
} //protocols
