// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/ConservativeDesignOperation.cc
/// @brief TaskOperation to allow only conservative mutations at designable residues
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/ConservativeDesignOperation.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/sql_utils.hh>

#include <cppdb/frontend.h>
#include <basic/database/open.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR("protocols.antibody.design.ConservativeDesignOperation");

namespace protocols {
namespace antibody {
namespace design {

ConservativeDesignOperation::ConservativeDesignOperation() : core::pack::task::operation::TaskOperation(){
	conserved_mutations_.resize(20);
	load_data_from_db();	
}

ConservativeDesignOperation::~ConservativeDesignOperation(){}

core::pack::task::operation::TaskOperationOP
ConservativeDesignOperation::clone() const
{
	return new ConservativeDesignOperation( *this );
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
ConservativeDesignOperation::load_data_from_db() {
	TR << "Loading conservative mutational data " <<  std::endl;
	std::string default_path = basic::database::full_name("/sequence/resinfo.db");
	utility::sql_database::sessionOP session = basic::database::get_db_session(default_path);
	
	std::string base_statement = 
		"SELECT \n"
		"	id,\n"
		"	conserved\n"
		"FROM \n"
		"	resinfo;";
	
	cppdb::statement select_statement(basic::database::safely_prepare_statement(base_statement, session));
	cppdb::result amino_result(basic::database::safely_read_from_database(select_statement));
	
	while(amino_result.next()){
		core::Size aa_index;
		std::string conserved_aminos;
		amino_result >> aa_index >> conserved_aminos;
		
		vector1<char> conserved_char(conserved_aminos.begin(), conserved_aminos.end());
		vector1<bool> allowed_aminos(20, false);
		
		conserved_mutations_[aa_index] = allowed_aminos;
	}
}

void
ConservativeDesignOperation::apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const {
	
	vector1< bool > design_positions = task.designing_residues();
	for (core::Size i = 1; i <=pose.total_residue(); ++i){
		
		
		//If its not a design position or the residue doesn't exist in set positions, keep going
		std::string seq = pose_sequence_; //dangerous.
		if (pose_sequence_.empty()){
			seq = pose.sequence();
		}
		
		if (! design_positions[i] || (! positions_.empty() && std::find(positions_.begin(), positions_.end(), i) == positions_.end() ) ){continue;}
		
		//Skip any unrecognized residues
		if (! core::chemical::oneletter_code_specifies_aa(seq[i-1])){ continue; };
		
		//task.nonconst_residue_task(i).get_original_residue()
		//TR << "Using conservative mutations for position: " << i << std::endl;
		core::chemical::AA native_amino = core::chemical::aa_from_oneletter_code(seq[i-1]);
		
		vector1<bool> allowed_aminos = conserved_mutations_[native_amino];
		
		if (include_native_aa_){
			allowed_aminos[native_amino] = true;	
		}
		
		//Add the residues to the allowed list in task, or replace it.
		if (add_to_allowed_aas_){
			for (core::Size aa_num; aa_num <= 20; ++aa_num){
				core::chemical::AA amino = static_cast<core::chemical::AA>(aa_num);
				if (allowed_aminos[aa_num]){
					task.nonconst_residue_task(i).allow_aa(amino);
				}
			}
		}
		else{
			task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aminos);
		}
		//task.show_residue_task(std::cout, i);
	}
}



} //design
} //antibody
} //protocols
