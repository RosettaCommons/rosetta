// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/TCRmodel.cc
/// @brief Class for the tcr modeling protocol
/// @author Ragul Gowthaman (ragul@umd.edu)

// Project Headers
#include <protocols/tcr/TCRmodel.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/TCRloopRefine.hh>
#include <protocols/tcr/util.hh>
#include <protocols/tcr/grafting_util.hh>
#include <protocols/tcr/template_util.hh>
#include <protocols/tcr/modeling_util.hh>
#include <protocols/loops/loops_main.hh>
// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
// Utility Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/tcrmodel.OptionKeys.gen.hh>

///////////////////////////////////////////////////////////////////////////////


static basic::Tracer TR( "tcr.TCRmodel" );

namespace protocols {
namespace tcr {
using namespace utility;
using namespace basic::options;
using namespace basic::options::OptionKeys;

TCRmodel::TCRmodel( TCRseqInfoOP tcrinfo) {
	tsinfo = tcrinfo;
	init_from_options();
	set_default();
	setup_templates();
	make_model();
}

TCRmodel::TCRmodel( std::string const & aseq, std::string const & bseq ):
	utility::pointer::ReferenceCount()
{
	TCRseqInfoOP tcrinfo( new TCRseqInfo( aseq, bseq) );
	tsinfo = tcrinfo;
	init_from_options();
	set_default();
	setup_templates();
	make_model();
}

TCRmodelOP
TCRmodel::clone() const{
	return utility::pointer::make_shared< TCRmodel >( * this);
}

TCRmodel::~TCRmodel() = default;

void TCRmodel::init_from_options(){
	skip_modeling_ = option[OptionKeys::tcrmodel::skip_modeling];
	minimize_model_ = option[OptionKeys::tcrmodel::minimize_model];
	relax_model_ = option[OptionKeys::tcrmodel::relax_model];
	include_ab_templates_ = option[OptionKeys::tcrmodel::include_ab_templates];
	dump_templates_ = option[OptionKeys::tcrmodel::dump_templates];
	blastp_identity_cutoff_ = option[OptionKeys::tcrmodel::blastp_identity_cutoff];
	template_identity_cutoff_ = option[OptionKeys::tcrmodel::template_identity_cutoff];
	use_user_templates_ = option[OptionKeys::tcrmodel::use_user_templates];
	atmplt_.ori.tpdb = option[ OptionKeys::tcrmodel::alpha_orientation_template_pdb ];
	atmplt_.ori.tid = option[ OptionKeys::tcrmodel::alpha_orientation_template_id ];
	btmplt_.ori.tpdb = option[ OptionKeys::tcrmodel::beta_orientation_template_pdb ];
	btmplt_.ori.tid = option[ OptionKeys::tcrmodel::beta_orientation_template_id ];
	atmplt_.gm.tpdb = option[ OptionKeys::tcrmodel::alpha_germline_template_pdb ];
	atmplt_.gm.tid = option[ OptionKeys::tcrmodel::alpha_germline_template_id ];
	atmplt_.fr.tpdb = option[ OptionKeys::tcrmodel::alpha_framework_template_pdb ];
	atmplt_.fr.tid = option[ OptionKeys::tcrmodel::alpha_framework_template_id ];
	btmplt_.gm.tpdb = option[ OptionKeys::tcrmodel::beta_germline_template_pdb ];
	btmplt_.gm.tid = option[ OptionKeys::tcrmodel::beta_germline_template_id ];
	btmplt_.fr.tpdb = option[ OptionKeys::tcrmodel::beta_framework_template_pdb ];
	btmplt_.fr.tid = option[ OptionKeys::tcrmodel::beta_framework_template_id ];
	atmplt_.cdr1.tpdb = option[ OptionKeys::tcrmodel::alpha_cdr1_template_pdb ];
	atmplt_.cdr1.tid = option[ OptionKeys::tcrmodel::alpha_cdr1_template_id ];
	atmplt_.cdr2hv4.tpdb = option[ OptionKeys::tcrmodel::alpha_cdr2_template_pdb ];
	atmplt_.cdr2hv4.tid = option[ OptionKeys::tcrmodel::alpha_cdr2_template_id ];
	atmplt_.cdr3.tpdb = option[ OptionKeys::tcrmodel::alpha_cdr3_template_pdb ];
	atmplt_.cdr3.tid = option[ OptionKeys::tcrmodel::alpha_cdr3_template_id ];
	btmplt_.cdr1.tpdb = option[ OptionKeys::tcrmodel::beta_cdr1_template_pdb ];
	btmplt_.cdr1.tid = option[ OptionKeys::tcrmodel::beta_cdr1_template_id ];
	btmplt_.cdr2hv4.tpdb = option[ OptionKeys::tcrmodel::beta_cdr2_template_pdb ];
	btmplt_.cdr2hv4.tid = option[ OptionKeys::tcrmodel::beta_cdr2_template_id ];
	btmplt_.cdr3.tpdb = option[ OptionKeys::tcrmodel::beta_cdr3_template_pdb ];
	btmplt_.cdr3.tid = option[ OptionKeys::tcrmodel::beta_cdr3_template_id ];
	return;
}

void TCRmodel::set_default() {
	using namespace basic::options;
	std::string tcr_db_path = option[ OptionKeys::tcrmodel::tcr_template_db_path ];
	if ( tcr_db_path == "" ) {
		TR <<"Using TCR template from database/additional_protocol_data"<<std::endl;
		TR <<"Using TCR template from database/additional_protocol_data"<<std::endl;
		TR <<"You may have to clone or download the template database separately within the Rosetta/database/additional_protocol_data"<<std::endl;
		TR <<"https://github.com/RosettaCommons/additional_protocol_data"<<std::endl;
		tseq_db_ = basic::database::full_name("additional_protocol_data/tcr/seq/");
		tpdb_db_ = basic::database::full_name("additional_protocol_data/tcr/pdb/");
	} else {
		TR <<"Using TCR template databse provided by the user : "<<tcr_db_path<<std::endl;
		tseq_db_ = tcr_db_path + "/seq/";
		tpdb_db_ = tcr_db_path + "/pdb/";
	}
	if ( include_ab_templates_ ) {
		ab_db_path_ = option[ OptionKeys::tcrmodel::ab_template_db_path ];
		if ( ab_db_path_.empty() ) {
			utility_exit_with_message("Use the flag '-ab_template_db_path' to provide the antibody template database path");
		}
	}
	nter_overhang_ = 3;// = option[ num_nter_overhang_res ];
	cter_overhang_ = 3;//= option[ num_cter_overhang_res ];
	scorefxn_ = core::scoring::get_score_function();
	bool gma_tmplt = option[OptionKeys::tcrmodel::use_alpha_germline_templates];
	set_use_gma_templates(gma_tmplt);
	bool gmb_tmplt = option[OptionKeys::tcrmodel::use_beta_germline_templates];
	set_use_gmb_templates(gmb_tmplt);
	core::pose::Pose init_pose;
	init_pose.clear();
	tcr_model_ = core::pose::PoseOP( new core::pose::Pose(init_pose) );
	tcr_graft_model_ = core::pose::PoseOP( new core::pose::Pose(init_pose) );
	tcr_loop_model_ = core::pose::PoseOP( new core::pose::Pose(init_pose) );
	return;
}

void TCRmodel::setup_templates() {
	if ( !use_user_templates_ ) {
		initialize_template_db_files(tsinfo->atcr(), tsinfo->btcr(), atmplt_, btmplt_, tseq_db_);
		initialize_template_ignore_list(tsinfo->atcr(), tsinfo->btcr(), ignore_lists_, blastp_identity_cutoff_, tseq_db_);
		if ( !option[ OptionKeys::tcrmodel::use_alpha_germline_templates ].value() ) {
			if ( check_seq_match_from_multiple_input_db(tsinfo->atcr().gm, atmplt_.gm.tdb, ignore_lists_) ) {
				set_use_gma_templates(true);
			}
		}
		if ( !option[ OptionKeys::tcrmodel::use_beta_germline_templates ].value() ) {
			if ( check_seq_match_from_multiple_input_db(tsinfo->btcr().gm, btmplt_.gm.tdb, ignore_lists_) ) {
				set_use_gmb_templates(true);
			}
		}
	}
	//setup scoring scheme for scoring alignment
	utility::file::FileName pam30( basic::database::full_name("sequence/substitution_matrix/PAM30"));
	core::sequence::ScoringSchemeOP tcr_ss( new core::sequence::MatrixScoringScheme( -11, -1, pam30 ) );
	core ::Size startres;
	core::Size endres;
	setup_orientation_template(tsinfo->atcr(), tsinfo->btcr(), atmplt_.ori, btmplt_.ori, tsinfo->aaho_posi(), tsinfo->baho_posi(), tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
	if ( use_gma_templates_ ) {
		setup_fw_template(tsinfo->atcr().gm, tsinfo->aaho_posi(), atmplt_.gm, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
	} else {
		setup_fw_template(tsinfo->atcr().fr, tsinfo->aaho_posi(), atmplt_.fr, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
		if ( include_ab_templates_ ) {
			std::string l_chain_db = ab_db_path_ + "/ab_cdrseq_L1_" + std::to_string(tsinfo->atcr().cdr1.size()) + ".fasta";
			std::string h_chain_db = ab_db_path_ + "/ab_cdrseq_H1_" + std::to_string(tsinfo->atcr().cdr1.size()) + ".fasta";
			atmplt_.cdr1.tdb.push_back(l_chain_db);
			atmplt_.cdr1.tdb.push_back(h_chain_db);
		}
		startres = tsinfo->aaho_posi().cdr1.begin - nter_overhang_;
		endres = tsinfo->aaho_posi().cdr1.end + cter_overhang_;
		setup_cdr_template(tsinfo->atcr().cdr1, atmplt_.cdr1, startres, endres, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
		startres = tsinfo->aaho_posi().cdr2hv4.begin - nter_overhang_;
		endres = tsinfo->aaho_posi().cdr2hv4.end + cter_overhang_;
		setup_cdr_template(tsinfo->atcr().cdr2, atmplt_.cdr2hv4, startres, endres, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
	}
	if ( include_ab_templates_ ) {
		std::string l_chain_db = ab_db_path_ + "/ab_cdrseq_L3_" + std::to_string(tsinfo->atcr().cdr3.size()) + ".fasta";
		std::string h_chain_db = ab_db_path_ + "/ab_cdrseq_H3_" + std::to_string(tsinfo->atcr().cdr3.size()) + ".fasta";
		atmplt_.cdr3.tdb.push_back(l_chain_db);
		atmplt_.cdr3.tdb.push_back(h_chain_db);
	}
	startres = tsinfo->aaho_posi().cdr3.begin - nter_overhang_;
	endres = tsinfo->aaho_posi().cdr3.end + cter_overhang_;
	setup_cdr_template(tsinfo->atcr().cdr3, atmplt_.cdr3, startres, endres, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);

	if ( use_gmb_templates_ ) {
		setup_fw_template(tsinfo->btcr().gm, tsinfo->baho_posi(), btmplt_.gm, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
	} else {
		setup_fw_template(tsinfo->btcr().fr, tsinfo->baho_posi(), btmplt_.fr, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
		if ( include_ab_templates_ ) {
			std::string l_chain_db = ab_db_path_ + "/ab_cdrseq_L1_" + std::to_string(tsinfo->btcr().cdr1.size()) + ".fasta";
			std::string h_chain_db = ab_db_path_ + "/ab_cdrseq_H1_" + std::to_string(tsinfo->btcr().cdr1.size()) + ".fasta";
			btmplt_.cdr1.tdb.push_back(l_chain_db);
			btmplt_.cdr1.tdb.push_back(h_chain_db);
		}
		startres = tsinfo->baho_posi().cdr1.begin - nter_overhang_;
		endres = tsinfo->baho_posi().cdr1.end + cter_overhang_;
		setup_cdr_template(tsinfo->btcr().cdr1, btmplt_.cdr1, startres, endres, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
		startres = tsinfo->baho_posi().cdr2hv4.begin - nter_overhang_;
		endres = tsinfo->baho_posi().cdr2hv4.end + cter_overhang_;
		setup_cdr_template(tsinfo->btcr().cdr2, btmplt_.cdr2hv4, startres, endres, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);
	}
	if ( include_ab_templates_ ) {
		std::string l_chain_db = ab_db_path_ + "/ab_cdrseq_L3_" + std::to_string(tsinfo->btcr().cdr3.size()) + ".fasta";
		std::string h_chain_db = ab_db_path_ + "/ab_cdrseq_H3_" + std::to_string(tsinfo->btcr().cdr3.size()) + ".fasta";
		btmplt_.cdr3.tdb.push_back(l_chain_db);
		btmplt_.cdr3.tdb.push_back(h_chain_db);
	}
	startres = tsinfo->baho_posi().cdr3.begin - nter_overhang_;
	endres = tsinfo->baho_posi().cdr3.end + cter_overhang_;
	setup_cdr_template(tsinfo->btcr().cdr3, btmplt_.cdr3, startres, endres, tpdb_db_, ignore_lists_, template_identity_cutoff_, tcr_ss);

	if ( dump_templates_ ) {
		dump_templates(atmplt_, btmplt_, use_gma_templates_, use_gmb_templates_);
	}
	return;
}

void TCRmodel::build_graft_model() {
	core::pose::Pose scafolda;
	core::pose::Pose scafoldb;
	if ( use_gma_templates_ ) {
		orient_tcr_chain( *atmplt_.gm.tpiece, *atmplt_.ori.tpiece, tsinfo->aaho_posi());
		core::Size start;
		core::Size end;
		scafolda = *atmplt_.gm.tpiece;
		start = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr3.begin );
		end = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr3.end );
		graft_cdr_to_fw(scafolda, *atmplt_.cdr3.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
	} else {
		orient_tcr_chain( *atmplt_.fr.tpiece, *atmplt_.ori.tpiece, tsinfo->aaho_posi());
		core::Size start;
		core::Size end;
		scafolda = *atmplt_.fr.tpiece;
		start = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr1.begin );
		end = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr1.end );
		graft_cdr_to_fw(scafolda, *atmplt_.cdr1.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
		start = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr2hv4.begin );
		end = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr2hv4.end );
		graft_cdr_to_fw(scafolda, *atmplt_.cdr2hv4.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
		start = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr3.begin );
		end = scafolda.pdb_info()->pdb2pose( scafolda.pdb_info()->chain(1), tsinfo->aaho_posi().cdr3.end );
		graft_cdr_to_fw(scafolda, *atmplt_.cdr3.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
	}
	if ( use_gmb_templates_ ) {
		orient_tcr_chain( *btmplt_.gm.tpiece, *btmplt_.ori.tpiece, tsinfo->baho_posi());
		core::Size start;
		core::Size end;
		scafoldb = *btmplt_.gm.tpiece;
		start = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr3.begin );
		end = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr3.end );
		graft_cdr_to_fw(scafoldb, *btmplt_.cdr3.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
	} else {
		orient_tcr_chain( *btmplt_.fr.tpiece, *btmplt_.ori.tpiece, tsinfo->baho_posi());
		core::Size start;
		core::Size end;
		scafoldb = *btmplt_.fr.tpiece;
		start = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr1.begin );
		end = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr1.end );
		graft_cdr_to_fw(scafoldb, *btmplt_.cdr1.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
		start = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr2hv4.begin );
		end = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr2hv4.end );
		graft_cdr_to_fw(scafoldb, *btmplt_.cdr2hv4.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
		start = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr3.begin );
		end = scafoldb.pdb_info()->pdb2pose( scafoldb.pdb_info()->chain(1), tsinfo->baho_posi().cdr3.end );
		graft_cdr_to_fw(scafoldb, *btmplt_.cdr3.tpiece, start-1, end+1, cter_overhang_, nter_overhang_);
	}
	core::pose::Pose posecopya( scafolda );
	match_template_and_target_sequence(scafolda, tsinfo->atcr().truncdomain, posecopya);
	core::pose::Pose posecopyb( scafoldb );
	match_template_and_target_sequence(scafoldb, tsinfo->btcr().truncdomain, posecopyb);
	core::pose::Pose graft_pose = posecopya;
	graft_pose.append_pose_by_jump(posecopyb,1);
	tcr_graft_model_ = core::pose::PoseOP( new core::pose::Pose(graft_pose) );
	tcr_model_ = core::pose::PoseOP( new core::pose::Pose(graft_pose) );
	(*scorefxn_)(*tcr_graft_model_);
	TR << "graft pose score: " << tcr_graft_model_->energies().total_energies()[ core::scoring::total_score ] << std::endl;
	return;
}

void TCRmodel::make_model() {
	if ( skip_modeling_ ) return;
	build_graft_model();
	//setup loop model/refinement
	TCRloopRefineOP tcr_looprefine( new TCRloopRefine( tsinfo->atcr(), tsinfo->btcr(), tsinfo->aseq_posi(), tsinfo->bseq_posi()) );
	if ( tcr_looprefine->tcr_loop_model() ) {
		tcr_looprefine->apply( *tcr_model_ );
		tcr_loop_model_ = tcr_model_;
	} else {
		//refine_graft_model
		core::pose::Pose graft_model = *tcr_graft_model_;
		if ( relax_model_ ) {
			relax_tcr_model(graft_model);
			(*scorefxn_)(graft_model);
			TR << "relaxed_pose score: " << graft_model.energies().total_energies()[ core::scoring::total_score ] << std::endl;
		} else if ( minimize_model_ ) {
			minimize_tcr_model(graft_model);
			(*scorefxn_)(graft_model);
			TR << "minimized_pose score: " << graft_model.energies().total_energies()[ core::scoring::total_score ] << std::endl;
		} else { //just repack rotomers
			repack_tcr_model(graft_model);
			(*scorefxn_)(graft_model);
			TR << "repacked_pose score: " << graft_model.energies().total_energies()[ core::scoring::total_score ] << std::endl;
		}
		tcr_model_ = core::pose::PoseOP( new core::pose::Pose(graft_model) );
	}//end refine_graft_model
	(*scorefxn_)(*tcr_model_);
	TR << "tcrmodel score: " << tcr_model_->energies().total_energies()[ core::scoring::total_score ] << std::endl;
	return;
}

} // namespace tcr
} // namespace protocols
