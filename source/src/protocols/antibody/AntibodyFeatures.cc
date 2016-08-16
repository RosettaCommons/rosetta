// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyFeatures.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/AntibodyFeatures.hh>
#include <protocols/antibody/metrics.hh>

#include <protocols/features/InterfaceFeatures.hh>
#include <protocols/features/FeaturesReporter.hh>
//#include <protocols/analysis/InterfaceAnalyzerMover.hh>

//Core Headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/util.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/id/AtomID.hh>

//Protocols Headers
#include <protocols/features/util.hh>

//Basic Headers
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/sql_utils.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <cppdb/frontend.h>


#include <boost/foreach.hpp>
#include <algorithm>

#define boost_foreach BOOST_FOREACH

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.AntibodyFeatures" );

namespace protocols {
namespace antibody {
using namespace protocols::features;
using namespace protocols::analysis;
using namespace core::scoring;

AntibodyFeatures::AntibodyFeatures():
	InterfaceFeatures()
{
	ab_info_ = NULL;
	skip_antigen_reports_ = false;
	include_proto_cdr4_ = false;
}

AntibodyFeatures::AntibodyFeatures(AntibodyInfoOP ab_info):
	InterfaceFeatures()
{
	ab_info_ = ab_info;
	skip_antigen_reports_ = false;
	include_proto_cdr4_ = false;
}

AntibodyFeatures::AntibodyFeatures(AntibodyInfoOP ab_info, ScoreFunctionCOP scorefxn):
	InterfaceFeatures(scorefxn)
{
	ab_info_ = ab_info;
	skip_antigen_reports_ = false;
	include_proto_cdr4_ = false;
}

void
AntibodyFeatures::set_ab_info(AntibodyInfoOP ab_info){
	ab_info_ = ab_info;
}

void
AntibodyFeatures::set_interface_chains(utility::vector1<std::string> const & intermediate_interfaces){
	intermediate_interfaces_ = intermediate_interfaces;
}

std::string
AntibodyFeatures::type_name() const {
	return "AntibodyFeatures";
}

void
AntibodyFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const {
	InterfaceFeatures::write_interface_schema_to_db(db_session);
	InterfaceFeatures::write_interface_side_schema_to_db(db_session);
	InterfaceFeatures::write_interface_residues_schema_to_db(db_session);

	write_ab_metrics_schema_to_db(db_session);
	write_cdr_metrics_schema_to_db(db_session);
	write_ab_H3_kink_metrics_schema_to_db(db_session);

	//write_cdr_definitions_schema_to_db(db_session);
	write_cdr_residue_schema_to_db(db_session);
}

core::Size
AntibodyFeatures::report_features(
	core::pose::Pose const & pose,
	const utility::vector1<bool>& relevant_residues,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session) {

	if ( ! ab_info_ ) {
		ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose) );
	}

	std::map<std::string, std::string > db_interfaces;


	if ( intermediate_interfaces_.size() == 0 ) {

		if ( ! ab_info_->is_camelid() ) {
			db_interfaces["L_H"] = "L_H";
		}

		if ( ab_info_->antigen_present() && ! skip_antigen_reports_ ) {
			utility::vector1<char> antigen_chains = ab_info_->get_antigen_chains();
			std::string antigen(antigen_chains.begin(), antigen_chains.end());

			db_interfaces["H_"+antigen] = "H_A";

			if ( ! ab_info_->is_camelid() ) {
				db_interfaces["LH_"+antigen] = "LH_A";
				db_interfaces["L_"+antigen] = "L_A";
			}
		}
	} else {
		//Could use Enums for this:
		for ( core::Size i = 1; i <= intermediate_interfaces_.size(); ++i ) {
			if ( intermediate_interfaces_[i] == "L_H" && ! ab_info_->is_camelid() ) {
				db_interfaces["L_H"] = "L_H";
			}

			if ( ab_info_->antigen_present() ) {
				utility::vector1<char> antigen_chains = ab_info_->get_antigen_chains();
				std::string antigen(antigen_chains.begin(), antigen_chains.end());

				if ( intermediate_interfaces_[i] == "L_A" && ! skip_antigen_reports_ && ! ab_info_->is_camelid() ) {
					db_interfaces["L_"+antigen] = "L_A";
				} else if ( intermediate_interfaces_[i] == "H_A" && ! skip_antigen_reports_ ) {
					db_interfaces["H_"+antigen] = "H_A";
				} else if ( intermediate_interfaces_[i] == "LH_A" && ! skip_antigen_reports_ && !ab_info_->is_camelid() ) {
					db_interfaces["LH_"+antigen] = "LH_A";
				} else if ( intermediate_interfaces_[i] == "LH_A" && ! skip_antigen_reports_ && ab_info_->is_camelid() ) {
					TR << "Replacing LH_A with H_A for camelid antibody.  Your welcome." << std::endl;
					db_interfaces["H_"+antigen] ="H_A";
				} else {
					//else if (intermediate_interfaces_[i] == "NONE"){
					// TR << "Not reporting any InterfaceFeatures." << std::endl;
					//}
					TR << "Skipping: " << intermediate_interfaces_[i] <<" please use InterfaceFeatures for general interfaces" << std::endl;
					continue;
				}
			}
		}
	}

	//Find which per residue data to use from IAM.  We have camelid and regular LH antibody.  an LA antibody wouldn't be kosher in AntibodyInfo right now anyway.
	//If its not in the interface, and not set to analyze, then we don't do anything (For speed issues?).
	std::string match_interface_type;
	if ( ab_info_->antigen_present() ) {
		if ( ab_info_->is_camelid() ) {
			match_interface_type =  "H_A";
		} else {
			match_interface_type = "LH_A";
		}

	}


	bool have_interface_data = false;
	typedef std::map<std::string, std::string> map_type;
	boost_foreach ( const map_type::value_type& interface_pair, db_interfaces ) {
		report_all_interface_features(pose, relevant_residues, struct_id, db_session, interface_pair.first, interface_pair.second);

		if ( interface_pair.second == match_interface_type && ! skip_antigen_reports_ ) {
			interface_data_ = interface_analyzer_->get_all_data();
			interface_data_res_ = interface_analyzer_->get_all_per_residue_data();
			have_interface_data = true;
		}

	}

	//We need zero values for interface_data and interface_data_res or we need to get data
	if ( ! have_interface_data ) {
		protocols::analysis::InterfaceAnalyzerMoverOP IAM( new protocols::analysis::InterfaceAnalyzerMover(match_interface_type, true, scorefxn_, compute_packstat_, pack_together_, pack_separated_) );
		if ( ! skip_antigen_reports_ ) {
			IAM->init_on_new_input(pose); //Zeros and initializes all values so we can use this as zero.
			interface_data_ = IAM->get_all_data();
			interface_data_res_ = IAM->get_all_per_residue_data();
		} else {
			IAM->apply_const(pose);
			interface_data_ = IAM->get_all_data();
			interface_data_res_ = IAM->get_all_per_residue_data();
		}

	}


	utility::vector1<bool> cdr_residues(pose.total_residue(), false);
	utility::vector1<bool> antigen_residues(pose.total_residue(), false);

	//Setup CDR residues:
	for ( core::SSize i = 1; i <= ab_info_->get_total_num_CDRs(include_proto_cdr4_); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		//TR << "Setting up CDR residues: " << ab_info_->get_CDR_name(cdr) << std::endl;

		for ( core::Size i = ab_info_->get_CDR_start(cdr, pose); i <= ab_info_->get_CDR_end(cdr, pose); ++i ) {
			//TR << "Residues: " << i << std::endl;
			cdr_residues[i] = true;
		}
	}

	//Setup antigen residues:
	vector1<core::Size> ag_chains = ab_info_->get_antigen_chain_ids(pose);
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( std::find(ag_chains.begin(), ag_chains.end(), pose.residue(i).chain()) != ag_chains.end() ) {
			antigen_residues[i] = true;
		}
	}

	//Get antigen - antibody contacts
	calculate_residue_atomic_contacts(pose, cdr_residues, antigen_residues);

	paratope_sasa_ = paratope_sasa(pose, *ab_info_, include_proto_cdr4_);
	paratope_charge_ = paratope_charge(pose, *ab_info_, include_proto_cdr4_);

	//Antibody-specific tables:



	for ( core::SSize i = 1; i <= ab_info_->get_total_num_CDRs(include_proto_cdr4_); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		report_cdr_metrics_features(pose, struct_id, db_session, cdr);
		report_cdr_residue_features(pose, struct_id, db_session, cdr, relevant_residues);

	}
	report_ab_H3_kink_metrics_features(pose, struct_id, db_session);
	TR << "reported kink metrics" << std::endl;
	report_ab_metrics_features(pose, struct_id, db_session);
	TR << "reported ab metrics " << std::endl;

	return 0;
}


void
AntibodyFeatures::write_ab_metrics_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;



	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column num_scheme("numbering_scheme", DbDataTypeOP( new DbText() ));
	Column cdr_def("cdr_definition", DbDataTypeOP( new DbText() ));
	Column cdr_res("cdr_residues", DbDataTypeOP( new DbInteger() ));
	Column vl_vh_pack("VL_VH_packing_angle", DbDataTypeOP( new DbReal() ));
	Column vl_vh_distance("VL_VH_distance", DbDataTypeOP( new DbReal() ));
	Column vl_vh_open("VL_VH_opening_angle", DbDataTypeOP( new DbReal() ));
	Column vl_vh_opp("VL_VH_opposite_opening_angle", DbDataTypeOP( new DbReal() ));
	Column antigen_present("antigen_present", DbDataTypeOP( new DbInteger() ));
	Column antigen_chains("antigen_chains", DbDataTypeOP( new DbText() ));
	Column antigen_charge("net_charge", DbDataTypeOP( new DbInteger() ));
	Column paratope_charge("paratope_charge", DbDataTypeOP( new DbInteger() ));
	Column paratope_SASA("paratope_SASA", DbDataTypeOP( new DbReal() ));
	Column paratope_hSASA("paratope_hSASA", DbDataTypeOP( new DbReal() ));
	Column paratope_pSASA("paratope_pSASA", DbDataTypeOP( new DbReal() ));
	Column is_camelid("is_camelid", DbDataTypeOP( new DbInteger() ));


	PrimaryKey primary_key(struct_id);
	ForeignKey foreign_key(struct_id, "structures", "struct_id", true);



	Schema table("ab_metrics", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(num_scheme);
	table.add_column(cdr_def);
	table.add_column(cdr_res);
	table.add_column(antigen_present);
	table.add_column(antigen_chains);
	table.add_column(antigen_charge);

	table.add_column(paratope_charge);
	table.add_column(paratope_SASA);
	table.add_column(paratope_hSASA);
	table.add_column(paratope_pSASA);
	table.add_column(vl_vh_pack);
	table.add_column(vl_vh_distance);
	table.add_column(vl_vh_open);
	table.add_column(vl_vh_opp);
	table.add_column(is_camelid);

	table.write(db_session);

}

void
AntibodyFeatures::report_ab_metrics_features(
	core::pose::Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	std::string stmt_string = "INSERT INTO ab_metrics("
		"struct_id,"
		"numbering_scheme,"
		"cdr_definition,"
		"cdr_residues,"
		"antigen_present,"
		"antigen_chains,"
		"net_charge,"
		"paratope_charge,"
		"paratope_SASA,"
		"paratope_hSASA,"
		"paratope_pSASA,"
		"VL_VH_opening_angle,"
		"VL_VH_distance,"
		"VL_VH_opposite_opening_angle,"
		"VL_VH_packing_angle,"
		"is_camelid"
		") VALUES "+get_question_mark_string(16);

	cppdb::statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));

	//Get CDR residue total - don't think this would be useful anywhere
	core::Size cdr_residues = 0;
	for ( core::SSize i = 1; i <= ab_info_->get_total_num_CDRs(include_proto_cdr4_); ++i ) {
		CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
		cdr_residues = cdr_residues + ab_info_->get_CDR_length(cdr_name, pose);
	}

	//Get packing angle info
	utility::vector1<core::Real> vl_vh_coords = vl_vh_orientation_coords(pose, *ab_info_);

	//Get chains
	utility::vector1<char> antigen_chains = ab_info_->get_antigen_chains();
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	if ( ! ab_info_->antigen_present() ) {
		antigen = "NA";
	}

	//TR << "camelid" <<  int(ab_info_->is_camelid()) << std::endl;
	core::Size i = 0;
	stmt.bind(i+=1, struct_id);

	stmt.bind(i+=1, ab_info_->get_current_AntibodyNumberingScheme());
	stmt.bind(i+=1, ab_info_->get_current_CDRDefinition());
	stmt.bind(i+=1, cdr_residues);
	stmt.bind(i+=1, int(ab_info_->antigen_present()));
	stmt.bind(i+=1, antigen);
	stmt.bind(i+=1, pose_charge(pose));
	stmt.bind(i+=1, paratope_charge_.paratope);
	stmt.bind(i+=1, paratope_sasa_.first.paratope);
	stmt.bind(i+=1, paratope_sasa_.second.paratope);
	//stmt.bind(i+=1, paratope_sasa_no_h_.first.paratope);
	//stmt.bind(i+=1, paratope_sasa_no_h_.second.paratope);
	stmt.bind(i+=1, paratope_sasa_.first.paratope - paratope_sasa_.second.paratope);
	stmt.bind(i+=1, vl_vh_coords[4]);
	stmt.bind(i+=1, vl_vh_coords[1]);
	stmt.bind(i+=1, vl_vh_coords[2]);
	stmt.bind(i+=1, vl_vh_coords[3]);
	stmt.bind(i+=1, int(ab_info_->is_camelid()));

	basic::database::safely_write_to_database(stmt);

}

void
AntibodyFeatures::write_cdr_metrics_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;


	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column cdr("CDR", DbDataTypeOP( new DbText() ));
	Column length("length", DbDataTypeOP( new DbInteger() ));
	Column start("start", DbDataTypeOP( new DbInteger() ));
	Column end("end", DbDataTypeOP( new DbInteger() ));
	Column ag_contact("ag_ab_contacts_total", DbDataTypeOP( new DbInteger() ));
	Column ag_contact_res("ag_ab_contacts_nres", DbDataTypeOP( new DbInteger() ));
	Column ag_ab_dsasa("ag_ab_dSASA", DbDataTypeOP( new DbReal() ));
	Column ag_ab_dsasa_sc("ag_ab_dSASA_sc", DbDataTypeOP( new DbReal() ));
	Column ag_ab_dhsasa("ag_ab_dhSASA", DbDataTypeOP( new DbReal() ));
	Column ag_ab_dhsasa_sc("ag_ab_dhSASA_sc", DbDataTypeOP( new DbReal() ));
	Column ag_ab_dhsasa_rel("ag_ab_dhSASA_rel_by_charge", DbDataTypeOP( new DbReal() ));

	Column SASA("SASA", DbDataTypeOP( new DbReal() ));
	Column energy("energy", DbDataTypeOP( new DbReal() ));
	Column charge("charge", DbDataTypeOP( new DbInteger() ));
	Column dG("ag_ab_dG", DbDataTypeOP( new DbReal() ));
	Column anchor_dis("anchor_CN_distance", DbDataTypeOP( new DbReal() ));
	Column aromatic_nres("aromatic_nres", DbDataTypeOP( new DbReal() ));

	Columns primary_keys;
	primary_keys.push_back(struct_id);
	primary_keys.push_back(cdr);
	PrimaryKey primary(primary_keys);
	ForeignKey foreign_key(struct_id, "structures", "struct_id", true);

	Schema table("cdr_metrics", primary);
	table.add_foreign_key(foreign_key);
	table.add_column(length);
	table.add_column(start);
	table.add_column(end);
	table.add_column(ag_contact);
	table.add_column(ag_contact_res);
	table.add_column(ag_ab_dsasa);
	table.add_column(ag_ab_dsasa_sc);
	table.add_column(ag_ab_dhsasa);
	table.add_column(ag_ab_dhsasa_sc);
	table.add_column(ag_ab_dhsasa_rel);
	table.add_column(dG);
	table.add_column(SASA);
	table.add_column(charge);
	table.add_column(energy);
	table.add_column(anchor_dis);
	table.add_column(aromatic_nres);
	table.write(db_session);
	//CDR_contacts_to_other_cdrs->hbonds->R
}

void
AntibodyFeatures::report_cdr_metrics_features(
	core::pose::Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	CDRNameEnum const & cdr)
{


	std::string stmt_string = "INSERT INTO cdr_metrics("
		"struct_id,"
		"CDR,"
		"length,"
		"start,"
		"end,"
		"ag_ab_contacts_total,"
		"ag_ab_contacts_nres,"
		"ag_ab_dSASA,"
		"ag_ab_dSASA_sc,"
		"ag_ab_dhSASA,"
		"ag_ab_dhSASA_sc,"
		"ag_ab_dhSASA_rel_by_charge,"
		"ag_ab_dG,"
		"SASA,"
		"charge,"
		"energy,"
		"anchor_CN_distance,"
		"aromatic_nres"
		") VALUES " + get_question_mark_string(18);
	cppdb::statement stmt = basic::database::safely_prepare_statement(stmt_string, db_session);

	std::string cdr_name = ab_info_->get_CDR_name(cdr);
	//TR << "cdr_name: "<<cdr_name<<std::endl;
	//TR << struct_id << std::endl;

	core::Size i = 0;

	stmt.bind(i+=1, struct_id);
	stmt.bind(i+=1, cdr_name);
	stmt.bind(i+=1, ab_info_->get_CDR_length(cdr, pose));
	stmt.bind(i+=1, ab_info_->get_CDR_start(cdr, pose));
	stmt.bind(i+=1, ab_info_->get_CDR_end(cdr, pose));
	stmt.bind(i+=1, calculate_cdr_contacts_total(pose, cdr));
	stmt.bind(i+=1, calculate_cdr_contacts_nres(pose, cdr));
	stmt.bind(i+=1, calculate_cdr_totals(cdr, pose, interface_data_res_.dSASA));
	stmt.bind(i+=1, calculate_cdr_totals(cdr, pose, interface_data_res_.dSASA_sc));
	stmt.bind(i+=1, calculate_cdr_totals(cdr, pose, interface_data_res_.dhSASA));
	stmt.bind(i+=1, calculate_cdr_totals(cdr, pose, interface_data_res_.dhSASA_sc));
	stmt.bind(i+=1, calculate_cdr_totals(cdr, pose, interface_data_res_.dhSASA_rel_by_charge));
	stmt.bind(i+=1, calculate_cdr_totals(cdr, pose, interface_data_res_.dG));
	stmt.bind(i+=1, paratope_sasa_.first.cdr[cdr]);
	stmt.bind(i+=1, paratope_charge_.cdr[cdr]);
	stmt.bind(i+=1, cdr_energy(pose, ab_info_, scorefxn_, cdr));
	stmt.bind(i+=1, cdr_CN_anchor_distance(pose, ab_info_, cdr));
	stmt.bind(i+=1, calculate_cdr_aromatic_nres(pose, cdr));
	//aromatic_fraction??
	basic::database::safely_write_to_database(stmt);

}

void
AntibodyFeatures::write_cdr_residue_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resnum("resNum", DbDataTypeOP( new DbInteger() ));
	Column cdr("CDR", DbDataTypeOP( new DbText() ));
	Column cdr_position("position", DbDataTypeOP( new DbInteger() ));
	Column ag_contacts("ag_contacts", DbDataTypeOP( new DbInteger() ));

	Columns primary_keys;
	primary_keys.push_back(struct_id);
	primary_keys.push_back(resnum);

	PrimaryKey primary_key(primary_keys);
	ForeignKey struct_foreign_key(struct_id, "structures", "struct_id", true);

	vector1<Column> res_keys;
	res_keys.push_back(struct_id);
	res_keys.push_back(resnum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey res_foreign_keys(res_keys, "residues", reference_columns, true);

	Schema table("cdr_residues", primary_key);
	table.add_foreign_key(res_foreign_keys);
	table.add_foreign_key(struct_foreign_key);
	table.add_column(cdr);
	table.add_column(cdr_position);
	table.add_column(ag_contacts);

	table.write(db_session);


}

void
AntibodyFeatures::report_cdr_residue_features(
	core::pose::Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	CDRNameEnum const & cdr,
	utility::vector1<bool> const & relevant_residues)
{
	for ( core::Size i = 1; i <= ab_info_->get_CDR_length(cdr, pose); ++i ) {
		core::Size resnum = ab_info_->get_CDR_start(cdr, pose) + i - 1;
		if ( relevant_residues[resnum] ) {
			report_cdr_residue_features_row(pose, struct_id, db_session, cdr, resnum, i);
		}
	}

}

void
AntibodyFeatures::report_cdr_residue_features_row(
	core::pose::Pose const & ,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	CDRNameEnum const & cdr,
	core::Size resnum,
	core::Size position)
{
	std::string stmt_string = "INSERT INTO cdr_residues("
		"struct_id,"
		"resNum,"
		"CDR,"
		"position,"
		"ag_contacts"
		") VALUES" + get_question_mark_string(5);

	cppdb::statement stmt = basic::database::safely_prepare_statement(stmt_string, db_session);

	core::Size i = 0;
	stmt.bind(i+=1, struct_id);
	stmt.bind(i+=1, resnum);
	stmt.bind(i+=1, ab_info_->get_CDR_name(cdr));
	stmt.bind(i+=1, position);
	stmt.bind(i+=1, ag_ab_atomic_contacts_[resnum]);

	basic::database::safely_write_to_database(stmt);

}

void
AntibodyFeatures::write_ab_H3_kink_metrics_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column H3_kink("kink_type", DbDataTypeOP( new DbText() ));
	Column kink_begin("begin", DbDataTypeOP( new DbInteger() ));
	Column kink_end("end", DbDataTypeOP( new DbInteger() ));
	Column anion_res("anion_res", DbDataTypeOP( new DbInteger() ));
	Column cation_res("cation_res", DbDataTypeOP( new DbInteger() ));
	Column kink_Rd("RD_Hbond_dis", DbDataTypeOP( new DbReal() ));
	Column kink_bb("bb_Hbond_dis", DbDataTypeOP( new DbReal() ));
	Column kink_trp("Trp_Hbond_dis", DbDataTypeOP( new DbReal() ));
	Column qdis("qdis", DbDataTypeOP( new DbReal() ));
	Column qdih("qdih", DbDataTypeOP( new DbReal() ));

	PrimaryKey primary(struct_id);
	ForeignKey foreign_key(struct_id, "structures", "struct_id", true);

	Schema table("ab_H3_kink_metrics", primary);
	table.add_foreign_key(foreign_key);
	table.add_column(H3_kink);
	table.add_column(kink_begin);
	table.add_column(kink_end);
	table.add_column(anion_res);
	table.add_column(cation_res);
	table.add_column(kink_Rd);
	table.add_column(kink_bb);
	table.add_column(kink_trp);
	table.add_column(qdis);
	table.add_column(qdih);

	table.write(db_session);

}

void
AntibodyFeatures::report_ab_H3_kink_metrics_features(
	const core::pose::Pose& pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	std::string stmt_string = "INSERT INTO ab_H3_kink_metrics("
		"struct_id,"
		"kink_type,"
		"begin,"
		"end,"
		"anion_res,"
		"cation_res,"
		"RD_Hbond_dis,"
		"bb_Hbond_dis,"
		"Trp_Hbond_dis,"
		"qdis,"
		"qdih"
		") VALUES "+get_question_mark_string(11);
	cppdb::statement stmt = basic::database::safely_prepare_statement(stmt_string, db_session);

	std::pair<core::Real, core::Real > k_dih = kink_dihedral(pose, *ab_info_);

	core::Size i = 0;
	stmt.bind(i+=1, struct_id);
	stmt.bind(i+=1, ab_info_->get_H3_kink_type_name());
	stmt.bind(i+=1, ab_info_->kink_begin(pose));
	stmt.bind(i+=1, ab_info_->kink_end(pose));
	stmt.bind(i+=1, ab_info_->kink_anion_residue(pose));
	stmt.bind(i+=1, ab_info_->kink_cation_residue(pose));
	stmt.bind(i+=1, kink_RD_Hbond(pose, *ab_info_));
	stmt.bind(i+=1, kink_bb_Hbond(pose, *ab_info_));
	stmt.bind(i+=1, kink_Trp_Hbond(pose, *ab_info_));
	stmt.bind(i+=1, k_dih.first);
	stmt.bind(i+=1, k_dih.second);

	basic::database::safely_write_to_database(stmt);
}


/* Not sure if this is useful.
void
AntibodyFeatures::write_cdr_definitions_schema_to_db(utility::sql_database::sessionOP db_session) const {
using namespace basic::database::schema_generator;
//struct_id
//cdr
//length in different definitions
//start and end in different definitions
}
*/


void
AntibodyFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & /*data*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose)
{

	pack_separated_ = tag->getOption<bool>("pack_separated", true);
	pack_together_ = tag->getOption<bool>("pack_together", false);
	dSASA_cutoff_ = tag->getOption<core::Real>("dSASA_cutoff", 100);
	compute_packstat_ = tag->getOption<bool>("compute_packstat", true);
	skip_antigen_reports_ = tag->getOption<bool>("skip_all_antigen_analysis", false); //By default if antigen is present, we will report some data on it.

	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_name(tag->getOption<std::string>("scorefxn"));
		scorefxn_ = data.get_ptr<core::scoring::ScoreFunction>("scorefxns", scorefxn_name);
	}

	if ( tag->hasOption("interfaces") && tag->hasOption("interface") ) {
		utility_exit_with_message("Cannot specify both interface and interfaces option for InterfaceFeatures Reporter");
	}

	if ( tag->hasOption("interface") ) {
		std::string interface = tag->getOption<std::string>("interface");
		if ( interface.find(",") != std::string::npos ) { utility_exit_with_message("Only one interface should be specified using the interface option");}
		intermediate_interfaces_.push_back(interface);
	}

	if ( tag->hasOption("interfaces") ) {
		std::string interfaces = tag->getOption<std::string>("interfaces");
		intermediate_interfaces_ = utility::string_split_multi_delim(interfaces, ":,'`~+*&|;.");//Why not?
	}


	if ( tag->hasOption("cdr_definition") && tag->hasOption("numbering_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();

		CDRDefinitionEnum definition = manager.cdr_definition_string_to_enum(tag->getOption<std::string>("cdr_definition"));
		AntibodyNumberingSchemeEnum scheme = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));

		ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose, scheme, definition) );

	} else if ( tag->hasOption("cdr_definition") || tag->hasOption("numbering_scheme") ) {
		TR <<"Please pass both cdr_definition and numbering_scheme.  These can also be set via cmd line options of the same name." << std::endl;

	}

	include_proto_cdr4_ = tag->getOption("include_proto_cdr4", include_proto_cdr4_);

}

core::Real
AntibodyFeatures::calculate_cdr_totals(const CDRNameEnum cdr, const core::pose::Pose & pose, const utility::vector1<core::Real> & data) const {
	core::Real result= 0.0;
	core::Size cdr_start = ab_info_->get_CDR_start(cdr, pose);
	core::Size cdr_end = ab_info_->get_CDR_end(cdr, pose);

	for ( core::Size i = cdr_start; i <= cdr_end; ++i ) {
		result = result + data[i];
	}
	return result;
}

core::Size
AntibodyFeatures::calculate_cdr_aromatic_nres(const core::pose::Pose & pose, const CDRNameEnum cdr) {
	core::Size nres= 0;
	core::Size cdr_start = ab_info_->get_CDR_start(cdr, pose);
	core::Size cdr_end = ab_info_->get_CDR_end(cdr, pose);

	for ( core::Size i = cdr_start; i <= cdr_end; ++i ) {
		if ( pose.residue(i).type().is_aromatic() ) {
			nres+=1;
		}
	}
	return nres;
}

void
AntibodyFeatures::calculate_residue_atomic_contacts(const core::pose::Pose& pose, const utility::vector1<bool> & residues_to_match, const utility::vector1<bool> & antigen_residues) {
	using namespace core::conformation;

	//We want to do this for atom-atom neighbors only,
	// the index of each vertex would need to be the AtomID to tie it back to antigen-cdr contact.
	// In order to do this we map index to AtomId and then trace it back when needed.

	core::Real dist_cutoff = 5;
	core::Size contact_cutoff = 5;

	ag_ab_atomic_contacts_.clear();
	ag_ab_atomic_contacts_.resize(pose.total_residue(), 0);


	//semi Slow N2 workaround till graph can be debugged:

	//CDR Residues:
	TR << "Calculating antigen contacts" << std::endl;
	for ( core::Size ci = 1; ci <= pose.total_residue(); ++ci ) {
		if ( ! residues_to_match[ci] ) continue;


		//CDR Residue Atoms:
		for ( core::Size cx = 1; cx <= pose.residue(ci).natoms(); ++cx ) {
			core::Size contacts = 0;

			//Antigen Residues:
			for ( core::Size ai = 1; ai <= pose.total_residue(); ++ai ) {
				if ( ! antigen_residues[ai] ) continue;

				//Antigen Residue Atoms:
				for ( core::Size ax = 1; ax <= pose.residue(ai).natoms(); ++ ax ) {
					numeric::xyzVector<core::Real> cx_xyz= pose.residue(ci).xyz(cx);
					numeric::xyzVector<core::Real> ax_xyz= pose.residue(ai).xyz(ax);
					core::Real cx_ax_dis = cx_xyz.distance(ax_xyz);
					if ( cx_ax_dis <= dist_cutoff ) {
						++contacts;
						//TR << "Contacts: " << contacts << std::endl;
					}
				}
			}

			if ( contacts > contact_cutoff ) ++ag_ab_atomic_contacts_[ci];
		}

	}

	/* Cannot get this to work correctly.  Giving answers for some cases, but not others, when they are clearly antigen-ab contacts.
	std::map<core::Size, core::id::AtomID> vertex_map;
	core::Size v = 1;

	for (core::Size i = 1; i <= pose.total_residue(); ++i){
	//Skip residues not needed at all.
	if (residues_to_match[i] || antigen_residues[i]){
	for (core::Size x = 1; x <= pose.residue(i).natoms(); ++x){
	vertex_map[v] = core::id::AtomID(x, i);
	++v;
	//TR << "V: " << v << " "<<x<<": "<<i << std::endl;
	}
	}
	}

	//Setup the PointGraph to be used in find neighbors algorithms
	conformation::PointGraphOP pg( new conformation::PointGraph );
	pg->set_num_vertices(vertex_map.size());

	for (core::Size i = 1; i <= vertex_map.size(); ++i){
	pg->get_vertex(i).data().xyz() = pose.residue(vertex_map[i].rsd()).xyz(vertex_map[i].atomno());
	}

	find_neighbors<PointGraphVertexData,PointGraphEdgeData>(pg, dist_cutoff);

	///Now we iterate over the atoms in our vector. For each atom, we get its edges.  If the vertex connected
	// corresponds to an antigen residue, we count.
	for (core::Size i = 1; i <= vertex_map.size(); ++i){
	if (residues_to_match[vertex_map[i].rsd()]){
	TR <<"matching res: " << vertex_map[i].rsd() << std::endl;
	core::Size contacts = 0;
	for ( PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( i ).upper_edge_list_begin(),
	edge_end_iter = pg->get_vertex( i ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {

	core::Size connected_res = vertex_map[edge_iter->upper_vertex()].rsd();
	//if (connected_res != vertex_map[i].rsd()){
	// TR << "contacted res: " << connected_res << std::endl;
	//}
	if (antigen_residues[connected_res]){
	++contacts;
	TR << "Contacts: " << contacts << std::endl;
	}
	}

	if (contacts > contact_cutoff) ++ag_ab_atomic_contacts_[vertex_map[i].rsd()];
	}
	}
	*/

}

core::Size
AntibodyFeatures::calculate_cdr_contacts_total(const core::pose::Pose & pose, CDRNameEnum const cdr){
	core::Size cdr_start = ab_info_->get_CDR_start(cdr, pose);
	core::Size cdr_end = ab_info_->get_CDR_end(cdr, pose);

	core::Size counts = 0;
	for ( core::Size i = cdr_start; i <= cdr_end; ++i ) {
		counts = counts + ag_ab_atomic_contacts_[i];
	}

	return counts;
}

core::Size
AntibodyFeatures::calculate_cdr_contacts_nres(const core::pose::Pose& pose, CDRNameEnum const cdr){


	core::Size cdr_start = ab_info_->get_CDR_start(cdr, pose);
	core::Size cdr_end = ab_info_->get_CDR_end(cdr, pose);

	core::Size counts = 0;
	for ( core::Size i = cdr_start; i <= cdr_end; ++i ) {
		if ( ag_ab_atomic_contacts_[i] > 0 ) ++counts;
	}

	TR << "CDR contacts: "<< ab_info_->get_CDR_name(cdr) << " "<< counts << std::endl;
	return counts;
}

} //antibody
} //protocols
