// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/StrandBundleFeatures.cc
/// @brief extract beta strand, strand pairs, sandwiches in pdb file, see wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StrandBundleFeatures for detail
/// @author Doo Nam Kim (based on Tim Jacobs' helix_assembly)
/// @overview
///  @ task 1: Identify all beta-strands
///   @ task 1-1: Write beta-strands into database
///  @ task 2: Identify all beta-sheets with these strands
///   @ task 2-1: Identify beta-sheets if their strands' two consecutive N-O pairs H-bond to each other
///   @ task 2-2: Write beta-sheets into database
///  @ task 3: Identify all beta-sandwiches with these sheets
///   @ task 3-1: Write beta-sandwiches into database

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//External
#include <boost/uuid/uuid.hpp>

//Devel
#include <protocols/features/strand_assembly/StrandBundleFeatures.hh>
#include <protocols/features/strand_assembly/StrandFragment.hh>

//Utility and basic
#include <basic/options/option.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <numeric/xyz.functions.hh> // for torsion calculations
#include <utility/vector1.hh> // for utility::vector1<Column> primary_key_columns;

//C library
#include <string>
#include <math.h> // for round

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/strand_assembly.OptionKeys.gen.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.hh> // for SASA
#include <core/scoring/ScoreFunction.hh> // ScoreFunction.hh seems required for compilation of InterfaceAnalyzerMover.hh
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/strand_assembly/StrandBundleFeaturesCreator.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

static THREAD_LOCAL basic::Tracer TR( "protocols.features.strand_assembly.StrandBundleFeatures" );

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StrandBundleFeatures::StrandBundleFeatures() :
	min_num_strands_to_deal_(5), // it should be at least 4
	max_num_strands_to_deal_(13), // (in 1LD9 chain A) 16 is too many number of strands, takes too long time
	extract_native_only_(false), // if true, extract native full strands only
	min_res_in_strand_(3), // min_res_in_strand_ used to be 4, but 3 would be OK
	max_res_in_strand_(22),
	min_O_N_dis_(2.4), // 1KIT shows (renumbered residues 178-181 and residues 76-81) show that 2.5 A exist!
	max_O_N_dis_(3.1),
	min_sheet_dis_(7.0), // 7 Angstrom may seem OK though
	max_sheet_dis_(15.0), // 15 Angstrom may seem OK though
	min_sheet_torsion_(-60.0), // although swissmodel.expasy.org/course/text/chapter4.htm mentions that -20 < torsion < -50
	max_sheet_torsion_(-5.0), // [1A5D] have -7.6 and -13.8 sheet torsion angles (measured by strands)
	min_sheet_angle_(30.0),
	max_sheet_angle_(150.0), // as Doonam observed, even 155 degree comes from same sheet (1ten)!
	min_shortest_dis_sidechain_inter_sheet_(2.0)
{
	init_from_options();
}


void StrandBundleFeatures::init_from_options(){
	using namespace basic::options;
	if ( option[OptionKeys::strand_assembly::min_num_strands_to_deal].user() ) {
		min_num_strands_to_deal_ = option[OptionKeys::strand_assembly::min_num_strands_to_deal];
	}
	if ( option[OptionKeys::strand_assembly::max_num_strands_to_deal].user() ) {
		max_num_strands_to_deal_ = option[OptionKeys::strand_assembly::max_num_strands_to_deal];
	}
	if ( option[OptionKeys::strand_assembly::extract_native_only].user() ) {
		extract_native_only_ = option[OptionKeys::strand_assembly::extract_native_only];
	}
	if ( option[OptionKeys::strand_assembly::min_res_in_strand].user() ) {
		min_res_in_strand_ = option[OptionKeys::strand_assembly::min_res_in_strand];
	}
	if ( option[OptionKeys::strand_assembly::max_res_in_strand].user() ) {
		max_res_in_strand_ = option[OptionKeys::strand_assembly::max_res_in_strand];
	}
	if ( option[OptionKeys::strand_assembly::min_O_N_dis].user() ) {
		min_O_N_dis_ = option[OptionKeys::strand_assembly::min_O_N_dis];
	}
	if ( option[OptionKeys::strand_assembly::max_O_N_dis].user() ) {
		max_O_N_dis_ = option[OptionKeys::strand_assembly::max_O_N_dis];
	}
	if ( option[OptionKeys::strand_assembly::min_sheet_dis].user() ) {
		min_sheet_dis_ = option[OptionKeys::strand_assembly::min_sheet_dis];
	}
	if ( option[OptionKeys::strand_assembly::max_sheet_dis].user() ) {
		max_sheet_dis_ = option[OptionKeys::strand_assembly::max_sheet_dis];
	}
	if ( option[OptionKeys::strand_assembly::min_sheet_torsion].user() ) {
		min_sheet_torsion_ = option[OptionKeys::strand_assembly::min_sheet_torsion];
	}
	if ( option[OptionKeys::strand_assembly::max_sheet_torsion].user() ) {
		max_sheet_torsion_ = option[OptionKeys::strand_assembly::max_sheet_torsion];
	}
	if ( option[OptionKeys::strand_assembly::min_sheet_angle].user() ) {
		min_sheet_angle_ = option[OptionKeys::strand_assembly::min_sheet_angle];
	}
	if ( option[OptionKeys::strand_assembly::max_sheet_angle].user() ) {
		max_sheet_angle_ = option[OptionKeys::strand_assembly::max_sheet_angle];
	}
	if ( option[OptionKeys::strand_assembly::min_shortest_dis_sidechain_inter_sheet].user() ) {
		min_shortest_dis_sidechain_inter_sheet_ = option[OptionKeys::strand_assembly::min_shortest_dis_sidechain_inter_sheet];
	}
} // init_from_options()

utility::vector1<std::string>
StrandBundleFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
} //features_reporter_dependencies()

void
StrandBundleFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	/****** <begin> write beta_selected_segments ******/

	// PrimaryKey
	// id of beta_selected_segments
	Column beta_selected_segments_id ("beta_selected_segments_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

	// unique key of original PDB file
	Column struct_id             ("struct_id", DbDataTypeOP( new DbBigInt() ),    false /*not null*/, false /*don't autoincrement*/);

	// ForeignKey
	Column residue_begin("residue_begin", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);
	Column residue_end  ("residue_end", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);

	utility::vector1<Column> primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(beta_selected_segments_id);

	// Schema
	// PrimaryKey
	Schema beta_selected_segments("beta_selected_segments",  PrimaryKey(primary_key_columns));

	// ForeignKey
	beta_selected_segments.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	// fkey_reference_cols -> one definition of fkey_reference_cols is enough since it will be used repeatedly
	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	// residue_begin_fkey_cols
	utility::vector1<Column> residue_begin_fkey_cols;
	residue_begin_fkey_cols.push_back(struct_id);
	residue_begin_fkey_cols.push_back(residue_begin);

	beta_selected_segments.add_foreign_key(ForeignKey(residue_begin_fkey_cols, "residues", fkey_reference_cols, true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/OneBodyFeaturesReporters#ResidueFeatures

	// residue_end_fkey_cols
	utility::vector1<Column> residue_end_fkey_cols;
	residue_end_fkey_cols.push_back(struct_id);
	residue_end_fkey_cols.push_back(residue_end);

	beta_selected_segments.add_foreign_key(ForeignKey(residue_end_fkey_cols, "residues", fkey_reference_cols, true /*defer*/));

	beta_selected_segments.write(db_session);

	/****** <end> write beta_selected_segments ******/


	/****** <begin> writing strand_pairs ******/

	// Columns
	// id of strand_pairs
	Column strand_pairs_id              ("strand_pairs_id",               DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

	// bool_parallel is PrimaryKey just because it doesn't point to any foreign values
	Column bool_parallel         ("bool_parallel",          DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);

	// unique key of original PDB file
	// Column struct_id             ("struct_id",              new DbBigInt(),    false /*not null*/, false /*don't autoincrement*/);

	// ForeignKey
	Column beta_select_id_i ("beta_select_id_i", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);
	Column beta_select_id_j ("beta_select_id_j", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);

	// Schema - strand_pairs
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sp;
	primary_key_columns_sp.push_back(struct_id);
	primary_key_columns_sp.push_back(strand_pairs_id);

	Schema strand_pairs("strand_pairs",  PrimaryKey(primary_key_columns_sp));

	// add column which is not PrimaryKey nor ForeignKey
	strand_pairs.add_column(bool_parallel);

	// ForeignKey

	strand_pairs.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols_beta_select;
	fkey_reference_cols_beta_select.push_back("struct_id");
	fkey_reference_cols_beta_select.push_back("beta_selected_segments_id");

	utility::vector1<Column> i_fkey_cols;
	i_fkey_cols.push_back(struct_id);
	i_fkey_cols.push_back(beta_select_id_i);

	strand_pairs.add_foreign_key(ForeignKey(i_fkey_cols, "beta_selected_segments", fkey_reference_cols_beta_select, true /*defer*/));

	utility::vector1<Column> j_fkey_cols;
	j_fkey_cols.push_back(struct_id);
	j_fkey_cols.push_back(beta_select_id_j);

	strand_pairs.add_foreign_key(ForeignKey(j_fkey_cols, "beta_selected_segments", fkey_reference_cols_beta_select, true /*defer*/));

	strand_pairs.write(db_session);

	/****** <end> writing strand_pairs ******/


	/****** <begin> writing sandwich ******/

	// PrimaryKey
	Column sandwich_id ("sandwich_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no-autoincrement*/);

	Column sp_id_1 ("sp_id_1", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);
	Column sp_id_2 ("sp_id_2", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);

	Column shortest_sc_dis ("shortest_sc_dis", DbDataTypeOP( new DbDouble() ), false /*not null*/, false /*don't autoincrement*/);

	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sw;
	primary_key_columns_sw.push_back(struct_id);
	primary_key_columns_sw.push_back(sandwich_id);

	Schema sandwich("sandwich",  PrimaryKey(primary_key_columns_sw));

	// ForeignKey
	sandwich.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	utility::vector1<std::string> fkey_reference_cols_strand_pairs;
	fkey_reference_cols_strand_pairs.push_back("struct_id");
	fkey_reference_cols_strand_pairs.push_back("strand_pairs_id");


	utility::vector1<Column> sp_1_fkey_cols;
	sp_1_fkey_cols.push_back(struct_id);
	sp_1_fkey_cols.push_back(sp_id_1);

	sandwich.add_foreign_key(ForeignKey(sp_1_fkey_cols, "strand_pairs", fkey_reference_cols_strand_pairs, true /*defer*/));

	utility::vector1<Column> sp_2_fkey_cols;
	sp_2_fkey_cols.push_back(struct_id);
	sp_2_fkey_cols.push_back(sp_id_2);

	sandwich.add_foreign_key(ForeignKey(sp_2_fkey_cols, "strand_pairs", fkey_reference_cols_strand_pairs, true /*defer*/));


	sandwich.add_column(shortest_sc_dis);

	sandwich.write(db_session);

	/****** <end> writing sandwich ******/



	/****** <begin> writing node *****/

	// Columns
	// id of node
	Column node_id ("node_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /* not autoincrement*/);

	Column bss_id_1 ("bss_id_1", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);
	Column bss_id_2 ("bss_id_2", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);

	// Schema - node
	// PrimaryKey
	Columns primary_key_columns_node;
	primary_key_columns_node.push_back(node_id);
	primary_key_columns_node.push_back(struct_id);

	Schema node("node",  PrimaryKey(primary_key_columns_node));


	// ForeignKey
	node.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	// bss_id_1_fkey_cols
	utility::vector1<Column> bss_id_1_fkey_cols;
	bss_id_1_fkey_cols.push_back(struct_id);
	bss_id_1_fkey_cols.push_back(bss_id_1); // not residues, bss_id

	// fkey_reference_cols_node -> one time definition
	utility::vector1<std::string> fkey_reference_cols_node;
	fkey_reference_cols_node.push_back("struct_id");
	fkey_reference_cols_node.push_back("beta_selected_segments_id");

	node.add_foreign_key(ForeignKey(bss_id_1_fkey_cols, "beta_selected_segments", fkey_reference_cols_node, true /*defer*/));

	// bss_id_2_fkey_cols
	utility::vector1<Column> bss_id_2_fkey_cols;
	bss_id_2_fkey_cols.push_back(struct_id);
	bss_id_2_fkey_cols.push_back(bss_id_2);

	node.add_foreign_key(ForeignKey(bss_id_2_fkey_cols, "beta_selected_segments", fkey_reference_cols_node, true /*defer*/));


	// bss_id_3_fkey_cols
	utility::vector1<Column> bss_id_3_fkey_cols;
	bss_id_3_fkey_cols.push_back(struct_id);
	bss_id_3_fkey_cols.push_back(sandwich_id);

	// fkey_reference_cols_node_sw -> one time definition
	utility::vector1<std::string> fkey_reference_cols_node_sw;
	fkey_reference_cols_node_sw.push_back("struct_id");
	fkey_reference_cols_node_sw.push_back("sandwich_id");

	node.add_foreign_key(ForeignKey(bss_id_3_fkey_cols, "sandwich", fkey_reference_cols_node_sw, true /*defer*/));

	node.write(db_session);

	/****** <end> writing node ****/

}

//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<StrandFragment> StrandBundleFeatures::get_full_strands(StructureID struct_id, sessionOP db_session)
{
	std::string select_string =
		"SELECT\n"
		"\tsegment_id,\n"
		"\tresidue_begin,\n"
		"\tresidue_end\n"
		"FROM\n"
		"\tsecondary_structure_segments\n"
		"WHERE\n"
		"\tstruct_id = ? AND (dssp = 'E');";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> all_strands;
	while ( res.next() )
			{
		Size strand_id,     residue_begin,   residue_end;
		res >> strand_id >> residue_begin >> residue_end;
		all_strands.push_back(StrandFragment(residue_begin, residue_end));
	}

	return all_strands;
}

//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<StrandFragment> StrandBundleFeatures::get_selected_strands(StructureID struct_id, sessionOP db_session)
{
	std::string select_string =
		"SELECT\n"
		"\tbeta_selected_segments_id,\n"
		"\tresidue_begin,\n"
		"\tresidue_end\n"
		"FROM\n"
		"\tbeta_selected_segments \n"
		"WHERE\n"
		"\tstruct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> all_strands;
	while ( res.next() )
			{
		Size bss_id,     residue_begin,   residue_end;
		res >> bss_id >> residue_begin >> residue_end;
		all_strands.push_back(StrandFragment(bss_id, residue_begin, residue_end));
	}

	return all_strands;
}

//Select all strand_pairs reported by the StrandBundleFeatures::get_strand_pairs and save them in a vector
utility::vector1<StrandFragment> StrandBundleFeatures::get_strand_pairs(StructureID struct_id, sessionOP db_session)
{
	std::string select_string =
		"SELECT\n"
		"\tstrand_pairs_id,\n"
		"\tbeta_select_id_i,\n"
		"\tbeta_select_id_j\n"
		"FROM\n"
		"\tstrand_pairs \n"
		"WHERE\n"
		"\tstruct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> retrieved_sp;
	while ( res.next() )
			{
		Size   strand_pairs_id,   beta_select_id_i,   beta_select_id_j;
		res >> strand_pairs_id >> beta_select_id_i >> beta_select_id_j;
		retrieved_sp.push_back(StrandFragment(strand_pairs_id, beta_select_id_i, beta_select_id_j));
	}
	return retrieved_sp;
}


utility::vector1<StrandFragment> StrandBundleFeatures::get_strand_from_bss_id(StructureID struct_id, sessionOP db_session, core::Size beta_selected_segments_id){
	std::string select_string =
		"SELECT\n"
		"\tbeta_selected_segments_id,\n"
		"\tresidue_begin,\n"
		"\tresidue_end\n"
		"FROM\n"
		"\tbeta_selected_segments \n"
		"WHERE\n"
		"\tstruct_id = ? AND beta_selected_segments_id = ? ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,beta_selected_segments_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> strand_from_bss_id;
	while ( res.next() ) {
		Size bss_id,     residue_begin,   residue_end;
		res >> bss_id >> residue_begin >> residue_end;
		strand_from_bss_id.push_back(StrandFragment(bss_id, residue_begin, residue_end));
	}
	return strand_from_bss_id;
}


bool
StrandBundleFeatures::find_antiparallel(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j)
{
	// seeing distances between 'O' of strand "i" and 'N' of strand "j"
	for ( Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++ ) {
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for ( Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++ ) {

			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_N_O = pose.residue(i_resnum).atom("N").xyz().distance(pose.residue(j_resnum).atom("O").xyz());
			if ( dis_N_O > min_O_N_dis_ && dis_N_O < max_O_N_dis_ ) {
				Real dis_O_N = pose.residue(i_resnum).atom("O").xyz().distance(pose.residue(j_resnum).atom("N").xyz());
				if ( dis_O_N > min_O_N_dis_ && dis_O_N < max_O_N_dis_ ) {
					if ( i_resnum+2 > strand_i.get_end() ) {
						continue; // I want to extract strand_pairs with only given ranges
					}
					if ( j_resnum-2 < strand_j.get_start() ) {
						continue; // I want to extract strand_pairs with only given ranges
					}
					Real dis_N_O = pose.residue(i_resnum+2).atom("N").xyz().distance(pose.residue(j_resnum-2).atom("O").xyz());
					if ( dis_N_O > min_O_N_dis_ && dis_N_O < max_O_N_dis_ ) {
						Real dis_O_N = pose.residue(i_resnum+2).atom("O").xyz().distance(pose.residue(j_resnum-2).atom("N").xyz());
						if ( dis_O_N > min_O_N_dis_ && dis_O_N < max_O_N_dis_ ) {
							return true;
						}
					}
				}
			}
		}
	}
	return false;
} //StrandBundleFeatures::find_antiparallel


bool
StrandBundleFeatures::find_parallel(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j)
{
	// seeing distances between 'O' of strand "i" and 'N' of strand "j"
	for ( Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++ ) {
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for ( Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++ ) {
			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_O_N = pose.residue(i_resnum).atom("O").xyz().distance(pose.residue(j_resnum).atom("N").xyz());
			if ( dis_O_N > min_O_N_dis_ && dis_O_N < max_O_N_dis_ ) {

				if ( i_resnum+2 > strand_i.get_end() ) {
					continue; // I want to extract strand_pairs with only given ranges
				}
				if ( j_resnum-2 < strand_j.get_start() ) {
					continue; // I want to extract strand_pairs with only given ranges
				}
				Real dis_N_O = pose.residue(i_resnum+2).atom("N").xyz().distance(pose.residue(j_resnum).atom("O").xyz());
				if ( dis_N_O > min_O_N_dis_ && dis_N_O < max_O_N_dis_ ) {
					Real dis_N_O = pose.residue(i_resnum).atom("N").xyz().distance(pose.residue(j_resnum-2).atom("O").xyz());
					if ( dis_N_O > min_O_N_dis_ && dis_N_O < max_O_N_dis_ ) {
						Real dis_O_N = pose.residue(i_resnum+2).atom("O").xyz().distance(pose.residue(j_resnum+2).atom("N").xyz());
						if ( dis_O_N > min_O_N_dis_ && dis_O_N < max_O_N_dis_ ) {
							return true;
						}
					}
				}
			}
		}
	}
	return false;
} //StrandBundleFeatures::find_parallel


Real
StrandBundleFeatures::get_avg_dis_CA_CA(
	Pose const & pose,
	Size i_resnum,
	Size i_resnum_1,
	Size i_resnum_2,
	Size i_resnum_3,
	Size j_resnum,
	Size j_resnum_1,
	Size j_resnum_2,
	Size j_resnum_3)
{
	Real dis_CA_CA_0 = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
	if ( dis_CA_CA_0 > 40 ) {
		return -999; // these sheets will not be sandwich ever, since these two sheets are too distant!
	}
	if ( dis_CA_CA_0 < min_sheet_dis_ || dis_CA_CA_0 > max_sheet_dis_ ) {
		return -99;
	}

	Real dis_CA_CA_1 = pose.residue(i_resnum_1).atom("CA").xyz().distance(pose.residue(j_resnum_1).atom("CA").xyz());

	if ( dis_CA_CA_1 < min_sheet_dis_ || dis_CA_CA_1 > max_sheet_dis_ ) {
		return -99;
	}

	Real dis_CA_CA_2 = pose.residue(i_resnum_2).atom("CA").xyz().distance(pose.residue(j_resnum_2).atom("CA").xyz());

	if ( dis_CA_CA_2 < min_sheet_dis_ || dis_CA_CA_2 > max_sheet_dis_ ) {
		return -99;
	}

	Real dis_CA_CA_3 = pose.residue(i_resnum_3).atom("CA").xyz().distance(pose.residue(j_resnum_3).atom("CA").xyz());

	if ( dis_CA_CA_3 < min_sheet_dis_ || dis_CA_CA_3 > max_sheet_dis_ ) {
		return -99;
	}

	Real avg_dis_CA_CA = (dis_CA_CA_0 + dis_CA_CA_1 + dis_CA_CA_2 + dis_CA_CA_3)/4;
	return avg_dis_CA_CA;
} // StrandBundleFeatures::get_avg_dis_CA_CA

Real
StrandBundleFeatures::check_sw_by_dis(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j,
	bool antiparalell // if false, find parallel way
)
{
	// TR.Info << "let me see distances between strands " << endl;
	Size i_resnum_1;
	Size j_resnum_1;
	Size i_resnum_2;
	Size j_resnum_2;
	Size i_resnum_3;
	Size j_resnum_3;
	for ( Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++ ) {
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for ( Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++ ) {
			Size j_resnum = strand_j.get_start()+strand_j_res;
			if ( antiparalell ) {
				i_resnum_1 = i_resnum+1;
				j_resnum_1 = j_resnum-1;

				i_resnum_2 = i_resnum+2;
				j_resnum_2 = j_resnum-2;

				i_resnum_3 = i_resnum+3;
				j_resnum_3 = j_resnum-3;

				if ( j_resnum_3 <= 0
						|| i_resnum_3 > pose.size()
						|| j_resnum_3 > pose.size() ) { // sometimes, j_resnum_3 becomes 18446744073709551615 where it should be -1
					continue;
				}
			} else { // paralell
				i_resnum_1 = i_resnum+1;
				j_resnum_1 = j_resnum+1;

				i_resnum_2 = i_resnum+2;
				j_resnum_2 = j_resnum+2;

				i_resnum_3 = i_resnum+3;
				j_resnum_3 = j_resnum+3;

				if ( i_resnum_3 > pose.size() || j_resnum_3 > pose.size() ) {
					continue;
				}
			}
			Real avg_dis_CA_CA = get_avg_dis_CA_CA(pose, i_resnum, i_resnum_1, i_resnum_2, i_resnum_3, j_resnum, j_resnum_1, j_resnum_2, j_resnum_3);
			if ( avg_dis_CA_CA == -999 ) {
				break; // these sheets will not be sandwich ever, since these two sheets are too distant!
			}
			if ( avg_dis_CA_CA == -99 ) {
				continue;
			}
			return avg_dis_CA_CA;
		} //for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
	} //for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	return -99; // these sheets are not sandwich with these strands
} //StrandBundleFeatures::check_sw_by_dis

Real
StrandBundleFeatures::sheet_torsion(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j)
{
	// check anti-parallel sheet

	//I need to define 4 terminal residues correctly for torsion calculation
	Real dis_i_end_and_j_start = pose.residue(strand_i.get_end()).atom("CA").xyz().distance(pose.residue(strand_j.get_start()).atom("CA").xyz());
	Real dis_i_end_and_j_end   = pose.residue(strand_i.get_end()).atom("CA").xyz().distance(pose.residue(strand_j.get_end()).atom("CA").xyz());

	if ( dis_i_end_and_j_start < dis_i_end_and_j_end ) {
		Vector const& first_xyz    ( pose.residue(strand_i.get_start()).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(strand_i.get_end()).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(strand_j.get_start()).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(strand_j.get_end()).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		Real torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
		return torsion_i_j;
	} else {
		Vector const& first_xyz    ( pose.residue(strand_i.get_start()).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(strand_i.get_end()).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(strand_j.get_end()).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(strand_j.get_start()).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		Real torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
		return torsion_i_j;
	}
} //StrandBundleFeatures::sheet_torsion


// check whether these sheets are too close
bool
StrandBundleFeatures::check_strand_too_closeness(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j)
{
	// check anti-parallel sheet distance

	// first, check the shortest distance between the two strand_pairs
	// seeing distances between 'CA' of strand "i" and 'CA' of strand "j"
	for ( Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++ ) {
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for ( Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++ ) {
			Size j_resnum = strand_j.get_start()+strand_j_res;

			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());

			if ( dis_CA_CA < 6.0 ) {
				return true; //TR.Info << "these two pair of strands are within 6 Angstrom according to check_strand_too_closeness" << endl;
			}
		}
	}
	return false;   //TR.Info << "OK, these two strand_pairs are farther than 6 Angstrom according to check_strand_too_closeness" << endl;
} //StrandBundleFeatures::check_strand_too_closeness


bool
StrandBundleFeatures::judge_sw_dis_too_close(
	Pose const & pose,
	StrandFragment temp_strand_ii_i,
	StrandFragment temp_strand_ii_j,
	StrandFragment temp_strand_jj_i,
	StrandFragment temp_strand_jj_j)
{
	Real return_of_sw_dis_ii_i_jj_i = check_strand_too_closeness (pose, temp_strand_ii_i, temp_strand_jj_i);
	// check of distance between strand ii_i and strand jj_i
	Real return_of_sw_dis_ii_j_jj_j = check_strand_too_closeness (pose, temp_strand_ii_j, temp_strand_jj_j);
	Real return_of_sw_dis_ii_i_jj_j = check_strand_too_closeness (pose, temp_strand_ii_i, temp_strand_jj_j);
	// check of distance between strand ii_i and strand jj_j
	Real return_of_sw_dis_ii_j_jj_i = check_strand_too_closeness (pose, temp_strand_ii_j, temp_strand_jj_i);

	if ( return_of_sw_dis_ii_i_jj_i == true || return_of_sw_dis_ii_j_jj_j == true ||  return_of_sw_dis_ii_i_jj_j == true || return_of_sw_dis_ii_j_jj_i == true ) {
		return true;
	}
	return false;
} //judge_sw_dis_too_close


bool
StrandBundleFeatures::final_check_sw_by_dis (
	Pose const & pose,
	StrandFragment temp_strand_ii_i,
	StrandFragment temp_strand_ii_j,
	StrandFragment temp_strand_jj_i,
	StrandFragment temp_strand_jj_j,
	bool antiparalell // if false, find parallel way
)
{
	Real return_of_sw_dis_ii_i_jj_i = check_sw_by_dis (pose, temp_strand_ii_i, temp_strand_jj_i, antiparalell);   // check of distance between strand ii_i and strand jj_i
	Real return_of_sw_dis_ii_j_jj_j = check_sw_by_dis (pose, temp_strand_ii_j, temp_strand_jj_j, antiparalell);
	Real return_of_sw_dis_ii_i_jj_j = check_sw_by_dis (pose, temp_strand_ii_i, temp_strand_jj_j, antiparalell);
	Real return_of_sw_dis_ii_j_jj_i = check_sw_by_dis (pose, temp_strand_ii_j, temp_strand_jj_i, antiparalell);
	if ( return_of_sw_dis_ii_i_jj_i == -999 || return_of_sw_dis_ii_j_jj_j == -999 ||  return_of_sw_dis_ii_i_jj_j == -999 || return_of_sw_dis_ii_j_jj_i == -999 ) {
		//TR.Info << "these sheets will not be sandwich ever since these are too distant!" << endl;
		return false;
	}

	if ( (return_of_sw_dis_ii_i_jj_i != -99 && return_of_sw_dis_ii_j_jj_j != -99) ||   (return_of_sw_dis_ii_i_jj_j != -99 && return_of_sw_dis_ii_j_jj_i != -99) ) {
		return true;
	}
	return false;
} //final_check_sw_by_dis

bool
StrandBundleFeatures::judge_sw_torsion(
	Pose const & pose,
	StrandFragment temp_strand_ii_i,
	StrandFragment temp_strand_ii_j,
	StrandFragment temp_strand_jj_i,
	StrandFragment temp_strand_jj_j)
{
	Real sheet_torsion_1 = sheet_torsion(pose, temp_strand_ii_i, temp_strand_jj_i);
	Real sheet_torsion_2 = sheet_torsion(pose, temp_strand_ii_i, temp_strand_jj_j);
	Real sheet_torsion_3 = sheet_torsion(pose, temp_strand_ii_j, temp_strand_jj_i);
	Real sheet_torsion_4 = sheet_torsion(pose, temp_strand_ii_j, temp_strand_jj_j);
	Real sheet_torsion_avg = (sheet_torsion_1 + sheet_torsion_2 + sheet_torsion_3 + sheet_torsion_4)/4;
	if ( sheet_torsion_avg > min_sheet_torsion_ && sheet_torsion_avg < max_sheet_torsion_ ) {
		return true;
	}
	return false;
}  // judge_sw_torsion


Real
StrandBundleFeatures::judge_sw_inter_dis (
	Pose const & pose,
	StrandFragment temp_strand_ii_i,
	StrandFragment temp_strand_ii_j,
	StrandFragment temp_strand_jj_i,
	StrandFragment temp_strand_jj_j)
{
	Real shortest_dis_sc_ii_i_jj_i = shortest_dis_sidechain (pose, temp_strand_ii_i, temp_strand_jj_i);
	Real shortest_dis_sc_ii_i_jj_j = shortest_dis_sidechain (pose, temp_strand_ii_i, temp_strand_jj_j);
	Real shortest_dis_sc_ii_j_jj_i = shortest_dis_sidechain (pose, temp_strand_ii_j, temp_strand_jj_i);
	Real shortest_dis_sc_ii_j_jj_j = shortest_dis_sidechain (pose, temp_strand_ii_j, temp_strand_jj_j);

	Real val_shortest_dis_sidechain = get_shortest_among_4(shortest_dis_sc_ii_i_jj_i, shortest_dis_sc_ii_i_jj_j, shortest_dis_sc_ii_j_jj_i, shortest_dis_sc_ii_j_jj_j);

	if ( val_shortest_dis_sidechain > min_shortest_dis_sidechain_inter_sheet_ ) {
		return val_shortest_dis_sidechain;
	}
	return -999;
}  //judge_sw_inter_dis

Size
StrandBundleFeatures::round(
	Real x)
{
	Size rounded = static_cast <Size> (floor(x+.5));
	return rounded;
} //round


Size
StrandBundleFeatures::get_nearest_res_from_strand(
	Pose const & pose,
	StrandFragment strand_to_be_searched,
	Size cen_resnum)
{
	Real shortest_dis = 999; //temp
	Size nearest_resnum = 999; //temp
	for ( Size strand_res=0; strand_res < strand_to_be_searched.get_size(); strand_res++ ) {
		Size j_resnum = strand_to_be_searched.get_start()+strand_res;
		Real dis_CA_CA = pose.residue(cen_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
		if ( shortest_dis > dis_CA_CA ) {
			shortest_dis = dis_CA_CA;
			nearest_resnum = j_resnum;
		}
	}
	return nearest_resnum;
} //get_nearest_res_from_strand



bool
StrandBundleFeatures::judge_facing(
	Pose const & pose,
	StrandFragment strand_ii_i,
	StrandFragment strand_ii_j,
	StrandFragment strand_jj_i,
	StrandFragment strand_jj_j)
{
	Real to_be_rounded_ii_i = (
		strand_ii_i.get_start() + strand_ii_i.get_end()
		)/(2.0);

	Size cen_resnum_ii_i = round(to_be_rounded_ii_i);

	Size cen_resnum_ii_j = get_nearest_res_from_strand(
		pose,
		strand_ii_j,
		cen_resnum_ii_i);

	Real to_be_rounded_jj_i = (
		strand_jj_i.get_start() + strand_jj_i.get_end()
		)/(2.0);

	Size cen_resnum_jj_i = round(to_be_rounded_jj_i);
	Size cen_resnum_jj_j = get_nearest_res_from_strand(
		pose,
		strand_jj_j,
		cen_resnum_jj_i);
	Real dis_ii_i_and_jj_i = pose.residue(cen_resnum_ii_i).atom("CA").xyz().distance(pose.residue(cen_resnum_jj_i).atom("CA").xyz());
	Real dis_ii_i_and_jj_j = pose.residue(cen_resnum_ii_i).atom("CA").xyz().distance(pose.residue(cen_resnum_jj_j).atom("CA").xyz());
	Real dis_ii_j_and_jj_i = pose.residue(cen_resnum_ii_j).atom("CA").xyz().distance(pose.residue(cen_resnum_jj_i).atom("CA").xyz());
	Real dis_ii_j_and_jj_j = pose.residue(cen_resnum_ii_j).atom("CA").xyz().distance(pose.residue(cen_resnum_jj_j).atom("CA").xyz());

	Real
		shortest_dis_inter_strand = get_shortest_among_4(dis_ii_i_and_jj_i, dis_ii_i_and_jj_j, dis_ii_j_and_jj_i, dis_ii_j_and_jj_j);

	//temp
	Real angle_1 = 999.9;
	Real angle_2 = 999.9;

	if ( shortest_dis_inter_strand == dis_ii_i_and_jj_i ) {
		Vector const& first_1_xyz    ( pose.residue(cen_resnum_ii_j).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(cen_resnum_ii_i).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(cen_resnum_jj_i).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(cen_resnum_jj_j).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(cen_resnum_jj_i).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(cen_resnum_ii_i).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);
	} else if ( shortest_dis_inter_strand == dis_ii_i_and_jj_j ) {
		Vector const& first_1_xyz    ( pose.residue(cen_resnum_ii_j).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(cen_resnum_ii_i).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(cen_resnum_jj_j).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(cen_resnum_jj_i).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(cen_resnum_jj_j).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(cen_resnum_ii_i).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);
	} else if ( shortest_dis_inter_strand == dis_ii_j_and_jj_i ) {
		Vector const& first_1_xyz    ( pose.residue(cen_resnum_ii_i).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(cen_resnum_ii_j).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(cen_resnum_jj_i).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(cen_resnum_jj_j).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(cen_resnum_jj_i).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(cen_resnum_ii_j).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);
	} else { // shortest_dis_inter_strand == dis_ii_j_and_jj_j
		Vector const& first_1_xyz    ( pose.residue(cen_resnum_ii_i).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(cen_resnum_ii_j).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(cen_resnum_jj_j).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(cen_resnum_jj_i).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(cen_resnum_jj_j).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(cen_resnum_ii_j).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);
	}

	if ( angle_1 > min_sheet_angle_ && angle_1 < max_sheet_angle_ && angle_2 > min_sheet_angle_ && angle_2 < max_sheet_angle_ ) {
		return true; // these two strand_pairs face to each other properly, so constitute a sandwich
	} else {
		return false; // these two strand_pairs are linear or do not face to each other properly!
	}
} // StrandBundleFeatures::judge_facing



Real
StrandBundleFeatures::shortest_dis_sidechain(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j) // calculate shortest distances for each pair
{
	Real temp_shortest_dis = 9999;
	for ( Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++ ) {
		Size i_resnum = strand_i.get_start()+strand_i_res;

		for ( Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++ ) {
			Size j_resnum = strand_j.get_start()+strand_j_res;
			for ( Size i_AtomNum=1; i_AtomNum < pose.residue(i_resnum).natoms(); i_AtomNum++ ) {
				for ( Size j_AtomNum=1; j_AtomNum < pose.residue(j_resnum).natoms(); j_AtomNum++ ) {
					Real dis_sc_sc = pose.residue(i_resnum).xyz(i_AtomNum).distance(pose.residue(j_resnum).xyz(j_AtomNum));
					if ( temp_shortest_dis > dis_sc_sc ) {
						temp_shortest_dis = dis_sc_sc;
					}
				}
			}
		}
	}
	return temp_shortest_dis;
} //StrandBundleFeatures::shortest_dis_sidechain (former one which is complicated)


Real
StrandBundleFeatures::get_shortest_among_4(
	Real val_shortest_dis_sidechain_1,
	Real val_shortest_dis_sidechain_2,
	Real val_shortest_dis_sidechain_3,
	Real val_shortest_dis_sidechain_4) // shortest of shortest distances
{
	Real temp_shortest_dis = val_shortest_dis_sidechain_1;

	if ( temp_shortest_dis > val_shortest_dis_sidechain_2 ) {
		{
			temp_shortest_dis = val_shortest_dis_sidechain_2;
			if ( temp_shortest_dis > val_shortest_dis_sidechain_3 ) {
				temp_shortest_dis = val_shortest_dis_sidechain_3;
				if ( temp_shortest_dis > val_shortest_dis_sidechain_4 ) {
					temp_shortest_dis = val_shortest_dis_sidechain_4;
				}
			} else {
				if ( temp_shortest_dis > val_shortest_dis_sidechain_4 ) {
					temp_shortest_dis = val_shortest_dis_sidechain_4;
				}
			}
		}
	} else {
		if ( temp_shortest_dis > val_shortest_dis_sidechain_3 ) {
			temp_shortest_dis = val_shortest_dis_sidechain_3;
			if ( temp_shortest_dis > val_shortest_dis_sidechain_4 ) {
				temp_shortest_dis = val_shortest_dis_sidechain_4;
			}
		} else {
			if ( temp_shortest_dis > val_shortest_dis_sidechain_4 ) {
				temp_shortest_dis = val_shortest_dis_sidechain_4;
			}
		}
	}
	return temp_shortest_dis;
} //StrandBundleFeatures::get_shortest_among_4 (simple one with just four parameters)


Real
StrandBundleFeatures::sheet_dis_by_terminals(
	Pose const & pose,
	StrandFragment strand_i,
	StrandFragment strand_j)
{
	// check anti-parallel sheet distance
	Real dis_strand_s_s_CA_CA = pose.residue(strand_i.get_start()).atom("CA").xyz().distance(pose.residue(strand_j.get_start()).atom("CA").xyz());
	Real dis_strand_s_e_CA_CA = pose.residue(strand_i.get_start()).atom("CA").xyz().distance(pose.residue(strand_j.get_end()).atom("CA").xyz());
	Real dis_strand_e_s_CA_CA = pose.residue(strand_i.get_end()).atom("CA").xyz().distance(pose.residue(strand_j.get_start()).atom("CA").xyz());
	Real dis_strand_e_e_CA_CA = pose.residue(strand_i.get_end()).atom("CA").xyz().distance(pose.residue(strand_j.get_end()).atom("CA").xyz());

	// consider only shortest distance
	if ( dis_strand_s_s_CA_CA < dis_strand_s_e_CA_CA ) {
		Real dis_strand_s_CA_CA = dis_strand_s_s_CA_CA;
		if ( dis_strand_e_s_CA_CA < dis_strand_e_e_CA_CA ) {
			Real dis_strand_e_CA_CA = dis_strand_e_s_CA_CA;
			Real val_sheet_dis_by_terminals = dis_strand_s_CA_CA + dis_strand_e_CA_CA;
			return val_sheet_dis_by_terminals;
		} else {
			Real dis_strand_e_CA_CA = dis_strand_e_e_CA_CA;
			Real val_sheet_dis_by_terminals = dis_strand_s_CA_CA + dis_strand_e_CA_CA;
			return val_sheet_dis_by_terminals;
		}
	} else {
		Real dis_strand_s_CA_CA = dis_strand_s_e_CA_CA;
		if ( dis_strand_e_s_CA_CA < dis_strand_e_e_CA_CA ) {
			Real dis_strand_e_CA_CA = dis_strand_e_s_CA_CA;
			Real val_sheet_dis_by_terminals = dis_strand_s_CA_CA + dis_strand_e_CA_CA;
			return val_sheet_dis_by_terminals;
		} else {
			Real dis_strand_e_CA_CA = dis_strand_e_e_CA_CA;
			Real val_sheet_dis_by_terminals = dis_strand_s_CA_CA + dis_strand_e_CA_CA;
			return val_sheet_dis_by_terminals;
		}
	}
} //StrandBundleFeatures::sheet_dis_by_terminals



/// @brief collect all the feature data for the pose
core::Size
StrandBundleFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	//TR.Info << "======================= <begin> report_features =========================" << endl;
	Size beta_selected_segments_id_counter=1; //initial value
	Size strand_pairs_id_counter=1; //initial value
	Size sandwich_id_counter=1; //initial value
	Size node_id_counter=1; //initial value

	utility::vector1<StrandFragment> all_strands = get_full_strands(struct_id, db_session);

	if ( all_strands.size() < min_num_strands_to_deal_ || all_strands.size() > max_num_strands_to_deal_ ) {
		TR.Info << "Exit early since all_strands.size(): " << all_strands.size() << endl;
		return 0;
	}

	if ( !extract_native_only_ ) {
		for ( Size ii=1; ii<=all_strands.size(); ++ii ) {
			for ( Size temp_res_num = all_strands[ii].get_start(); temp_res_num <= all_strands[ii].get_end()-3; ++temp_res_num ) {
				for ( Size jj = min_res_in_strand_-1; jj <= max_res_in_strand_; jj++ ) {
					if ( temp_res_num + jj <= all_strands[ii].get_end() ) {
						string beta_selected_insert =  "INSERT INTO beta_selected_segments (beta_selected_segments_id, struct_id, residue_begin, residue_end)  VALUES (?,?,?,?);";
						statement beta_selected_insert_stmt(basic::database::safely_prepare_statement(beta_selected_insert, db_session));
						beta_selected_insert_stmt.bind(1,beta_selected_segments_id_counter);
						beta_selected_insert_stmt.bind(2,struct_id);
						beta_selected_insert_stmt.bind(3,temp_res_num);
						beta_selected_insert_stmt.bind(4,temp_res_num+jj);
						basic::database::safely_write_to_database(beta_selected_insert_stmt);
						beta_selected_segments_id_counter++;
					} //if (temp_res_num + jj =< all_strands[i].get_end())
				} //for(Size jj = min_res_in_strand_; jj <= max_res_in_strand_; ++jj)
			}
		}
		all_strands = get_selected_strands(struct_id, db_session);
	} //if (!extract_native_only_)


	// strand pairing begins
	for ( Size i=1; i<all_strands.size(); ++i ) { // I don't need the last strand since this double for loops are exhaustive search for all pairs of strands
		if ( all_strands[i].get_size() >= min_res_in_strand_ && all_strands[i].get_size() <= max_res_in_strand_ ) {
			// the legnth of this beta strand is between min_res_in_strand_ and max_res_in_strand_
			for ( Size j=i+1; j<=all_strands.size(); ++j ) { // I need the last strand for this second for loop
				if ( all_strands[j].get_size() >= min_res_in_strand_ && all_strands[j].get_size() <= max_res_in_strand_ ) {
					// the length of this beta strand is between min_res_in_strand_ and max_res_in_strand_ too
					//TR.Info << "-- both strands are long enough to be considered for strands pairing --" << endl;

					StrandFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
					StrandFragment temp_strand_j(all_strands[j].get_start(), all_strands[j].get_end());

					//TR.Info << "-- anti-parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand begins --- " << endl;

					bool return_of_find_antiparallel(false); // temporary 'false' designation
					bool return_of_find_parallel(false); // temporary 'false' designation

					return_of_find_antiparallel = find_antiparallel (pose, temp_strand_i, temp_strand_j);

					if ( !return_of_find_antiparallel ) {
						return_of_find_parallel = find_parallel (pose, temp_strand_i, temp_strand_j);
					}

					Size bool_parallel = 0; // temp 0

					if ( return_of_find_parallel ) {
						bool_parallel = 1;
					}

					if ( return_of_find_antiparallel || return_of_find_parallel ) {
						string pair_insert =
							"INSERT INTO strand_pairs (strand_pairs_id, struct_id, bool_parallel, beta_select_id_i, beta_select_id_j)  VALUES (?,?, ?,?,?);";
						statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,db_session));
						pair_insert_stmt.bind(1, strand_pairs_id_counter);
						pair_insert_stmt.bind(2, struct_id);
						pair_insert_stmt.bind(3, bool_parallel);
						pair_insert_stmt.bind(4, i); // bss_id of the first strand
						pair_insert_stmt.bind(5, j); // bss_id of the second strand
						basic::database::safely_write_to_database(pair_insert_stmt);
						strand_pairs_id_counter++;
					}
				}
			}
		}
	} // end of 'for' loop for i strand
	//TR << "============== <Done> saving sheets (pairs of strands) ==========" << endl;
	// <end> strand pairing


	// <begin> sheet pairing

	utility::vector1<StrandFragment> retrieved_sp = get_strand_pairs(struct_id, db_session);
	for ( Size ii=1; ii<retrieved_sp.size(); ++ii ) { // I don't need the last pair of strands in this
		Size ii_i = retrieved_sp[ii].get_i(); // beta_selected_segments_id_i
		Size ii_j = retrieved_sp[ii].get_j(); // beta_selected_segments_id_j

		for ( Size jj=ii+1; jj<=retrieved_sp.size(); ++jj ) { // I need the last pair of strands in this second 'for' loop
			Size jj_i = retrieved_sp[jj].get_i(); // beta_selected_segments_id_i
			Size jj_j = retrieved_sp[jj].get_j(); // beta_selected_segments_id_j
			if ( ii_i != jj_i && ii_i != jj_j && ii_j != jj_i && ii_j != jj_j ) {
				//TR.Info << "<begins> sandwich check by distance in anti-parallel direction between " << ii << " th strand_pair and " << jj << " th strand_pair" << endl;
				utility::vector1<StrandFragment> ii_i_strand = get_strand_from_bss_id(struct_id, db_session, ii_i);
				utility::vector1<StrandFragment> ii_j_strand = get_strand_from_bss_id(struct_id, db_session, ii_j);
				utility::vector1<StrandFragment> jj_i_strand = get_strand_from_bss_id(struct_id, db_session, jj_i);
				utility::vector1<StrandFragment> jj_j_strand = get_strand_from_bss_id(struct_id, db_session, jj_j);

				StrandFragment temp_strand_ii_i(ii_i_strand[1].get_start(), ii_i_strand[1].get_end());
				StrandFragment temp_strand_ii_j(ii_j_strand[1].get_start(), ii_j_strand[1].get_end());
				StrandFragment temp_strand_jj_i(jj_i_strand[1].get_start(), jj_i_strand[1].get_end());
				StrandFragment temp_strand_jj_j(jj_j_strand[1].get_start(), jj_j_strand[1].get_end());

				bool too_close_sheet = judge_sw_dis_too_close (pose, temp_strand_ii_i, temp_strand_ii_j, temp_strand_jj_i, temp_strand_jj_j);
				if ( too_close_sheet ) {
					continue;
				}

				bool sw_dis_anti = final_check_sw_by_dis (pose, temp_strand_ii_i, temp_strand_ii_j, temp_strand_jj_i, temp_strand_jj_j, true);
				bool sw_dis_parallel = final_check_sw_by_dis (pose, temp_strand_ii_i, temp_strand_ii_j, temp_strand_jj_i, temp_strand_jj_j, false);
				if ( !sw_dis_anti && !sw_dis_parallel ) {
					continue;
				}

				bool sw_torsion = judge_sw_torsion (pose, temp_strand_ii_i, temp_strand_ii_j, temp_strand_jj_i, temp_strand_jj_j);
				Real val_shortest_dis_sidechain = judge_sw_inter_dis (pose, temp_strand_ii_i, temp_strand_ii_j, temp_strand_jj_i, temp_strand_jj_j);
				bool facing = judge_facing(pose, temp_strand_ii_i, temp_strand_ii_j, temp_strand_jj_i, temp_strand_jj_j);
				// if false, these two strand_pairs are linear to each other or do not face properly to each other!

				if ( !sw_torsion || val_shortest_dis_sidechain == -999 || !facing ) {
					continue;
				}

				// TR.Info << "===== <begin> saving sandwich ===== " << endl;
				string sandwich_insert =
					"INSERT INTO sandwich (sandwich_id, struct_id, sp_id_1, sp_id_2, shortest_sc_dis)  VALUES (?,?,?,?,?);";

				statement sandwich_insert_stmt(basic::database::safely_prepare_statement(sandwich_insert, db_session));
				sandwich_insert_stmt.bind(1, sandwich_id_counter);
				// I should increase sandwich_id_counter later     // sandwich_id_counter++;
				sandwich_insert_stmt.bind(2, struct_id);
				sandwich_insert_stmt.bind(3, ii); // sp_id of the first strand_pairs
				sandwich_insert_stmt.bind(4, jj); // sp_id of the second strand_pairs
				sandwich_insert_stmt.bind(5, val_shortest_dis_sidechain);
				basic::database::safely_write_to_database(sandwich_insert_stmt);

				Real sheet_dis_by_term_ii_i_jj_i = sheet_dis_by_terminals (pose, temp_strand_ii_i, temp_strand_jj_i);
				Real sheet_dis_by_term_ii_i_jj_j = sheet_dis_by_terminals (pose, temp_strand_ii_i, temp_strand_jj_j);
				Real sheet_dis_by_term_ii_j_jj_i = sheet_dis_by_terminals (pose, temp_strand_ii_j, temp_strand_jj_i);
				Real sheet_dis_by_term_ii_j_jj_j = sheet_dis_by_terminals (pose, temp_strand_ii_j, temp_strand_jj_j);

				// I need to store closer strands together as a group for Tim's Sewing
				Real cross_dis = sheet_dis_by_term_ii_i_jj_j + sheet_dis_by_term_ii_j_jj_i;
				Real non_cross_dis = sheet_dis_by_term_ii_i_jj_i + sheet_dis_by_term_ii_j_jj_j;

				if ( cross_dis > non_cross_dis ) {
					string node_insert_1_1 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";
					statement node_insert_1_1_stmt(basic::database::safely_prepare_statement(node_insert_1_1, db_session));
					node_insert_1_1_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_1_1_stmt.bind(2, sandwich_id_counter);
					node_insert_1_1_stmt.bind(3, struct_id);
					node_insert_1_1_stmt.bind(4, ii_i);
					node_insert_1_1_stmt.bind(5, jj_i);
					basic::database::safely_write_to_database(node_insert_1_1_stmt);

					string node_insert_1_2 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";
					statement node_insert_1_2_stmt(basic::database::safely_prepare_statement(node_insert_1_2, db_session));
					node_insert_1_2_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_1_2_stmt.bind(2, sandwich_id_counter);
					node_insert_1_2_stmt.bind(3, struct_id);
					node_insert_1_2_stmt.bind(4, jj_i);
					node_insert_1_2_stmt.bind(5, ii_i);
					basic::database::safely_write_to_database(node_insert_1_2_stmt);


					string node_insert_2_1 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";

					statement node_insert_2_1_stmt(basic::database::safely_prepare_statement(node_insert_2_1, db_session));

					node_insert_2_1_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_2_1_stmt.bind(2, sandwich_id_counter);
					node_insert_2_1_stmt.bind(3, struct_id);
					node_insert_2_1_stmt.bind(4, ii_j);
					node_insert_2_1_stmt.bind(5, jj_j);

					basic::database::safely_write_to_database(node_insert_2_1_stmt);


					string node_insert_2_2 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";

					statement node_insert_2_2_stmt(basic::database::safely_prepare_statement(node_insert_2_2, db_session));

					node_insert_2_2_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_2_2_stmt.bind(2, sandwich_id_counter);
					node_insert_2_2_stmt.bind(3, struct_id);
					node_insert_2_2_stmt.bind(4, jj_j);
					node_insert_2_2_stmt.bind(5, ii_j);

					basic::database::safely_write_to_database(node_insert_2_2_stmt);

					sandwich_id_counter++; // I should increase sandwich_id_counter now
				} else { //(cross_dis > non_cross_dis) // cross_dis =< non_cross_dis
					string node_insert_1_1 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";

					statement node_insert_1_1_stmt(basic::database::safely_prepare_statement(node_insert_1_1, db_session));

					node_insert_1_1_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_1_1_stmt.bind(2, sandwich_id_counter);
					node_insert_1_1_stmt.bind(3, struct_id);
					node_insert_1_1_stmt.bind(4, ii_i);
					node_insert_1_1_stmt.bind(5, jj_j);

					basic::database::safely_write_to_database(node_insert_1_1_stmt);


					string node_insert_1_2 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";

					statement node_insert_1_2_stmt(basic::database::safely_prepare_statement(node_insert_1_2, db_session));

					node_insert_1_2_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_1_2_stmt.bind(2, sandwich_id_counter);
					node_insert_1_2_stmt.bind(3, struct_id);
					node_insert_1_2_stmt.bind(4, jj_j);
					node_insert_1_2_stmt.bind(5, ii_i);

					basic::database::safely_write_to_database(node_insert_1_2_stmt);


					string node_insert_2_1 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";

					statement node_insert_2_1_stmt(basic::database::safely_prepare_statement(node_insert_2_1, db_session));

					node_insert_2_1_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_2_1_stmt.bind(2, sandwich_id_counter);
					node_insert_2_1_stmt.bind(3, struct_id);
					node_insert_2_1_stmt.bind(4, ii_j);
					node_insert_2_1_stmt.bind(5, jj_i);

					basic::database::safely_write_to_database(node_insert_2_1_stmt);


					string node_insert_2_2 =
						"INSERT INTO node (node_id, sandwich_id, struct_id, bss_id_1, bss_id_2)  VALUES (?,?,?,?,?);";
					statement node_insert_2_2_stmt(basic::database::safely_prepare_statement(node_insert_2_2, db_session));
					node_insert_2_2_stmt.bind(1, node_id_counter);
					node_id_counter++;
					node_insert_2_2_stmt.bind(2, sandwich_id_counter);
					node_insert_2_2_stmt.bind(3, struct_id);
					node_insert_2_2_stmt.bind(4, jj_i);
					node_insert_2_2_stmt.bind(5, ii_j);
					basic::database::safely_write_to_database(node_insert_2_2_stmt);

					sandwich_id_counter++; // I should increase sandwich_id_counter now
				} // cross_dis =< non_cross_dis

			} // if bss are not same
		} // end of 'jj' for loop
	} // end of 'ii' for loop
	TR << "=== <Done> saving sandwichs (pairs of sheets) ===" << endl;
	// <end> sheet pairing
	return 0;
} //StrandBundleFeatures::report_features

std::string StrandBundleFeatures::type_name() const {
	return class_name();
}

std::string StrandBundleFeatures::class_name() {
	return "StrandBundleFeatures";
}

void StrandBundleFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Extract beta strand, strand pairs, sandwiches in pdb file and writes it into database",
		attlist );
}

std::string StrandBundleFeaturesCreator::type_name() const {
	return StrandBundleFeatures::class_name();
}

protocols::features::FeaturesReporterOP
StrandBundleFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new StrandBundleFeatures );
}

void StrandBundleFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StrandBundleFeatures::provide_xml_schema( xsd );
}



} //namespace strand_assembly
} //namespace features
} //namespace protocols
