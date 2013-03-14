// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/SandwichFeatures.cc
/// @brief analyze beta-sandwich features
/// @author Doo Nam Kim (based on Tim Jacobs' helix_assembly)


//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh> // for dssp application

//External
#include <boost/uuid/uuid.hpp>

//Devel
#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFragment.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <numeric/xyz.functions.hh> // for torsion calculations
#include <utility/vector1.hh> // for utility::vector1<Column> primary_key_columns;

//for vector
#include <numeric/xyzVector.hh>
#include <core/id/NamedAtomID.hh>

//C library
#include <string>
#include <math.h> // for round and sqrt

//External Headers
#include <cppdb/frontend.h>

//Basic rosetta
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/strand_assembly.OptionKeys.gen.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

// basic c++
#include <iostream>
#include <stdio.h>     //for remove( ) and rename( )
#include <fstream>
#include <string>
#include <vector> // for get_sw_can_by_sh_id, get_two_central_residues

// exception handling
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.hh> // for SASA
#include <core/scoring/ScoreFunction.hh> // ScoreFunction.hh seems required for compilation of InterfaceAnalyzerMover.hh

//DSSP
#include <core/scoring/dssp/Dssp.hh>

static basic::Tracer TR("protocols.features.strand_assembly.SandwichFeatures");

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

SandwichFeatures::SandwichFeatures() :
min_num_strands_to_deal_(4), // At least 4 strands should be in pdb file
max_num_strands_to_deal_(60),
min_strand_size_(2), // minimum number of residues in a strand, for edge strand definition & analysis, 4=< is recommended (in 1A8M) min_strand_size = 2, (in 1PMY) min_strand_size = 3
min_CA_CA_dis_(3.5), // definition: minimum CA_CA_distance between strands in same sheet
					 // example: (in 1A8M) 'min_CA_CA_dis_= 3.5', (in 1KIT) 'min_CA_CA_dis_= 4.0'
max_CA_CA_dis_(6.2), // (in 1A8M) 'max_CA_CA_dis_= 6.2', (in 1KIT) 'max_CA_CA_dis_= 5.7'
min_C_O_N_angle_(120), // (in 1L0Q chain A), 138 is the smallest C_O_N_angle (C and O from one sheet, N from other sheet)
min_sheet_dis_(7.0),	// definition: minimum CA_CA_distance between strands in different sheets to constitute a sandwich
						// 7 Angstrom seems OK
max_sheet_dis_(15.0),	// definition: maximum CA_CA_distance between strands in different sheets to constitute a sandwich
						// 15 Angstrom seems OK
min_sheet_angle_(30.0), // definition: Minimum angle between sheets (CA and CA)
						//	usage: used in judge_facing "angle_1 > min_sheet_angle_"
max_sheet_angle_(150.0), // (in 1TEN) even 155 degree comes from same sheet!
min_sheet_torsion_cen_res_(-150.0), //	definition: "Minimum torsion between sheets (CA and CA) with respect to terminal central residues in each beta-sheet
									//	explanation: with respect to central residues, one torsion angles of 1TEN is 84.9
max_sheet_torsion_cen_res_(150.0),  //	definition: "maximum torsion between sheets (CA and CA) with respect to terminal central residues in each beta-sheet
								    //	usage: used in judge_facing "torsion_i_j < max_sheet_torsion_cen_res_"
min_num_strands_in_sheet_(3),	//  definition: a sheet with < 3 strands will be ignored
								//	usage: if (num_strands_i < min_num_strands_in_sheet_)
min_inter_sheet_dis_CA_CA_(4.0),	//	example:	(in 12E8) the distance between S52 and G64 is 4.2 A
									//	example:	(in 2WYF) the distance between val and gln is 5.8 A
max_inter_sheet_dis_CA_CA_(24.0),	//	usage:	shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_
									//	example:	(in 1TEN) shortest_avg_dis_inter_sheet between sheet 1 and 2 = 11.6 A (and these sheets should be a sandwich)
									//	example:	(in 1A64_chain_A) shortest_avg_dis_inter_sheet between sheet 1 and 2 = 25 A (and these sheets should not be a sandwich)
									//	example:	(in 1A1M) the average distance between sheet 1 and 4 > 20 A (but these sheets should be a sandwich)
									//	example:	(in 1ASQ) the average distance between sheet 1 and 4 > 20 A (and these sheets should not be a sandwich)
extract_sandwich_(true),
write_chain_B_resnum_(true), // if true, write chain_B_resnum file for InterfaceAnalyzer
no_helix_in_pdb_(false), // if true, ignore any pdb that has helix
no_helix_in_extracted_sw_(true), // if true, ignore any pdb that has helix in extracted sandwich
no_strand_in_loop_in_extracted_sw_(false), // if true, ignore any pdb that has a strand as a loop in extracted sandwich
exclude_sandwich_that_is_linked_w_same_direction_strand_(true), // if true, exclude_sandwich_that_is_linked with same_direction_strand
write_phi_psi_(false), // if true, write phi_psi_file
max_num_sw_per_pdb_(50), // maximum number of sandwiches to be extracted per a pdb file
check_N_to_C_direction_by_("PE"), // check N->C going direction by option
check_canonicalness_cutoff_(80) //	definition:	cutoff to determine canonicalness of L/R, P/A and directionality
{
	init_from_options();
}

void
SandwichFeatures::init_from_options(){
	using namespace basic::options;

	if(option[OptionKeys::strand_assembly::min_num_strands_to_deal].user()){
		min_num_strands_to_deal_ = option[OptionKeys::strand_assembly::min_num_strands_to_deal];
	}
	if(option[OptionKeys::strand_assembly::max_num_strands_to_deal].user()){
		max_num_strands_to_deal_ = option[OptionKeys::strand_assembly::max_num_strands_to_deal];
	}
	if(option[OptionKeys::strand_assembly::min_strand_size].user()){
		min_strand_size_ = option[OptionKeys::strand_assembly::min_strand_size];
	}
	if(option[OptionKeys::strand_assembly::min_CA_CA_dis].user()){
		min_CA_CA_dis_ = option[OptionKeys::strand_assembly::min_CA_CA_dis];	
	}
	if(option[OptionKeys::strand_assembly::max_CA_CA_dis].user()){
		max_CA_CA_dis_ = option[OptionKeys::strand_assembly::max_CA_CA_dis];
	}
	if(option[OptionKeys::strand_assembly::min_C_O_N_angle].user()){
		min_C_O_N_angle_ = option[OptionKeys::strand_assembly::min_C_O_N_angle];
	}
	if(option[OptionKeys::strand_assembly::min_sheet_dis].user()){
		min_sheet_dis_ = option[OptionKeys::strand_assembly::min_sheet_dis];
	}
	if(option[OptionKeys::strand_assembly::max_sheet_dis].user()){
		max_sheet_dis_ = option[OptionKeys::strand_assembly::max_sheet_dis];
	}
	if(option[OptionKeys::strand_assembly::min_sheet_angle].user()){
		min_sheet_angle_ = option[OptionKeys::strand_assembly::min_sheet_angle];
	}
	if(option[OptionKeys::strand_assembly::max_sheet_angle].user()){
		max_sheet_angle_ = option[OptionKeys::strand_assembly::max_sheet_angle];
	}
	if(option[OptionKeys::strand_assembly::min_sheet_torsion_cen_res].user()){
		min_sheet_torsion_cen_res_ = option[OptionKeys::strand_assembly::min_sheet_torsion_cen_res];
	}
	if(option[OptionKeys::strand_assembly::max_sheet_torsion_cen_res].user()){
		max_sheet_torsion_cen_res_ = option[OptionKeys::strand_assembly::max_sheet_torsion_cen_res];
	}
	if(option[OptionKeys::strand_assembly::min_num_strands_in_sheet].user()){
		min_num_strands_in_sheet_ = option[OptionKeys::strand_assembly::min_num_strands_in_sheet];
	}
	if(option[OptionKeys::strand_assembly::min_inter_sheet_dis_CA_CA].user()){
		min_inter_sheet_dis_CA_CA_ = option[OptionKeys::strand_assembly::min_inter_sheet_dis_CA_CA];
	}
	if(option[OptionKeys::strand_assembly::max_inter_sheet_dis_CA_CA].user()){
		max_inter_sheet_dis_CA_CA_ = option[OptionKeys::strand_assembly::max_inter_sheet_dis_CA_CA];
	}
	if(option[OptionKeys::strand_assembly::extract_sandwich].user()){
		extract_sandwich_ = option[OptionKeys::strand_assembly::extract_sandwich];
	}
	if(option[OptionKeys::strand_assembly::write_chain_B_resnum].user()){
		write_chain_B_resnum_ = option[OptionKeys::strand_assembly::write_chain_B_resnum];
	}
	if(option[OptionKeys::strand_assembly::write_phi_psi].user()){
		write_phi_psi_ = option[OptionKeys::strand_assembly::write_phi_psi];
	}
	if(option[OptionKeys::strand_assembly::no_helix_in_pdb].user()){
		no_helix_in_pdb_ = option[OptionKeys::strand_assembly::no_helix_in_pdb];
	}
	if(option[OptionKeys::strand_assembly::no_helix_in_extracted_sw].user()){
		no_helix_in_extracted_sw_ = option[OptionKeys::strand_assembly::no_helix_in_extracted_sw];
	}
	if(option[OptionKeys::strand_assembly::no_strand_in_loop_in_extracted_sw].user()){
		no_strand_in_loop_in_extracted_sw_ = option[OptionKeys::strand_assembly::no_strand_in_loop_in_extracted_sw];
	}
	if(option[OptionKeys::strand_assembly::exclude_sandwich_that_is_linked_w_same_direction_strand].user()){
		exclude_sandwich_that_is_linked_w_same_direction_strand_ = option[OptionKeys::strand_assembly::exclude_sandwich_that_is_linked_w_same_direction_strand];
	}
	if(option[OptionKeys::strand_assembly::max_num_sw_per_pdb].user()){
		max_num_sw_per_pdb_ = option[OptionKeys::strand_assembly::max_num_sw_per_pdb];
	}
	if(option[OptionKeys::strand_assembly::check_N_to_C_direction_by].user()){
		check_N_to_C_direction_by_ = option[OptionKeys::strand_assembly::check_N_to_C_direction_by];
	}
	if(option[OptionKeys::strand_assembly::check_canonicalness_cutoff].user()){
		check_canonicalness_cutoff_ = option[OptionKeys::strand_assembly::check_canonicalness_cutoff];
	}
} // init_from_options()

utility::vector1<std::string>
SandwichFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
} //features_reporter_dependencies()

void
SandwichFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

/****** <begin> writing sheet ******/
	// Columns
	// id of sheet
	//unique
	Column sheet_PK_id	("sheet_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	//Column sheet_id	("sheet_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);
	Column sheet_id	("sheet_id",	new DbInteger(), true /*could be null*/, false /*no autoincrement*/);
	// changed into 'null-possible' because of sw_by_components

	Column sheet_antiparallel	("sheet_antiparallel",	new DbText(), true /* could be null at first, eventually it will not be null though*/, false /*no autoincrement*/);
	// A: antiparallel
	// P_or_mix: parallel or mix of antiparallel and parallel

	// unique key of original PDB file
	Column struct_id             ("struct_id",              new DbUUID(),    false /*not null*/, false /*don't autoincrement*/);

	// ForeignKey
	Column segment_id ("segment_id",	new DbInteger(), false /*not null*/, false /*don't autoincrement*/);

	// Schema - sheet
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sh;
	primary_key_columns_sh.push_back(struct_id);
	primary_key_columns_sh.push_back(sheet_PK_id);

	Schema sheet("sheet",  PrimaryKey(primary_key_columns_sh));

	// add column which is not PrimaryKey nor ForeignKey
	sheet.add_column(sheet_id);
	sheet.add_column(sheet_antiparallel);

	// ForeignKey
	sheet.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
		// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols_beta;
	fkey_reference_cols_beta.push_back("struct_id");
	fkey_reference_cols_beta.push_back("segment_id");

	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(segment_id);

	sheet.add_foreign_key(ForeignKey(fkey_cols,	"secondary_structure_segments",	fkey_reference_cols_beta,	true /*defer*/));

	sheet.write(db_session);
/****** <end> writing sheet ******/



/****** <begin> writing sw_can_by_sh (sandwich candidate by sheets) ******/
	// Columns
	// id of sw_can_by_sh
	//unique
	Column sw_can_by_sh_PK_id	("sw_can_by_sh_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column tag	("tag",	new DbText(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column sw_can_by_sh_id	("sw_can_by_sh_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column strand_num	("strand_num",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// Schema - sw_can_by_sh_id
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sw_can_by_sh;
	primary_key_columns_sw_can_by_sh.push_back(struct_id);
	primary_key_columns_sw_can_by_sh.push_back(sw_can_by_sh_PK_id);

	Schema sw_can_by_sh("sw_can_by_sh",  PrimaryKey(primary_key_columns_sw_can_by_sh));

	// add column which is not PrimaryKey nor ForeignKey
	sw_can_by_sh.add_column(tag);
	sw_can_by_sh.add_column(sw_can_by_sh_id);

	// ForeignKey
	sw_can_by_sh.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
		// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols_sh;
	fkey_reference_cols_sh.push_back("struct_id");
	fkey_reference_cols_sh.push_back("sheet_PK_id");
	//fkey_reference_cols_sh.push_back("sheet_id"); // sqlite3 OK but psql "ERROR: cppdb::posgresql: statement execution failed : ERROR:  there is no unique constraint matching given keys for referenced table "sheet""

	utility::vector1<Column> fkey_cols_sh;
	fkey_cols_sh.push_back(struct_id);
	fkey_cols_sh.push_back(sheet_id);

	sw_can_by_sh.add_foreign_key(ForeignKey(fkey_cols_sh,	"sheet",	fkey_reference_cols_sh,	true /*defer*/));

	// add column which is not PrimaryKey nor ForeignKey
	sw_can_by_sh.add_column(strand_num);

	sw_can_by_sh.write(db_session);
/****** <end> writing sw_can_by_sh ******/


/****** <begin> writing sw_by_components (sandwich candidate by components such as strands, loops, helices) ******/

	// Columns
	// id of sw_by_components

	//unique
	Column sw_by_components_PK_id	("sw_by_components_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column sw_by_components_bs_id	("sw_by_components_bs_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column sw_by_components_bs_edge	("sw_by_components_bs_edge",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column intra_sheet_con_id	("intra_sheet_con_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column inter_sheet_con_id	("inter_sheet_con_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column LR	("LR",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column cano_LR	("cano_LR",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	// T, -> true, canonical chiral
	// F, -> false, non-canonical chiral
	// U, -> uncertain, this loop-size with this condition has no definite canonical chiral reference in the first place!

	Column PA_by_preceding_E	("PA_by_preceding_E",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column PA_by_following_E	("PA_by_following_E",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column cano_PA	("cano_PA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	// T, -> true, canonical PA
	// F, -> false, non-canonical PA
	// U, -> uncertain, this loop-size with this condition has no definite canonical PA reference in the first place!

	Column heading_direction	("heading_direction",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column parallel_EE	("parallel_EE",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column cano_parallel_EE	("cano_parallel_EE",	new DbText(), true /* could be null*/, false /*no autoincrement*/);	
	Column loop_size	("loop_size",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column residue_begin("residue_begin", new DbInteger(), false /*not null*/, false /*don't autoincrement*/);
	Column residue_end  ("residue_end", new DbInteger(), false /*not null*/, false /*don't autoincrement*/);

	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sw_by_components;
	primary_key_columns_sw_by_components.push_back(struct_id);
	primary_key_columns_sw_by_components.push_back(sw_by_components_PK_id);

	Schema sw_by_components("sw_by_components",  PrimaryKey(primary_key_columns_sw_by_components));

	// add column which is not PrimaryKey nor ForeignKey
	sw_by_components.add_column(tag);
	sw_by_components.add_column(sw_can_by_sh_id);
	sw_by_components.add_column(sheet_id);
	sw_by_components.add_column(sheet_antiparallel);
	sw_by_components.add_column(sw_by_components_bs_id);
	sw_by_components.add_column(sw_by_components_bs_edge);
	sw_by_components.add_column(intra_sheet_con_id);
	sw_by_components.add_column(inter_sheet_con_id);
	sw_by_components.add_column(LR);
	sw_by_components.add_column(cano_LR);
	sw_by_components.add_column(PA_by_preceding_E);
	sw_by_components.add_column(PA_by_following_E);
	sw_by_components.add_column(cano_PA);
	sw_by_components.add_column(heading_direction);
	sw_by_components.add_column(parallel_EE);
	sw_by_components.add_column(cano_parallel_EE);
	sw_by_components.add_column(loop_size);

	// ForeignKey
	sw_by_components.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
		// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	utility::vector1<Column> residue_begin_fkey_cols;
	residue_begin_fkey_cols.push_back(struct_id);
	residue_begin_fkey_cols.push_back(residue_begin);

	sw_by_components.add_foreign_key(ForeignKey(residue_begin_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

	utility::vector1<Column> residue_end_fkey_cols;
	residue_end_fkey_cols.push_back(struct_id);
	residue_end_fkey_cols.push_back(residue_end);

	sw_by_components.add_foreign_key(ForeignKey(residue_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

	sw_by_components.write(db_session);
/****** <end> writing sw_by_components ******/
}


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
SandwichFeatures::get_full_strands(
   boost::uuids::uuid struct_id,
   sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	segment_id,\n"
	"	residue_begin,\n"
	"	residue_end\n"
	"FROM\n"
	"	secondary_structure_segments\n"
	"WHERE\n"
	"	dssp = 'E' \n"
	"	AND struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next()){
		Size strand_id,     residue_begin,   residue_end;
		res >> strand_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(residue_begin, residue_end));
	}
	return all_strands;
}


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
SandwichFeatures::get_full_strands_from_sheet(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	string select_string =
	"SELECT\n"
	"	bs.segment_id,\n"
	"	bs.residue_begin,\n"
	"	bs.residue_end\n"
	"FROM\n"
	"	secondary_structure_segments as bs, \n"
	"	sheet as sh \n"
	"WHERE\n"
	"	bs.struct_id = sh.struct_id \n"
	"	AND bs.dssp = 'E'\n"
	"	AND bs.struct_id = ? \n"
	"	AND sh.segment_id = bs.segment_id \n"
	"	AND sh.sheet_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sheet_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next())
	{
		Size strand_id,     residue_begin,   residue_end;
		res >> strand_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(residue_begin, residue_end));
	}
	return all_strands;
}


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
bool
SandwichFeatures::check_whether_strand_i_is_in_sheet(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size segment_id)
{
	string select_string =
	"SELECT\n"
	"	segment_id\n"
	"FROM\n"
	"	sheet\n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND segment_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,segment_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	bool strand_i_is_in_any_sheet = false;
	if (res.next())
	{
		strand_i_is_in_any_sheet = true;
	}
	return strand_i_is_in_any_sheet;
} //check_whether_strand_i_is_in_sheet


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
SandwichFeatures::get_current_strands_in_sheet(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	sh.sheet_id,\n"
	"	sh.segment_id,\n"
	"	bs.residue_begin,\n"
	"	bs.residue_end \n"
	"FROM\n"
	"	sheet as sh,\n"
	"	secondary_structure_segments AS bs\n"
	"WHERE\n"
	"	sh.segment_id = bs.segment_id \n"
	"	AND bs.dssp = 'E' \n"
	"	AND sh.struct_id = bs.struct_id \n"
	"	AND sh.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next()){
		Size sheet_id, segment_id,	residue_begin,	residue_end;
		res >> sheet_id >> segment_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(sheet_id, residue_begin, residue_end));
	}
	return all_strands;
} //get_current_strands_in_sheet


utility::vector1<Size>	
SandwichFeatures::get_distinct_sheet_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	distinct sheet_id\n"
	"FROM\n"
	"	sheet \n"
	"WHERE\n"
	"	struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> all_distinct_sheet_ids;
	while(res.next())
	{
		Size distinct_sheet_id;
		res >> distinct_sheet_id;
		all_distinct_sheet_ids.push_back(distinct_sheet_id);
	}
	return all_distinct_sheet_ids;
} //get_distinct_sheet_id


//get_max_sheet_id
Size
SandwichFeatures::get_max_sheet_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	max(sh.sheet_id) \n"
	"FROM\n"
	"	sheet AS sh \n"
	"WHERE\n"
	"	(sh.struct_id = ?) \n"
	"	AND (sh.sheet_id != 99999);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size max_sheet_id;
	while(res.next()){
		res >> max_sheet_id;
	}
	return max_sheet_id;
} //get_max_sheet_id


//update_sheet_id
Size
SandwichFeatures::update_sheet_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size new_sheet_id,
	Size old_sheet_id)
{
	string select_string =
	"UPDATE sheet set sheet_id = ?	"
	"WHERE\n"
	"	(sheet_id = ?) \n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,new_sheet_id);
	select_statement.bind(2,old_sheet_id);
	select_statement.bind(3,struct_id);
	basic::database::safely_write_to_database(select_statement);
	return 0;
} //update_sheet_id


//update_sheet_antiparallel
void	
SandwichFeatures::update_sheet_antiparallel(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sheet_id,
	string antiparallel)
{
	string select_string =
	"UPDATE sheet set sheet_antiparallel = ?	"
	"WHERE\n"
	"	(sheet_id = ?) \n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));

	select_statement.bind(1,antiparallel);
	select_statement.bind(2,sheet_id);
	select_statement.bind(3,struct_id);

	basic::database::safely_write_to_database(select_statement);
} //update_sheet_antiparallel


//get_chain_B_resNum
utility::vector1<SandwichFragment>
SandwichFeatures::get_chain_B_resNum(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	r.resNum \n"
	"FROM\n"
	"	sheet AS sh, \n"
	"	secondary_structure_segments AS bs, \n"
	"	residues AS r \n"
	"WHERE\n"
	"	(sh.sheet_id=1) \n"
	"	AND (bs.dssp = 'E') \n"
	"	AND (sh.segment_id=bs.segment_id) \n"
	"	AND (r.resNum >= bs.residue_begin AND r.resNum <= bs.residue_end) \n"
	"	AND (sh.struct_id = bs.struct_id) \n"
	"	AND (sh.struct_id = r.struct_id) \n"
	"	AND (sh.struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> return_chain_B_resNum;
	while(res.next())
	{
		Size chain_B_resNum;
		res >> chain_B_resNum;
		return_chain_B_resNum.push_back(SandwichFragment(chain_B_resNum));
	}
	return return_chain_B_resNum;
} //get_chain_B_resNum



//get_tag
string
SandwichFeatures::get_tag(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	tag \n"
	"FROM\n"
	"	structures AS s \n"
	"WHERE\n"
	"	s.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	string selected_tag;
	while (res.next())
	{
		res >> selected_tag;
	}
	return selected_tag;
} //get_tag

//get_num_strands_in_this_sheet
Size
SandwichFeatures::get_num_strands_in_this_sheet(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	string select_string =
	"SELECT\n"
	"	count(*) \n"
	"FROM\n"
	"	sheet AS sh \n"
	"WHERE\n"
	"	(sh.struct_id = ?) \n"
	"	AND (sh.sheet_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sheet_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_strands;
	while(res.next())
	{
		res >> num_strands;
	}
	return num_strands;
} //get_num_strands_in_this_sheet

//prepare_to_fill_sw_by_components
utility::vector1<SandwichFragment>
SandwichFeatures::prepare_to_fill_sw_by_components(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT sw_sh.sw_can_by_sh_id AS sw_can_by_sh_id, sh.sheet_id AS sheet_id, \n"
	"	bs.segment_id AS sw_by_components_bs_id, bs.residue_begin AS	residue_begin, \n"
	"	bs.residue_end AS residue_end \n"
	"FROM  secondary_structure_segments AS bs, structures AS s, sheet AS sh, sw_can_by_sh AS sw_sh\n"
	"WHERE (bs.struct_id = ?) AND (bs.struct_id = s.struct_id) AND (bs.struct_id = sh.struct_id) \n"
	"	AND (bs.struct_id = sw_sh.struct_id) \n"
	"	AND (bs.dssp = 'E') \n"
	"	AND (sw_sh.sheet_id = sh.sheet_id) AND (sh.segment_id = bs.segment_id);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));

	select_statement.bind(1,struct_id);

	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next())
	{
		Size sw_can_by_sh_id, sheet_id, sw_by_components_bs_id, residue_begin,	residue_end;
		res >> sw_can_by_sh_id >> sheet_id >> sw_by_components_bs_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(sw_can_by_sh_id, sheet_id, sw_by_components_bs_id, residue_begin, residue_end));
	}
	return all_strands;
} //prepare_to_fill_sw_by_components

Real
SandwichFeatures::abs (numeric::xyzVector<Real> vector)
{
	Real abs = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	return abs;
}


// I tried to slightly vary Nobu's recent nature paper
string
SandwichFeatures::check_LR ( // check L/R chirality of sidechain
	Pose & dssp_pose,
	Size residue_begin,
	Size residue_end)
{
	using core::id::NamedAtomID;
	using numeric::xyzVector;

	Size preceding_E = residue_begin-1;
	Size following_E = residue_end+1;

	xyzVector<Real> vector_u	=	dssp_pose.xyz(NamedAtomID("C", preceding_E)) - dssp_pose.xyz(NamedAtomID("N", preceding_E));
	xyzVector<Real> vector_v	=	dssp_pose.xyz(NamedAtomID("CA", following_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));
		// referred engineersphere.com/math/unit-vector-between-two-points.html
		// referred http://wiki.answers.com/Q/Can_dotproduct_of_two_vectors_be_negative
		// "If the dot-product is positive, then the angle between the two vectors is between 0 and 90 degrees. When the dot-product is negative, the angle is more than 90 degrees."

	numeric::xyzVector<Real> cross_product_u_v = cross_product(vector_u, vector_v);

	xyzVector<Real> vector_a_b; //initial

	if (dssp_pose.residue_type(following_E).name3() != "GLY")
	{
		vector_a_b = dssp_pose.xyz(NamedAtomID("CB", following_E)) - dssp_pose.xyz(NamedAtomID("CA", following_E));
	}
	else
	{	
		vector_a_b = dssp_pose.xyz(NamedAtomID("2HA", following_E)) - dssp_pose.xyz(NamedAtomID("CA", following_E));
	}
	
	Real dot_product_with_a_b = dot_product( cross_product_u_v, vector_a_b );
	Real cosine_theta = dot_product_with_a_b / (abs(cross_product_u_v))*(abs(vector_a_b));
		// referred http://www.mvps.org/directx/articles/math/dot/index.htm

	if (cosine_theta == 0) // theta = 90
	{
		return "B"; // borderline
	}
	if (cosine_theta <= 0.173 && cosine_theta > 0) // 80 <= theta < 90
	{
		return "BR"; // borderline, but could be classified as R;
	}
	else if (cosine_theta < 0 && cosine_theta >= -0.173) // 90 < theta <= 100
	{
		return "BL"; // borderline, but could be classified as L;
	}
	else if (cosine_theta > 0.173) // 0 <= theta < 80
	{
		return "R";
	}
	else // 100 < theta <= 180
	{
		return "L";
	}
} //check_LR

std::pair<string, string>
SandwichFeatures::check_PA( // parallel & anti-parallel
	Pose & dssp_pose,
	Size residue_begin,
	Size residue_end)
{
	using core::id::NamedAtomID;
	using numeric::xyzVector;

	Size preceding_E = residue_begin-1;
	Size following_E = residue_end+1;

	xyzVector<Real> preceding_E_vec_a_b; //initial
	xyzVector<Real> following_E_vec_a_b; //initial
	
	if (dssp_pose.residue_type(preceding_E).name3() != "GLY")
	{
		preceding_E_vec_a_b = dssp_pose.xyz(NamedAtomID("CB", preceding_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));
	}
	else
	{	
		preceding_E_vec_a_b = dssp_pose.xyz(NamedAtomID("2HA", preceding_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));
	}

	if (dssp_pose.residue_type(following_E).name3() != "GLY")
	{
		following_E_vec_a_b = dssp_pose.xyz(NamedAtomID("CB", following_E)) - dssp_pose.xyz(NamedAtomID("CA", following_E));
	}
	else
	{	
		following_E_vec_a_b = dssp_pose.xyz(NamedAtomID("2HA", following_E)) - dssp_pose.xyz(NamedAtomID("CA", following_E));
	}

	xyzVector<Real> vector_v	=	dssp_pose.xyz(NamedAtomID("CA", following_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));

	Real	dot_product_with_v_and_preceding_E = dot_product( preceding_E_vec_a_b, vector_v );
	Real	cosine_theta_between_v_and_preceding_E = dot_product_with_v_and_preceding_E / (abs(preceding_E_vec_a_b))*(abs(vector_v));

	Real	dot_product_with_v_and_following_E = dot_product( following_E_vec_a_b, vector_v );
	Real	cosine_theta_between_v_and_following_E = dot_product_with_v_and_following_E / (abs(following_E_vec_a_b))*(abs(vector_v));

	if (cosine_theta_between_v_and_preceding_E >=0 && cosine_theta_between_v_and_following_E >=0)
	{
		return std::make_pair("P", "P"); //"P_by_preceding_E", "P_by_following_E"
	}
	else if (cosine_theta_between_v_and_preceding_E >=0 && cosine_theta_between_v_and_following_E < 0)
	{
		return std::make_pair("P", "A");
	}
	else if (cosine_theta_between_v_and_preceding_E <0 && cosine_theta_between_v_and_following_E >= 0)
	{
		return std::make_pair("A", "P");
	}
	else
	{
		return std::make_pair("A", "A");
	}
} //check_PA


string
SandwichFeatures::check_heading_direction( // exclusively between preceding E and following E
	Pose & dssp_pose,
	Size residue_begin,
	Size residue_end,
	string check_N_to_C_direction_by_)
{
	using core::id::NamedAtomID;
	using numeric::xyzVector;

	Size preceding_E = residue_begin-1;
	Size following_E = residue_end+1;

	xyzVector<Real> vector_v	=	dssp_pose.xyz(NamedAtomID("CA", following_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));

	xyzVector<Real> preceding_E_vec_a_b; //initial
	xyzVector<Real> following_E_vec_a_b; //initial
	
	if (dssp_pose.residue_type(preceding_E).name3() != "GLY")
	{
		preceding_E_vec_a_b = dssp_pose.xyz(NamedAtomID("CB", preceding_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));
	}
	else
	{	
		preceding_E_vec_a_b = dssp_pose.xyz(NamedAtomID("2HA", preceding_E)) - dssp_pose.xyz(NamedAtomID("CA", preceding_E));
	}

	if (dssp_pose.residue_type(following_E).name3() != "GLY")
	{
		following_E_vec_a_b = dssp_pose.xyz(NamedAtomID("CB", following_E)) - dssp_pose.xyz(NamedAtomID("CA", following_E));
	}
	else
	{	
		following_E_vec_a_b = dssp_pose.xyz(NamedAtomID("2HA", following_E)) - dssp_pose.xyz(NamedAtomID("CA", following_E));
	}

	Real	dot_product_with_preceding_E_and_following_E = dot_product( preceding_E_vec_a_b, following_E_vec_a_b );
	Real	cosine_theta_between_sidechains = dot_product_with_preceding_E_and_following_E / (abs(preceding_E_vec_a_b))*(abs(following_E_vec_a_b));
		// referred http://www.mvps.org/directx/articles/math/dot/index.htm
		// cosine 80 = 0.173

	if (check_N_to_C_direction_by_ == "PE") // default
	{
		Real	dot_product_with_v_and_preceding_E = dot_product( preceding_E_vec_a_b, vector_v );
		Real	cosine_theta_between_v_and_preceding_E = dot_product_with_v_and_preceding_E / (abs(preceding_E_vec_a_b))*(abs(vector_v));

		if (cosine_theta_between_sidechains > 0) // 0 <= theta < 90
		{
			if (cosine_theta_between_v_and_preceding_E >=0) {	return "posi";	}
			else { return "nega";	}
		}
		else // 90 <= theta <= 180
		{
			if (cosine_theta_between_v_and_preceding_E >=0) {	return "away";	}
			else { return "meet";	}
		}
	}
	else if (check_N_to_C_direction_by_ == "FE")
	{
		Real	dot_product_with_v_and_following_E = dot_product( following_E_vec_a_b, vector_v );
		Real	cosine_theta_between_v_and_following_E = dot_product_with_v_and_following_E / (abs(following_E_vec_a_b))*(abs(vector_v));

		if (cosine_theta_between_sidechains > 0) // 0 <= theta < 90
		{
			if (cosine_theta_between_v_and_following_E >=0) {	return "posi";	}
			else { return "nega";	}
		}
		else // 90 <= theta <= 180
		{
			if (cosine_theta_between_v_and_following_E >=0) {	return "away";	}
			else { return "meet";	}
		}
	}
	else
	{
			TR.Info << "Exception:: check_N_to_C_direction_by should either PF or FE!!!" << endl;
		return "except";
	}
} //check_heading_direction

//find_sheet (assign a new strand into a sheet)
Size
SandwichFeatures::find_sheet(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell // if false, find parallel way
	)
{
	// seeing distances between 'O' of strand "i" and 'N' of strand "j"
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_CA_CA_0_0 = 0; // just initial assignment of value
			Real dis_CA_CA_1_1 = 0;
			Real dis_CA_CA_2_2 = 0;

			Real angle_C_O_N_0_0 = 0; // just initial assignment of value
			Real angle_C_O_N_1_1 = 0; // just initial assignment of value
			Real angle_C_O_N_2_2 = 0; // just initial assignment of value

			if (strand_i.get_size() == 2 || strand_j.get_size() == 2)
			{
				if (antiparalell)
				{
					if (i_resnum+1 > strand_i.get_end() || j_resnum-1 < strand_j.get_start())
					{
						continue; // I want to extract strand_pairs with only given ranges
					}
					dis_CA_CA_0_0 = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
						//	TR.Info << "distance between resnum("<< i_resnum << ")'s N and resnum(" << j_resnum << ")'s O = " << dis_N_O << endl;
					dis_CA_CA_1_1 = pose.residue(i_resnum+1).atom("CA").xyz().distance(pose.residue(j_resnum-1).atom("CA").xyz());
				}
				else // find a sheet in a parallel way
				{
					if (i_resnum+1 > strand_i.get_end() || j_resnum+1 > strand_j.get_end())
					{
						continue; // I want to extract strand_pairs with only given ranges
					}

					dis_CA_CA_0_0 = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
						//	TR.Info << "distance between resnum("<< i_resnum << ")'s N and resnum(" << j_resnum << ")'s O = " << dis_N_O << endl;
					dis_CA_CA_1_1 = pose.residue(i_resnum+1).atom("CA").xyz().distance(pose.residue(j_resnum+1).atom("CA").xyz());
				}

				if (dis_CA_CA_0_0 > 40)
				{
					return 999; // since these two strands are too distant to each other, there is virtually no chance to be sheet!
				}

				if (
					(dis_CA_CA_0_0 >= min_CA_CA_dis_ && dis_CA_CA_0_0 <= max_CA_CA_dis_)
					&& (dis_CA_CA_1_1 >= min_CA_CA_dis_ && dis_CA_CA_1_1 <= max_CA_CA_dis_)
					)
				{
					return 1; //  may have kinkness or not
				}

			}
			else // strand_i.get_size() >= 3 && strand_j.get_size() >= 3)
			{
				if (antiparalell)
				{
					if (i_resnum+2 > strand_i.get_end() || j_resnum-2 < strand_j.get_start())
					{
						continue; // I want to extract strand_pairs with only given ranges
					}

					vector<Real> dis_angle_inter_strands =
					cal_dis_angle_to_find_sheet(
												pose,
												i_resnum,
												i_resnum+1,
												i_resnum+2,
												j_resnum,
												j_resnum-1,
												j_resnum-2);

					dis_CA_CA_0_0 = dis_angle_inter_strands[0];
					dis_CA_CA_1_1 = dis_angle_inter_strands[1];
					dis_CA_CA_2_2 = dis_angle_inter_strands[2];

					angle_C_O_N_0_0 = dis_angle_inter_strands[3];
					angle_C_O_N_1_1 = dis_angle_inter_strands[4];
					angle_C_O_N_2_2 = dis_angle_inter_strands[5];

				}
				else // find a sheet in a parallel way
				{
					if (i_resnum+2 > strand_i.get_end() || j_resnum+2 > strand_j.get_end())
					{
						continue; // I want to extract strand_pairs with only given ranges
					}
					vector<Real> dis_angle_inter_strands =
					cal_dis_angle_to_find_sheet(
												pose,
												i_resnum,
												i_resnum+1,
												i_resnum+2,
												j_resnum,
												j_resnum+1,
												j_resnum+2);

					dis_CA_CA_0_0 = dis_angle_inter_strands[0];
					dis_CA_CA_1_1 = dis_angle_inter_strands[1];
					dis_CA_CA_2_2 = dis_angle_inter_strands[2];

					angle_C_O_N_0_0 = dis_angle_inter_strands[3];
					angle_C_O_N_1_1 = dis_angle_inter_strands[4];
					angle_C_O_N_2_2 = dis_angle_inter_strands[5];
				}

				if (dis_CA_CA_0_0 > 40)
				{
					return 999; // since these two strands are too distant to each other, there is virtually no chance to be sheet!
				}

				if (
					(dis_CA_CA_0_0 >= min_CA_CA_dis_ && dis_CA_CA_0_0 <= max_CA_CA_dis_)
					&& (dis_CA_CA_1_1 >= min_CA_CA_dis_ && dis_CA_CA_1_1 <= max_CA_CA_dis_)
					&& (dis_CA_CA_2_2 >= min_CA_CA_dis_ && dis_CA_CA_2_2 <= max_CA_CA_dis_)
					&& ((angle_C_O_N_0_0 >= min_C_O_N_angle_ && angle_C_O_N_2_2 >= min_C_O_N_angle_)
					|| (angle_C_O_N_1_1 >= min_C_O_N_angle_))
					)
				{
					return 1; //  may have kinkness or not, but these strands can be part of one sheet
				}
			} // strand_i.get_size() >= 3 && strand_j.get_size() >= 3)
		}
	}
	return 0; // these strands cannot be in one sheet
} //SandwichFeatures::find_sheet


//cal_dis_angle_to_find_sheet (assign a new strand into a sheet)
vector<Real>
SandwichFeatures::cal_dis_angle_to_find_sheet(
	Pose const & pose,
	Size res_i_0,
	Size res_i_1,
  	Size res_i_2,
	Size res_j_0,
	Size res_j_1,
  	Size res_j_2)
{
	Real dis_CA_CA_0_0 = pose.residue(res_i_0).atom("CA").xyz().distance(pose.residue(res_j_0).atom("CA").xyz());

	Vector const& first_0_xyz    ( pose.residue(res_i_0).xyz("C") );
	Vector const& middle_0_xyz   ( pose.residue(res_i_0).xyz("O") );
	Vector const& third_0_xyz    ( pose.residue(res_j_0).xyz("N") );

	Real angle_C_O_N_0_0 = numeric::angle_degrees(first_0_xyz, middle_0_xyz, third_0_xyz);

	Real dis_CA_CA_1_1 = pose.residue(res_i_1).atom("CA").xyz().distance(pose.residue(res_j_1).atom("CA").xyz());

	Vector const& first_1_xyz    ( pose.residue(res_i_1).xyz("C") );
	Vector const& middle_1_xyz   ( pose.residue(res_i_1).xyz("O") );
	Vector const& third_1_xyz    ( pose.residue(res_j_1).xyz("N") );

	Real angle_C_O_N_1_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);

	Real dis_CA_CA_2_2 = pose.residue(res_i_2).atom("CA").xyz().distance(pose.residue(res_j_2).atom("CA").xyz());

	Vector const& first_2_xyz    ( pose.residue(res_i_2).xyz("C") );
	Vector const& middle_2_xyz   ( pose.residue(res_i_2).xyz("O") );
	Vector const& third_2_xyz    ( pose.residue(res_j_2).xyz("N") );

	Real angle_C_O_N_2_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);

	vector<Real> dis_angle_inter_strands;

	dis_angle_inter_strands.push_back ( dis_CA_CA_0_0 );
	dis_angle_inter_strands.push_back ( dis_CA_CA_1_1 );
	dis_angle_inter_strands.push_back ( dis_CA_CA_2_2 );
	dis_angle_inter_strands.push_back ( angle_C_O_N_0_0 );
	dis_angle_inter_strands.push_back ( angle_C_O_N_1_1 );
	dis_angle_inter_strands.push_back ( angle_C_O_N_2_2 );

	return dis_angle_inter_strands;
} //cal_dis_angle_to_find_sheet

bool
SandwichFeatures::see_whether_sheets_can_be_combined(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size i_sheet,
	Size j_sheet)
{
	utility::vector1<SandwichFragment> strands_from_i = get_full_strands_from_sheet(struct_id, db_session, i_sheet);
	utility::vector1<SandwichFragment> strands_from_j = get_full_strands_from_sheet(struct_id, db_session, j_sheet);

	Size return_of_find_sheet_antiparallel(0); // temporary 'false' designation
	Size return_of_find_sheet_parallel(0); // temporary 'false' designation

	for(Size i=1; i<=strands_from_i.size(); ++i)
	{
		for(Size j=1; j<=strands_from_j.size(); ++j)
		{
			SandwichFragment temp_strand_i(strands_from_i[i].get_start(), strands_from_i[i].get_end());
			SandwichFragment temp_strand_j(strands_from_j[j].get_start(), strands_from_j[j].get_end());

			return_of_find_sheet_antiparallel = find_sheet (pose, temp_strand_i, temp_strand_j, true);

			if (return_of_find_sheet_antiparallel == 999)
			{
				break; // too distance strands
			}

			if (!return_of_find_sheet_antiparallel)
			{
				return_of_find_sheet_parallel = find_sheet (pose, temp_strand_i, temp_strand_j, false);
			}

			if (return_of_find_sheet_parallel == 999)
			{
				break;
			}

			if (return_of_find_sheet_antiparallel || return_of_find_sheet_parallel)
			{
				return true; // these two sheets should be combined
			}
		}
	}
	return false; // these two sheets should not be combined
} //SandwichFeatures::see_whether_sheets_can_be_combined


bool
SandwichFeatures::change_sheet_id_if_possible( // combine_sheets_if_possible
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose)
{
	bool	sheet_id_changed = false; // don't repeat change_sheet_id_if_possible
	Size max_sheet_id = get_max_sheet_id(struct_id, db_session);
	for(Size i=1; i<=max_sheet_id-1; ++i)
	{
		for(Size j=i+1; j<=max_sheet_id; ++j)
		{
			bool sheets_can_be_combined = see_whether_sheets_can_be_combined(
																			 struct_id,
																			 db_session,
																			 pose,
																			 i,
																			 j);
			if (sheets_can_be_combined)
			{
				if(i<j)
				{
					update_sheet_id(
									struct_id,
									db_session,
									i, //new_sheet_id
									j); //old_sheet_id
				}
				else
				{
					update_sheet_id(
									struct_id,
									db_session,
									j, //new_sheet_id
									i); //old_sheet_id
				}
				sheet_id_changed = true; // repeat change_sheet_id_if_possible
			}
		}
	}
	return sheet_id_changed;
}	//SandwichFeatures::change_sheet_id_if_possible


// In order to be antiparallel sheet, all strands in the sheet should be anti-parallel to each other
string	
SandwichFeatures::see_whether_sheet_is_antiparallel(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size i_sheet)
{
	utility::vector1<SandwichFragment> strands_from_i = get_full_strands_from_sheet(struct_id, db_session, i_sheet); // struct_id, db_session, sheet_id

	// <begin> get central residues
	vector<Real> vec_cen_resnum_i; // array of central residues
	vector<Size> vec_unrepresentative_strand;
	
	for(Size i=1; i<=strands_from_i.size(); ++i)
	{
		bool this_strand_represent_a_terminal = can_this_strand_represent_a_terminal (pose, strands_from_i, i); //pose, full_strands_from_this_sheet, each strand
		if (!this_strand_represent_a_terminal)
		{
			vec_cen_resnum_i.push_back(-999); // won't be used, but needed
			vec_unrepresentative_strand.push_back(i);
			continue;
		}// I ignore a unrepresentative strand
		
		Real to_be_rounded_i = (strands_from_i[i].get_start() + strands_from_i[i].get_end())/(2.0);
		Size cen_resnum_i = round(to_be_rounded_i);
		vec_cen_resnum_i.push_back(cen_resnum_i);
	}
	// <end> get central residues

	// <begin> get array of sum of all distances from central residues
	vector<Real> vec_dis_sum_from_cen_resnum_i; // array of sum of all distances from central residues
	for(Size ii=0; ii<=strands_from_i.size()-1; ++ii)
	{
		Real dis_from_ii = 0;
		
		if (vec_cen_resnum_i[ii] == -999) // I ignore a core strand
		{
			dis_from_ii = -999;
			vec_dis_sum_from_cen_resnum_i.push_back(dis_from_ii);
			continue;
		}

		for (Size jj=0; jj<=strands_from_i.size()-1; ++jj)
		{				
			if (vec_cen_resnum_i[jj] == -999) // I ignore a core strand
			{
				continue;
			}
			//Real dis_CA_CA = pose.residue(vec_cen_resnum_i[ii]).atom("CA").xyz().distance(pose.residue(vec_cen_resnum_i[jj]).atom("CA").xyz());
			Real dis_CA_CA	=	pose.residue(static_cast<Size>(vec_cen_resnum_i[ii])).atom("CA").xyz().distance(pose.residue(static_cast<Size>(vec_cen_resnum_i[jj])).atom("CA").xyz());
			dis_from_ii = dis_from_ii + dis_CA_CA;
		}
		vec_dis_sum_from_cen_resnum_i.push_back(dis_from_ii);
	}
	// <end> get array of sum of all distances from central residues


	// <begin> get the largest distance and its residue index
	Size res_index_having_the_largest_dis;
	Real largest_dis = -99;

	for(Size ii=0; ii<=strands_from_i.size()-1; ++ii)
	{
		if (vec_cen_resnum_i[ii] == -999) 
		{
			continue; // unrepresentative_strands
		}
		if (largest_dis < vec_dis_sum_from_cen_resnum_i[ii])
		{
			largest_dis = vec_dis_sum_from_cen_resnum_i[ii];
			res_index_having_the_largest_dis = ii;
		}
	}
	// <end> get largest distance and its residue index

	Size former_res_index_nearest_strand = res_index_having_the_largest_dis; //just for the first step of 'while' loop
	Size former_res_index_having_the_largest_dis = res_index_having_the_largest_dis; //just for the first step of 'while' loop
	Size exit_condition = strands_from_i.size()-1-vec_unrepresentative_strand.size();
	Size count_anti = 0;
	while (count_anti < exit_condition)
	{
		SandwichFragment temp_strand_i(strands_from_i[former_res_index_nearest_strand+1].get_start(), strands_from_i[former_res_index_nearest_strand+1].get_end());
			
		Real shortest_dis_inter_strand = 999;
		Size res_index_nearest_strand = 999;

		//<begin> search the nearest strand from the current strand
		for(Size ii=0; ii<=strands_from_i.size()-1; ++ii)
		{
			if (ii == former_res_index_nearest_strand || ii == former_res_index_having_the_largest_dis)
			{
				continue;
			}
			if (vec_cen_resnum_i[ii] == -999) // I ignore a core strand
			{
				continue; 
			}
			SandwichFragment temp_strand_j(strands_from_i[ii+1].get_start(), strands_from_i[ii+1].get_end());
			Real inter_strand_avg_dis = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
			if (inter_strand_avg_dis < shortest_dis_inter_strand)
			{
				shortest_dis_inter_strand = inter_strand_avg_dis;
				res_index_nearest_strand = ii;
			}
		}
		//<end> search the nearest strand from the current strand
		
		SandwichFragment temp_strand_j(strands_from_i[res_index_nearest_strand+1].get_start(), strands_from_i[res_index_nearest_strand+1].get_end());
		
		Real return_of_check_sw_by_dis_anti = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, true);
		if ( return_of_check_sw_by_dis_anti != -99 )
		{
			count_anti++;
			former_res_index_having_the_largest_dis = former_res_index_nearest_strand;
			former_res_index_nearest_strand = res_index_nearest_strand;
		}
		else
		{
			return "P_or_mix";	//	sheet is parallel or mixed form
		}
	}
	return "A"; //	sheet_is_antiparallel
} //SandwichFeatures::see_whether_sheet_is_antiparallel


// check whether these sheets are too close, the closeness is checked for every possible distances
bool
SandwichFeatures::check_strand_too_closeness(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j)
{
	// check anti-parallel sheet distance
	// first, check the shortest distance between the two strand_pairs
	// seeing distances between 'CA' of strand "i" and 'CA' of strand "j"
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
			if (dis_CA_CA < min_inter_sheet_dis_CA_CA_) // these two pair of strands are too close
			{
				return true;
			}
		}
	}
	return false; // OK, these two strand_pairs are farther enough according to this check_strand_too_closeness
} //SandwichFeatures::check_strand_too_closeness


// check whether these sheets are too close, the closeness is checked for every possible distances
Real
SandwichFeatures::get_avg_dis_strands(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j)
{
	Real sum_dis_CA_CA = 0;
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
			sum_dis_CA_CA = sum_dis_CA_CA + dis_CA_CA;
		}
	}
	return sum_dis_CA_CA/(strand_i.get_size()*strand_j.get_size());
} //SandwichFeatures::get_avg_dis_strands



Real
SandwichFeatures::get_avg_dis_CA_CA(
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
		//TR.Info << "dis_CA_CA_0: " << dis_CA_CA_0 << endl;
	if (dis_CA_CA_0 > 40)
	{
		return -999; // these sheets will not be sandwich ever, since these two sheets are too distant!
	}
	if (dis_CA_CA_0 < min_sheet_dis_ || dis_CA_CA_0 > max_sheet_dis_)
	{
		return -99;
	}

	Real dis_CA_CA_1 = pose.residue(i_resnum_1).atom("CA").xyz().distance(pose.residue(j_resnum_1).atom("CA").xyz());
		//TR.Info << "dis_CA_CA_1: " << dis_CA_CA_1 << endl;
	if (dis_CA_CA_1 < min_sheet_dis_ || dis_CA_CA_1 > max_sheet_dis_)
	{
		return -99;
	}

	Real dis_CA_CA_2 = pose.residue(i_resnum_2).atom("CA").xyz().distance(pose.residue(j_resnum_2).atom("CA").xyz());
		//TR.Info << "dis_CA_CA_2: " << dis_CA_CA_2 << endl;
	if (dis_CA_CA_2 < min_sheet_dis_ || dis_CA_CA_2 > max_sheet_dis_)
	{
		return -99;
	}

	Real dis_CA_CA_3 = pose.residue(i_resnum_3).atom("CA").xyz().distance(pose.residue(j_resnum_3).atom("CA").xyz());
		//TR.Info << "dis_CA_CA_3: " << dis_CA_CA_3 << endl;
	if (dis_CA_CA_3 < min_sheet_dis_ || dis_CA_CA_3 > max_sheet_dis_)
	{
		return -99;
	}

	Real avg_dis_CA_CA = (dis_CA_CA_0 + dis_CA_CA_1 + dis_CA_CA_2 + dis_CA_CA_3)/4;
	return avg_dis_CA_CA;
} // SandwichFeatures::get_avg_dis_CA_CA



Real
SandwichFeatures::check_sw_by_dis(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell // if false, find parallel way
	)
{
	Size i_resnum_1 = 0; // just initial temporary assignment
	Size j_resnum_1 = 0;

	Size i_resnum_2 = 0;
	Size j_resnum_2 = 0;

	Size i_resnum_3 = 0;
	Size j_resnum_3 = 0;

	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;

			if (antiparalell)
			{
				i_resnum_1 = i_resnum+1;
				j_resnum_1 = j_resnum-1;

				i_resnum_2 = i_resnum+2;
				j_resnum_2 = j_resnum-2;

				i_resnum_3 = i_resnum+3;
				j_resnum_3 = j_resnum-3;

				if (j_resnum_3 <= 0 
					|| i_resnum_3 > pose.total_residue() 
					|| j_resnum_3 > pose.total_residue()) // sometimes, j_resnum_3 becomes 18446744073709551615 where it should be -1
				{
					continue;
				}
			}

			else // paralell
			{
				i_resnum_1 = i_resnum+1;
				j_resnum_1 = j_resnum+1;

				i_resnum_2 = i_resnum+2;
				j_resnum_2 = j_resnum+2;

				i_resnum_3 = i_resnum+3;
				j_resnum_3 = j_resnum+3;

				if (i_resnum_3 > pose.total_residue() || j_resnum_3 > pose.total_residue())
				{
					continue;
				}
			}

			Real avg_dis_CA_CA = get_avg_dis_CA_CA(pose, i_resnum,	i_resnum_1, i_resnum_2, i_resnum_3, j_resnum, j_resnum_1, j_resnum_2, j_resnum_3);

			if (avg_dis_CA_CA == -999)
			{
				break; // these sheets will not be sandwich ever, since these two sheets are too distant!
			}

			if (avg_dis_CA_CA == -99) // dis_CA_CA_x < min_sheet_dis_ || dis_CA_CA_x > max_sheet_dis_
			{
				continue;
			}
			return avg_dis_CA_CA;
		} //for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
	} //for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	return -99; // these sheets are not sandwich with these strands
} //SandwichFeatures::check_sw_by_dis

Size
SandwichFeatures::round(
	Real x)
{
	Size rounded = static_cast <Size> (floor(x+.5));
	return rounded;
} //round

//can this strand represent a terminal strand for inter-sheet angle calculation?
bool
SandwichFeatures::can_this_strand_represent_a_terminal(
	Pose const & pose,
	utility::vector1<SandwichFragment> strands_from_sheet_i,
	Size current_strand_id_as_i)
{
	///////////////////// <begin> identify the closest strand from current_strand_id_as_i
	vector<Real> vec_inter_strand_avg_dis;
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i)
	{
		if (i == current_strand_id_as_i)
		{
			vec_inter_strand_avg_dis.push_back(999);
			continue; // we are calculating a distance between same strands
		}
		SandwichFragment temp_strand_i(strands_from_sheet_i[i].get_start(), strands_from_sheet_i[i].get_end());
		SandwichFragment temp_strand_j(strands_from_sheet_i[current_strand_id_as_i].get_start(), strands_from_sheet_i[current_strand_id_as_i].get_end());

		Real inter_strand_avg_dis = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
		vec_inter_strand_avg_dis.push_back(inter_strand_avg_dis);
	}

	Size array_size = vec_inter_strand_avg_dis.size();

	Real min_inter_strand_avg_dis = 9999;
	Size index_having_min_dis = 0;

	for(Size i=0; i<=array_size-1; ++i)
	{
		if (min_inter_strand_avg_dis > vec_inter_strand_avg_dis[i])
		{
			min_inter_strand_avg_dis = vec_inter_strand_avg_dis[i];
			index_having_min_dis = i+1; // index of vec_inter_strand_avg_dis starts with 0 while index of strands_from_sheet_i starts with 1
		}
	}
	///////////////////// <end> identify the closest strand from current_strand_id_as_i
	
	Real to_be_rounded_i = (strands_from_sheet_i[index_having_min_dis].get_start() + strands_from_sheet_i[index_having_min_dis].get_end())/(2.0);
	Size cen_resnum_i = round(to_be_rounded_i);

	for(Size strand_i_res=strands_from_sheet_i[current_strand_id_as_i].get_start(); 
		strand_i_res <= strands_from_sheet_i[current_strand_id_as_i].get_end(); 
		strand_i_res++)
	{
		Real dis_CA_CA = pose.residue(strand_i_res).atom("CA").xyz().distance(pose.residue(cen_resnum_i).atom("CA").xyz());
		if (dis_CA_CA >= min_CA_CA_dis_ && dis_CA_CA <= max_CA_CA_dis_)
		{
			return true; // this strand can represent a terminal strand for inter-sheet angle calculation
		}
	}
	return false; // this strand cannot represent a terminal strand for inter-sheet angle calculation
}//SandwichFeatures::can_this_strand_represent_a_terminal

//check whether this sheet is constituted with 2 (or less) residues long only
bool
SandwichFeatures::check_whether_this_sheet_is_too_short(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_i)
{
		//TR.Info << "sheet_i: " << sheet_i << endl;
	utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_i);
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i)
	{
		if (strands_from_sheet_i[i].get_size() > 2)
		{
			return false; // no, this sheet is not too short
		}
	}
	return true; // yes, this sheet is too short
} //check_whether_this_sheet_is_too_short

//return_two_central_residues
std::pair<Size, Size>
SandwichFeatures::get_two_central_residues(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i)
{
		//TR.Info << "sheet_i: " << sheet_i << endl;
	utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_i);

	// strands_from_sheet_i
	// get central residue numbers
	vector<Size> cen_res_arr_sheet_i;
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i)
	{
		if (strands_from_sheet_i[i].get_size() <= 2)
		{
				TR.Info << "strands_from_sheet_i[i].get_size() <= 2" << endl;
			continue;
		}
		bool this_strand_represent_a_terminal = can_this_strand_represent_a_terminal (pose, strands_from_sheet_i, i); // pose. strands_from_sheet_i, current_strand_id_as_i

		if (!this_strand_represent_a_terminal)
		{
				TR.Info << "!this_strand_represent_a_terminal" << endl;
			continue;
		}

		Real to_be_rounded_i = (strands_from_sheet_i[i].get_start() + strands_from_sheet_i[i].get_end())/(2.0);
		Size cen_resnum_i = round(to_be_rounded_i);
		cen_res_arr_sheet_i.push_back(cen_resnum_i);
	}
	Size array_size = cen_res_arr_sheet_i.size();

	// get sum of distances
	vector<Real> sum_dis_array_i;
	for(Size i=0; i<=array_size-1; ++i)
	{
		Real sum_dis_i_and_j = 0;
		for(Size j=0; (j<=array_size-1); ++j)
		{
			if (i == j)
			{
				continue;
			}
			Real dis_i_and_j = pose.residue(cen_res_arr_sheet_i[i]).atom("CA").xyz().distance(pose.residue(cen_res_arr_sheet_i[j]).atom("CA").xyz());
			sum_dis_i_and_j = sum_dis_i_and_j + dis_i_and_j;
		}
		sum_dis_array_i.push_back(sum_dis_i_and_j);
	}

	// pick two terminal central residue numbers

	// terminal central residue 1
	Real max_i_1 = -99;
	Size index_terminal_cen_res_pos_1 = 0;
	for(Size i=0; i<=array_size-1; ++i)
	{
		if (max_i_1 < sum_dis_array_i[i])
		{
			max_i_1 = sum_dis_array_i[i];
			index_terminal_cen_res_pos_1 = i;
		}
	}

	// terminal central residue 2
	Real max_i_2 = -99;
	Size index_terminal_cen_res_pos_2 = 0;
	for(Size i=0; i<=array_size-1; ++i)
	{
		if (i == index_terminal_cen_res_pos_1)
		{
			continue;
		}
		if (max_i_2 < sum_dis_array_i[i])
		{
			max_i_2 = sum_dis_array_i[i];
			index_terminal_cen_res_pos_2 = i;
		}
	}
	Size terminal_cen_res_pos_1 = cen_res_arr_sheet_i[index_terminal_cen_res_pos_1];	// index of the first central_residue
	Size terminal_cen_res_pos_2 = cen_res_arr_sheet_i[index_terminal_cen_res_pos_2];	// index of the second central_residue
	return std::make_pair(terminal_cen_res_pos_1, terminal_cen_res_pos_2);
} //get_two_central_residues

Real
SandwichFeatures::get_shortest_among_4(
	Real arr_dis_inter_sheet[])
{
	Real temp_shortest_dis = 9999;
	for(Size i=0; i<=3; ++i)
	{
		if (temp_shortest_dis > arr_dis_inter_sheet[i])
		{
			temp_shortest_dis = arr_dis_inter_sheet[i];
		}
	}
	return temp_shortest_dis;
} //SandwichFeatures::get_shortest_among_4 (simple one with just four parameters)

bool
SandwichFeatures::judge_facing(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i,
	Size sheet_j)
{
	bool this_strand_is_too_short = check_whether_this_sheet_is_too_short(
		struct_id,
		db_session,
		sheet_i);

	if (this_strand_is_too_short)
	{
		return false; // I can't choose two central residues since this sheet is constituted with 2 residues long strands only
	}

	this_strand_is_too_short = check_whether_this_sheet_is_too_short(
		struct_id,
		db_session,
		sheet_j);

	if (this_strand_is_too_short)
	{
		return false; // I can't choose two central residues since this sheet is constituted with 2 residues long strands only
	}

	// <begin> identify four terminal central residues
	std::pair<Size, Size>
	two_central_residues =	get_two_central_residues(
		struct_id,
		db_session,
		pose,
		sheet_i);
		
	Size i_ter_cen_1 = two_central_residues.first;
	Size i_ter_cen_2 = two_central_residues.second;
		///TR.Info << "i_ter_cen_1: " << i_ter_cen_1 << endl;
		//TR.Info << "i_ter_cen_2: " << i_ter_cen_2 << endl;

	two_central_residues =	get_two_central_residues(
		struct_id,
		db_session,
		pose,
		sheet_j);
		
	Size j_ter_cen_1 = two_central_residues.first;
	Size j_ter_cen_2 = two_central_residues.second;
		//TR.Info << "j_ter_cen_1: " << j_ter_cen_1 << endl;
		//TR.Info << "j_ter_cen_2: " << j_ter_cen_2 << endl;
		
	Real arr_dis_inter_sheet [4];
	arr_dis_inter_sheet[0] = pose.residue(i_ter_cen_1).atom("CA").xyz().distance(pose.residue(j_ter_cen_1).atom("CA").xyz());
	arr_dis_inter_sheet[1] = pose.residue(i_ter_cen_1).atom("CA").xyz().distance(pose.residue(j_ter_cen_2).atom("CA").xyz());
	arr_dis_inter_sheet[2] = pose.residue(i_ter_cen_2).atom("CA").xyz().distance(pose.residue(j_ter_cen_1).atom("CA").xyz());
	arr_dis_inter_sheet[3] = pose.residue(i_ter_cen_2).atom("CA").xyz().distance(pose.residue(j_ter_cen_2).atom("CA").xyz());

	Real
	shortest_dis_inter_sheet = get_shortest_among_4(arr_dis_inter_sheet);

	Real angle_1 = 999.9; //temp
	Real angle_2 = 999.9;
	Real torsion_i_j = 999.9;

	if (shortest_dis_inter_sheet == arr_dis_inter_sheet[0])
	{
		Vector const& first_1_xyz    ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(j_ter_cen_1).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(j_ter_cen_2).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(j_ter_cen_1).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(i_ter_cen_1).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);

		Vector const& first_xyz    ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(j_ter_cen_1).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(j_ter_cen_2).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}

	else if (shortest_dis_inter_sheet == arr_dis_inter_sheet[1])
	{
		Vector const& first_1_xyz    ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(j_ter_cen_2).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(j_ter_cen_1).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(j_ter_cen_2).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(i_ter_cen_1).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);

		Vector const& first_xyz    ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(j_ter_cen_2).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(j_ter_cen_1).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}

	else if (shortest_dis_inter_sheet == arr_dis_inter_sheet[2])
	{
		Vector const& first_1_xyz    ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(j_ter_cen_1).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(j_ter_cen_2).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(j_ter_cen_1).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(i_ter_cen_2).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);

		Vector const& first_xyz    ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(j_ter_cen_1).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(j_ter_cen_2).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}

	else // (shortest_dis_inter_sheet == arr_dis_inter_sheet[3])
	{
		Vector const& first_1_xyz    ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& middle_1_xyz   ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& third_1_xyz    ( pose.residue(j_ter_cen_2).xyz("CA") );

		Vector const& first_2_xyz    ( pose.residue(j_ter_cen_1).xyz("CA") );
		Vector const& middle_2_xyz   ( pose.residue(j_ter_cen_2).xyz("CA") );
		Vector const& third_2_xyz    ( pose.residue(i_ter_cen_2).xyz("CA") );

		angle_1 = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		angle_2 = numeric::angle_degrees(first_2_xyz, middle_2_xyz, third_2_xyz);

		Vector const& first_xyz    ( pose.residue(i_ter_cen_1).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(i_ter_cen_2).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(j_ter_cen_2).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(j_ter_cen_1).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}
		//TR.Info << "angle_1: " << angle_1 << endl;
		//TR.Info << "angle_2: " << angle_2 << endl;
		//TR.Info << "torsion_i_j: " << torsion_i_j << endl;

	if ((angle_1 > min_sheet_angle_) && (angle_1 < max_sheet_angle_) && (angle_2 > min_sheet_angle_) && (angle_2 < max_sheet_angle_) && (torsion_i_j > min_sheet_torsion_cen_res_) && (torsion_i_j < max_sheet_torsion_cen_res_))
	{
		return true; // these two strand_pairs face to each other properly, so constitute a sandwich
	}

	else
	{
		return false; // these two strand_pairs are linear or do not face to each other properly!
	}

} //SandwichFeatures::judge_facing


Size
SandwichFeatures::write_to_sheet (
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_PK_id_counter,
	Size sheet_id,
	Size segment_id)
{
	string sheet_insert_i =
	"INSERT INTO sheet (sheet_PK_id, sheet_id, struct_id, segment_id)  VALUES (?,?,?,?);";
	statement sheet_insert_i_stmt(basic::database::safely_prepare_statement(sheet_insert_i,db_session));

	sheet_insert_i_stmt.bind(1,	sheet_PK_id_counter);
	sheet_insert_i_stmt.bind(2,	sheet_id);
	sheet_insert_i_stmt.bind(3,	struct_id);
	sheet_insert_i_stmt.bind(4,	segment_id);
	basic::database::safely_write_to_database(sheet_insert_i_stmt);
	return 0;
} //write_to_sheet


Size
SandwichFeatures::write_to_sw_can_by_sh	(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id_counter,
	Size sheet_id,
	Size num_strands_from_sheet)
{
	string sw_can_by_sh_insert =
	 "INSERT INTO sw_can_by_sh (struct_id, sw_can_by_sh_PK_id, tag, sw_can_by_sh_id, sheet_id, strand_num)  VALUES (?,?,?,?,?,?);";

	statement sw_can_by_sh_insert_stmt(basic::database::safely_prepare_statement(sw_can_by_sh_insert,	db_session));

	sw_can_by_sh_insert_stmt.bind(1,	struct_id);
	sw_can_by_sh_insert_stmt.bind(2,	sw_can_by_sh_PK_id_counter);
	sw_can_by_sh_insert_stmt.bind(3,	tag);
	sw_can_by_sh_insert_stmt.bind(4,	sw_can_by_sh_id_counter);
	sw_can_by_sh_insert_stmt.bind(5,	sheet_id);
	sw_can_by_sh_insert_stmt.bind(6,	num_strands_from_sheet); // number of strands of sheet
	basic::database::safely_write_to_database(sw_can_by_sh_insert_stmt);
	return 0;
} //write_to_sw_can_by_sh

	
// <begin> see_whether_strand_is_at_edge
bool
SandwichFeatures::see_whether_strand_is_at_edge	(
	Pose const & pose,
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id,
	string sheet_antiparallel,
	Size residue_begin,
	Size residue_end)
{
// <begin> see two closest strands from temp_strand_i	
	utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_id);

	if (strands_from_sheet_i.size() < 3)
	{
		return true; // this strand is at edge
	}
	
	SandwichFragment temp_strand_i(residue_begin, residue_end);
	vector<Real> vec_inter_strand_avg_dis;
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i)
	{
		SandwichFragment temp_strand_j(strands_from_sheet_i[i].get_start(), strands_from_sheet_i[i].get_end());
		Real inter_strand_avg_dis = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
		vec_inter_strand_avg_dis.push_back(inter_strand_avg_dis);
	}
	
	Size array_size = vec_inter_strand_avg_dis.size();
	
	// <begin> exclude self-strand
	Real min_inter_strand_avg_dis = 9999;
	Size index_having_self_strand = 0;
	for(Size i=0; i<=array_size-1; ++i)
	{
		if (min_inter_strand_avg_dis > vec_inter_strand_avg_dis[i])
		{
			min_inter_strand_avg_dis = vec_inter_strand_avg_dis[i];
			index_having_self_strand = i+1; // index of vec_inter_strand_avg_dis starts with 0 while index of strands_from_sheet_i starts with 1
		}
	}
	// <end> exclude self-strand
	
	// <begin> find the closest strand
	min_inter_strand_avg_dis = 9999;
	Size index_having_min_dis = 0;
	for(Size i=0; i<=array_size-1; ++i)
	{
		if (i != index_having_self_strand-1 && (min_inter_strand_avg_dis > vec_inter_strand_avg_dis[i]))
		{
			min_inter_strand_avg_dis = vec_inter_strand_avg_dis[i];
			index_having_min_dis = i+1; // index of vec_inter_strand_avg_dis starts with 0 while index of strands_from_sheet_i starts with 1
		}
	}
	// <end> find the closest strand
	
	// <begin> find the 2nd closest strand
	min_inter_strand_avg_dis = 9999;
	Size index_having_second_min_dis = 0;
	for(Size i=0; i<=array_size-1; ++i)
	{
		if (i != index_having_self_strand-1 && i != index_having_min_dis-1 && (min_inter_strand_avg_dis > vec_inter_strand_avg_dis[i]))
		{
			min_inter_strand_avg_dis = vec_inter_strand_avg_dis[i];
			index_having_second_min_dis = i+1; // index of vec_inter_strand_avg_dis starts with 0 while index of strands_from_sheet_i starts with 1
		}
	}
	// <end> find the 2nd closest strand
// <end> see two closest strands from temp_strand_i
	
	SandwichFragment temp_strand_1(strands_from_sheet_i[index_having_min_dis].get_start(), strands_from_sheet_i[index_having_min_dis].get_end());
	SandwichFragment temp_strand_2(strands_from_sheet_i[index_having_second_min_dis].get_start(), strands_from_sheet_i[index_having_second_min_dis].get_end());
	
	if (sheet_antiparallel == "A")
	{
		bool	return_of_find_sheet_antiparallel_1 = find_sheet (pose, temp_strand_i, temp_strand_1, true);
		bool	return_of_find_sheet_antiparallel_2 = find_sheet (pose, temp_strand_i, temp_strand_2, true);		
		
		if (return_of_find_sheet_antiparallel_1 && return_of_find_sheet_antiparallel_2)
		{
			return false; // this strand is at core
		}
		return true; // this strand is at edge
	}	
	else // this sheet could be parallel or combination of anti-parallel and parallel
	{
		bool	return_of_find_sheet_antiparallel_1 = find_sheet (pose, temp_strand_i, temp_strand_1, false);
		bool	return_of_find_sheet_antiparallel_2 = find_sheet (pose, temp_strand_i, temp_strand_2, false);		
		
		if (return_of_find_sheet_antiparallel_1 && return_of_find_sheet_antiparallel_2)
		{
			return false; // this strand is at core		
		}
		
		return_of_find_sheet_antiparallel_1 = find_sheet (pose, temp_strand_i, temp_strand_1, true);
		return_of_find_sheet_antiparallel_2 = find_sheet (pose, temp_strand_i, temp_strand_2, false);		
		
		if (return_of_find_sheet_antiparallel_1 && return_of_find_sheet_antiparallel_2)
		{
			return false; // this strand is at core		
		}
		
		return_of_find_sheet_antiparallel_1 = find_sheet (pose, temp_strand_i, temp_strand_1, false);
		return_of_find_sheet_antiparallel_2 = find_sheet (pose, temp_strand_i, temp_strand_2, true);		
		
		if (return_of_find_sheet_antiparallel_1 && return_of_find_sheet_antiparallel_2)
		{
			return false; // this strand is at core		
		}

		return true; // this strand is at edge
	}
}
// <end> see_whether_strand_is_at_edge
	

void	
SandwichFeatures::fill_sw_by_components	(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_by_components_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size sheet_id,
	string sheet_antiparallel,
	Size sw_by_components_bs_id,
	Size strand_is_at_edge,
	Size residue_begin,
	Size residue_end)
{
	string insert =
	"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, sheet_id, sheet_antiparallel, sw_by_components_bs_id, sw_by_components_bs_edge, residue_begin, residue_end)  VALUES (?,?,?,?,?,?,?,?,?,?);";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	sw_by_components_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	sw_can_by_sh_id);
	insert_stmt.bind(5,	sheet_id);
	insert_stmt.bind(6,	sheet_antiparallel);
	insert_stmt.bind(7,	sw_by_components_bs_id); //bs_id
	insert_stmt.bind(8,	strand_is_at_edge);
	insert_stmt.bind(9,	residue_begin);
	insert_stmt.bind(10,	residue_end);
	basic::database::safely_write_to_database(insert_stmt);
} //fill_sw_by_components


//get_distinct(sw_can_by_sh_id)
utility::vector1<Size>
SandwichFeatures::get_vec_sw_can_by_sh_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	distinct(sw_can_by_sh_id) \n"
	"FROM\n"
	"	sw_can_by_sh\n"
	"WHERE\n"
	"	(struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> vec_sw_can_by_sh_id;
	while(res.next())
	{
		Size sw_can_by_sh_id;
		res >> sw_can_by_sh_id ;
		vec_sw_can_by_sh_id.push_back(sw_can_by_sh_id);
	}
	return vec_sw_can_by_sh_id;
} //get_sw_can_by_sh_id



//get_size_sw_by_components_PK_id
Size
SandwichFeatures::get_size_sw_by_components_PK_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	count(sw_by_components_PK_id) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size size_of_sw_by_components_PK_id;
	while(res.next())
	{
		res >> size_of_sw_by_components_PK_id;
	}
	return size_of_sw_by_components_PK_id;
} //get_size_sw_by_components_PK_id



//get_sheet_antiparallel_info
string
SandwichFeatures::get_sheet_antiparallel_info(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	string select_string =
	"SELECT\n"
	"	distinct sheet_antiparallel \n"
	"FROM\n"
	"	sheet \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sheet_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sheet_id);
	result res(basic::database::safely_read_from_database(select_statement));

	string sheet_is_antiparallel;
	while(res.next())
	{
		res >> sheet_is_antiparallel;
	}
	return sheet_is_antiparallel;
} //get_sheet_antiparallel_info


//get_starting_res_for_connecting_strands
std::pair<Size, Size>
SandwichFeatures::get_starting_res_for_connecting_strands(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Size former_res_end)
{
	string select_string =
	"SELECT\n"
	"	min(residue_end) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (residue_end > ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	select_statement.bind(3,former_res_end);
	result res(basic::database::safely_read_from_database(select_statement));

	Size starting_res_for_connecting_strands;
	while(res.next()){
		res >> starting_res_for_connecting_strands;
	}

	select_string =
	"SELECT\n"
	"	sheet_id \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (residue_end = ?);";

	statement select_sh_id_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_sh_id_statement.bind(1,struct_id);
	select_sh_id_statement.bind(2,sw_can_by_sh_id);
	select_sh_id_statement.bind(3,starting_res_for_connecting_strands);
	result res_sh_id(basic::database::safely_read_from_database(select_sh_id_statement));

	Size sheet_id;
	while(res_sh_id.next()){
		res_sh_id >> sheet_id;
	}
	return std::make_pair(starting_res_for_connecting_strands, sheet_id);
} //get_starting_res_for_connecting_strands


//get_next_starting_res_for_connecting_strands
std::pair<Size, Size> //Size
SandwichFeatures::get_next_starting_res_for_connecting_strands(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Size former_ending_res)
{
	string select_string =
	"SELECT\n"
	"	min(residue_begin) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (residue_end > ?);";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	select_statement.bind(3,former_ending_res);
	result res(basic::database::safely_read_from_database(select_statement));

	Size next_starting_res_for_connecting_strands;
	while(res.next()){
		res >> next_starting_res_for_connecting_strands;
	}

	select_string =
	"SELECT\n"
	"	sheet_id \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (residue_begin = ?);";
	statement select_sh_id_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_sh_id_statement.bind(1,struct_id);
	select_sh_id_statement.bind(2,sw_can_by_sh_id);
	select_sh_id_statement.bind(3,next_starting_res_for_connecting_strands);
	result res_sh_id(basic::database::safely_read_from_database(select_sh_id_statement));

	Size sh_id_of_next_start_res;
	while(res_sh_id.next()){
		res_sh_id >> sh_id_of_next_start_res;
	}
	return std::make_pair(next_starting_res_for_connecting_strands, sh_id_of_next_start_res);
} //get_next_starting_res_for_connecting_strands



//update intra/inter_sheet_connection part
Size
SandwichFeatures::update_sheet_con(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_by_components_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	bool intra_sheet_con, // if false, then inter_sheet_con
	Size intra_sheet_con_id,
	Size inter_sheet_con_id,
	string LR,
	string cano_LR,
	string PA_by_preceding_E,
	string PA_by_following_E,
	string cano_PA,
	string heading_direction,
	string parallel_EE,
	string cano_parallel_EE,
	Size loop_size,
	Size start_res,
	Size end_res)
{
	if (intra_sheet_con)
	{
		string insert =
		"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, LR, cano_LR, PA_by_preceding_E, PA_by_following_E,	cano_PA,	heading_direction, parallel_EE, cano_parallel_EE,	loop_size,	residue_begin, residue_end, intra_sheet_con_id)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
		statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
		insert_stmt.bind(1,	struct_id);
		insert_stmt.bind(2,	sw_by_components_PK_id_counter);
		insert_stmt.bind(3,	tag);
		insert_stmt.bind(4,	sw_can_by_sh_id);
		insert_stmt.bind(5,	LR);
		insert_stmt.bind(6,	cano_LR);
		insert_stmt.bind(7,	PA_by_preceding_E);
		insert_stmt.bind(8,	PA_by_following_E);
		insert_stmt.bind(9,	cano_PA);
		insert_stmt.bind(10,	heading_direction);
		insert_stmt.bind(11,	parallel_EE);
		insert_stmt.bind(12,	cano_parallel_EE);
		insert_stmt.bind(13,	loop_size);
		insert_stmt.bind(14,	start_res);
		insert_stmt.bind(15,	end_res);
		insert_stmt.bind(16,	intra_sheet_con_id);
		basic::database::safely_write_to_database(insert_stmt);
	}
	else
	{
		string insert =
		"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, LR, cano_LR, PA_by_preceding_E, PA_by_following_E,	cano_PA,	heading_direction, parallel_EE, cano_parallel_EE,	loop_size,	residue_begin, residue_end, inter_sheet_con_id)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
		statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
		insert_stmt.bind(1,	struct_id);
		insert_stmt.bind(2,	sw_by_components_PK_id_counter);
		insert_stmt.bind(3,	tag);
		insert_stmt.bind(4,	sw_can_by_sh_id);
		insert_stmt.bind(5,	LR);
		insert_stmt.bind(6,	cano_LR);
		insert_stmt.bind(7,	PA_by_preceding_E);
		insert_stmt.bind(8,	PA_by_following_E);
		insert_stmt.bind(9,	cano_PA);
		insert_stmt.bind(10,	heading_direction);
		insert_stmt.bind(11,	parallel_EE);
		insert_stmt.bind(12,	cano_parallel_EE);
		insert_stmt.bind(13,	loop_size);
		insert_stmt.bind(14,	start_res);
		insert_stmt.bind(15,	end_res);
		insert_stmt.bind(16,	inter_sheet_con_id);
		basic::database::safely_write_to_database(insert_stmt);
	}
	return 0;
} //update_intra_sheet_con



//see_whether_other_strands_are_contained
bool
SandwichFeatures::see_whether_other_strands_are_contained(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Size res_num)
{
	string select_string =
	"SELECT\n"
	"	bs.segment_id \n"
	"FROM\n"
	"	secondary_structure_segments AS bs \n"
	"WHERE\n"
	"	(bs.struct_id = ?) \n"
	"	AND (bs.dssp = 'E') \n"
	"	AND (? between bs.residue_begin and bs.residue_end);";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,res_num);
	result res(basic::database::safely_read_from_database(select_statement));

	bool check_bool = false;
	Size segment_id;
	while(res.next())
	{
		res >> segment_id; //segment_id is used below
		check_bool = true;
	}

	if (!check_bool)
	{
		return false; // so no other strands are contained between start_res and next_start_res
	}

	select_string =
	"SELECT\n"
	"	sw_com.sw_by_components_PK_id AS sw_com_PK \n"
	"FROM\n"
	"	sw_by_components AS	sw_com\n"
	"WHERE\n"
	"	(sw_com.struct_id = ?) \n"
	"	AND (sw_com.sw_can_by_sh_id = ?) \n"
	"	AND (sw_com.sw_by_components_bs_id = ?);";
	
	statement select_PK_id_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_PK_id_statement.bind(1,struct_id);
	select_PK_id_statement.bind(2,sw_can_by_sh_id);
	select_PK_id_statement.bind(3,segment_id);
	result res_PK_id(basic::database::safely_read_from_database(select_PK_id_statement));

	check_bool = false;
	while(res_PK_id.next())
	{
		check_bool = true;
	}

	if (!check_bool)
	{
		return true; // so other strands are contained between start_res and next_start_res
	}
	return false;
} //see_whether_other_strands_are_contained



//delete_this_sw_can_by_sh_id
Size
SandwichFeatures::delete_this_sw_can_by_sh_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"DELETE	\n"
	"FROM\n"
	"	sw_by_components	\n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	basic::database::safely_write_to_database(select_statement);
	return 0;
} //delete_this_sw_can_by_sh_id


//get_segment_id
Size
SandwichFeatures::get_segment_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size all_strands_index)
{
	string select_string =
	"SELECT	\n"
	"	segment_id \n"
	"FROM\n"
	"	secondary_structure_segments\n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"   AND (dssp = 'E') \n"
	"	limit ?;";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,all_strands_index);
	result res(basic::database::safely_read_from_database(select_statement));
	Size segment_id;
	while(res.next())
	{
		res >> segment_id;
	}
	return segment_id;
} //get_segment_id

	
Size
SandwichFeatures::get_num_of_distinct_sheet_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	count(distinct sheet_id)\n"
	"FROM\n"
	"	sheet \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"  AND sheet_id != 99999 ;";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));
	
	Size num_distinct_sheet_id;
	while(res.next())
	{
		res >> num_distinct_sheet_id;
	}
	return num_distinct_sheet_id;
} //get_num_of_distinct_sheet_id
	
	
bool
SandwichFeatures::check_helix_existence(
	Pose const & dssp_pose)
{
	for(core::Size ii=1; ii<=dssp_pose.total_residue(); ii++ )
	{
		char res_ss( dssp_pose.secstruct( ii ) ) ;
		if( res_ss == 'H')
		{
			return true;
		}
	}
	return false;
}


string
SandwichFeatures::check_canonicalness_of_LR(
	Size loop_size,
	bool intra_sheet,
	string LR)
{
	// T, -> true, canonical chiral
	// F, -> false, non-canonical chiral
	// U, -> uncertain, this loop-size with this condition has no definite canonical chiral reference in the first place!

	//check_canonicalness_of_LR is same whether canocheck_canonicalness_cutoff_ is 80% or 75%

	if (loop_size == 2)
	{
		if (LR=="L" || LR=="BL" )	{return "T";}
		else	{return "F";}
	}
	if (loop_size == 3)
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "T";}
			else	{return "F";}
		}
		else {return "U";}
	}
	if (loop_size == 4)
	{
		if (intra_sheet)	{return "U";}
		else
		{
			if (LR=="L" || LR=="BL" )	{return "T";}
			else	{return "F";}
		}
	}
	if (loop_size == 5)
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "F";}
			else	{return "T";}
		}
		else {return "U";}
	}
	if (loop_size == 6)
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "F";}
			else	{return "T";}
		}
		else	{return "U";}
	}
	if (loop_size == 11)
	{
		if (intra_sheet)	{return "U";}
		else
		{
			if (LR=="L" || LR=="BL" )	{return "T";}
			else	{return "F";}
		}
	}
	else // all else loop sizes
	{
		return "U";
	}
} //check_canonicalness_of_LR

string
SandwichFeatures::check_canonicalness_of_PA(
	Size loop_size,
	bool intra_sheet,
	string PA_by_preceding_E,
	string PA_by_following_E,
	Real canocheck_canonicalness_cutoff_)
{
	// T, -> true, canonical PA
	// F, -> false, non-canonical PA
	// U, -> uncertain, this loop-size with this condition has no definite canonical PA reference in the first place!

	if (canocheck_canonicalness_cutoff_ == 80)
	{
		if (loop_size == 3)
		{
			if (intra_sheet)
			{
				if (PA_by_following_E=="P")	{return "T";}
				else	{return "F";}
			}
			else
			{
				return "U";
			}
		}
		if (loop_size == 5)
		{
			if (intra_sheet)
			{
				if (PA_by_preceding_E=="A" && PA_by_following_E=="P")	{return "T";}
				else	{return "F";}
			}
			else {return "U";}
		}
		if (loop_size == 6)
		{
			if (intra_sheet)
			{
				if (PA_by_following_E=="P") {return "T";}
				else	{return "F";}
			}
			else	{return "U";}
		}
		if (loop_size == 10)
		{
			if (intra_sheet) {	return "U";	}
			else
			{
				if (PA_by_following_E=="P") {return "T";}
				else	{return "F";}
			}
		}
		if (loop_size == 11)
		{
			if (intra_sheet) {	return "U";	}
			else
			{
				if (PA_by_preceding_E=="A") {return "T";}
				else	{return "F";}
			}
		}
		else // all else loop sizes
		{
			return "U";
		}
	} //if (canocheck_canonicalness_cutoff_ == 80)

	else // (canocheck_canonicalness_cutoff_ == 75)
	{
		if (loop_size == 2)
		{
			if (intra_sheet)
			{
				if (PA_by_preceding_E=="A")	{return "T";}
				else	{return "F";}
			}
			else {return "U";}
		}
		if (loop_size == 3)
		{
			if (intra_sheet)
			{
				if (PA_by_following_E=="P")	{return "T";}
				else	{return "F";}
			}
			else
			{
				if (PA_by_following_E=="A")	{return "T";}
				else	{return "F";}
			}
		}
		if (loop_size == 5)
		{
			if (intra_sheet)
			{
				if (PA_by_preceding_E=="A" && PA_by_following_E=="P")	{return "T";}
				else	{return "F";}
			}
			else {return "U";}
		}
		if (loop_size == 6)
		{
			if (intra_sheet)
			{
				if (PA_by_following_E=="P") {return "T";}
				else	{return "F";}
			}
			else	{return "U";}
		}
		if (loop_size == 8)
		{
			if (intra_sheet)	{return "U";}
			else
			{
				if (PA_by_following_E=="A") {return "T";}
				else	{return "F";}
			}
		}
		if (loop_size == 9)
		{
			if (intra_sheet)
			{
				if (PA_by_preceding_E=="A") {return "T";}
				else	{return "F";}
			}
			else	{	return "U";	}
		}
		if (loop_size == 10)
		{
			if (intra_sheet) {	return "U";	}
			else
			{
				if (PA_by_following_E=="P") {return "T";}
				else	{return "F";}
			}
		}
		if (loop_size == 11)
		{
			if (intra_sheet) {	return "U";	}
			else
			{
				if (PA_by_preceding_E=="A") {return "T";}
				else	{return "F";}
			}
		}
		else // all else loop sizes
		{
			return "U";
		}
	} //if (canocheck_canonicalness_cutoff_ == 75)
} // check_canonicalness_of_PA


string
SandwichFeatures::check_canonicalness_of_parallel_EE(
	Size loop_size,
	bool intra_sheet,
	string parallel_EE)
{
	// T, -> true, canonical parallel_EE
	// F, -> false, non-canonical parallel_EE
	// U, -> uncertain, this loop-size with this condition has no definite canonical parallel_EE reference in the first place!

	if (loop_size == 2 || loop_size == 4)
	{
		if (intra_sheet)
		{
			if (parallel_EE=="P_EE")	{return "T";}
			else	{return "F";}
		}
		else
		{	return "U";	}
	}
	if (loop_size == 3)
	{
		if (intra_sheet)			{	return "U";	}
		else
		{
			if (parallel_EE=="A_EE")	{return "T";}
			else	{return "F";}
		}
	}
	if (loop_size == 6)
	{
		if (intra_sheet)
		{
			if (parallel_EE=="A_EE")	{return "T";}
			else	{return "F";}
		}
		else
		{	return "U";	}
	}
	if (loop_size == 11)
	{
		if (intra_sheet)			{	return "U";	}
		else
		{
			if (parallel_EE=="P_EE")	{return "T";}
			else	{return "F";}
		}
	}
	else // all else loop sizes
	{
		return "U";
	}
} //check_canonicalness_of_parallel_EE


//check_whether_sheets_are_connected_by_same_direction_strand
bool
SandwichFeatures::check_whether_sheets_are_connected_by_same_direction_strand(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Pose const & pose,
	Size start_res,
	Size next_start_res)
{
	//get other terminus of next_start_res
	string select_string =
	"SELECT\n"
	"	residue_end	\n"
	"FROM\n"
	"	secondary_structure_segments \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND dssp = 'E' \n"
	"	AND residue_begin = ?;";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,next_start_res);
	result res(basic::database::safely_read_from_database(select_statement));
	
	Size other_end_of_next_start_res;
	while(res.next())
	{
		res >> other_end_of_next_start_res;
	}

	// angle of start_res, next_start_res and other terminus of next_start_res
	Vector const& first_0_xyz    ( pose.residue(start_res).xyz("CA") );
	Vector const& middle_0_xyz   ( pose.residue(next_start_res).xyz("CA") );
	Vector const& third_0_xyz    ( pose.residue(other_end_of_next_start_res).xyz("CA") );

	Real angle_next_start_res_being_middle = numeric::angle_degrees(first_0_xyz, middle_0_xyz, third_0_xyz);
	if (angle_next_start_res_being_middle > 130)
	{
		return true; // don't use this sandwich, this sandwich is not our target
	}
	
	//get other terminus of start_res
	select_string =
	"SELECT\n"
	"	residue_begin	\n"
	"FROM\n"
	"	secondary_structure_segments \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND dssp = 'E' \n"
	"	AND residue_end = ?;";
	
	statement select_statement_start_res(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement_start_res.bind(1,struct_id);
	select_statement_start_res.bind(2,start_res);
	result res_start_res(basic::database::safely_read_from_database(select_statement_start_res));
	
	Size other_end_of_start_res;
	while(res_start_res.next())
	{
		res_start_res >> other_end_of_start_res;
	}

	// angle of other terminus of start_res, start_res, and next_start_res
	Vector const& first_1_xyz    ( pose.residue(other_end_of_start_res).xyz("CA") );
	Vector const& middle_1_xyz   ( pose.residue(start_res).xyz("CA") );
	Vector const& third_1_xyz    ( pose.residue(next_start_res).xyz("CA") );

	Real angle_start_res_being_middle = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
	if (angle_start_res_being_middle > 130)
	{
		return true; // don't use this sandwich, this sandwich is not our target
	}

	return false; // use this sandwich
} //check_whether_sheets_are_connected_by_same_direction_strand


///@brief collect all the feature data for the pose
core::Size
SandwichFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session)
{
		//	TR.Info << "======================= <begin> report_features =========================" << endl;
	pose::Pose dssp_pose ( pose ); //copy of pose, since the original pose is called as 'const'
	core::scoring::dssp::Dssp dssp( dssp_pose );
	dssp.insert_ss_into_pose( dssp_pose );

	if (no_helix_in_pdb_)
	{
		bool helix_existence = check_helix_existence(dssp_pose);
		if (helix_existence)
		{
				TR.Info << "Exit early since this pdb has helix " << endl;
			return 0;
		}
	}

	Size sheet_PK_id_counter=1; //initial value
	Size sw_can_by_sh_PK_id_counter=1; //initial value
	Size sw_can_by_sh_id_counter=1; //initial value
	Size sw_by_components_PK_id_counter=1; //initial value
	Size intra_sheet_con_id_counter=1; //initial value
	Size inter_sheet_con_id_counter=1; //initial value

	utility::vector1<SandwichFragment> all_strands = get_full_strands(struct_id, db_session);
		//TR.Info << "all_strands.size(): " << all_strands.size() << endl;
	
	if (all_strands.size() < min_num_strands_to_deal_ || all_strands.size() > max_num_strands_to_deal_ )
	{
			TR.Info << "Exit early since all_strands.size(): " << all_strands.size() << endl;
		return 0;
	}

// define very first sheet ("1")
// <begin> assignment of strands into the sheet
	bool first_sheet_assigned = false;
	for(Size i=1; i<all_strands.size() && !first_sheet_assigned; ++i) // I don't need the last strand since this double for loops are exhaustive search for all pairs of strands
	{
		if (all_strands[i].get_size() < min_strand_size_)
		{
			continue;
		}
		for(Size j=i+1; j<=all_strands.size() && !first_sheet_assigned; ++j) // I need the last strand for this second for loop
		{
			if (all_strands[j].get_size() < min_strand_size_)
			{
				continue;
			}
			SandwichFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
			SandwichFragment temp_strand_j(all_strands[j].get_start(), all_strands[j].get_end());
				//TR.Info << "--- anti-parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand begins --- " << endl;

			Size return_of_find_sheet_antiparallel(0); // temporary 'false' designation
			Size return_of_find_sheet_parallel(0); // temporary 'false' designation
			return_of_find_sheet_antiparallel = find_sheet (pose, temp_strand_i, temp_strand_j, true);
			if (return_of_find_sheet_antiparallel == 999)
			{
				break;
			}
			if (!return_of_find_sheet_antiparallel)
			{
				return_of_find_sheet_parallel = find_sheet (pose, temp_strand_i, temp_strand_j, false);
			}
			if (return_of_find_sheet_parallel == 999)
			{
				break;
			}
			if (return_of_find_sheet_antiparallel || return_of_find_sheet_parallel)
			{
				first_sheet_assigned = true;
				write_to_sheet (
					struct_id,
					db_session,
					sheet_PK_id_counter,
					1, // sheet_id
					get_segment_id( // segment_id
						struct_id,
						db_session,
						i)); 

				sheet_PK_id_counter++;
				write_to_sheet (
					struct_id,
					db_session,
					sheet_PK_id_counter,
					1, // sheet_id
					get_segment_id( // segment_id
						struct_id,
						db_session,
						j));
				sheet_PK_id_counter++;
			}
		}
	}
// <end> assignment of strand into sheet



// <begin> assignment of strand into rest sheets (other than "1")
	for(Size i=1; i<=all_strands.size(); ++i)
	{
		if (all_strands[i].get_size() >= min_strand_size_) // the length of this beta strand > min_strand_size_
		{
			bool strand_i_is_in_any_sheet = check_whether_strand_i_is_in_sheet(struct_id, db_session, get_segment_id(
				struct_id,
				db_session,
				i));
			if (!strand_i_is_in_any_sheet) //  this strand is not in any sheet
			{
					//TR.Info << "then this strand needs to be assigned into any sheet" << endl; ;
				utility::vector1<SandwichFragment> cur_strands = get_current_strands_in_sheet(struct_id, db_session);
				bool sheet_unassigned = true;
				for(Size j=1; j<=cur_strands.size() && sheet_unassigned == true; ++j)
				{
					SandwichFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
					SandwichFragment temp_strand_j(cur_strands[j].get_start(), cur_strands[j].get_end());

					Size return_of_find_sheet_antiparallel(0); // temporary 'false' designation
					Size return_of_find_sheet_parallel(0); // temporary 'false' designation

					return_of_find_sheet_antiparallel = find_sheet (pose, temp_strand_i, temp_strand_j, true);

					if (return_of_find_sheet_antiparallel == 999)
					{
						break; // too distant strands
					}

					if (!return_of_find_sheet_antiparallel)
					{
						return_of_find_sheet_parallel = find_sheet (pose, temp_strand_i, temp_strand_j, false);
					}

					if (return_of_find_sheet_parallel == 999)
					{
						break; // too distant strands
					}

					if (return_of_find_sheet_antiparallel || return_of_find_sheet_parallel)
					{
						sheet_unassigned = false;
						write_to_sheet (struct_id, db_session, sheet_PK_id_counter, cur_strands[j].get_sheet_id(), get_segment_id(
																																  struct_id,
																																  db_session,
																																  i));	// struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
						sheet_PK_id_counter++;
					}
				} //for(Size j=1; j<=cur_strands.size(); ++j)

				if (sheet_unassigned)
				{
					Size max_sheet_id = get_max_sheet_id(struct_id, db_session);
					write_to_sheet (struct_id, db_session, sheet_PK_id_counter, max_sheet_id+1, get_segment_id(
																											   struct_id,
																											   db_session,
																											   i));	// struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
					sheet_PK_id_counter++;
				}
			}
		} //all_strands[i].get_size() >= min_strand_size_

		else // all_strands[i].get_size() < min_strand_size_ 				TR.Info << "this strand is too small, assign it into '99999' sheet" << endl;
		{
			write_to_sheet (struct_id, db_session, sheet_PK_id_counter, 99999, get_segment_id(
																							  struct_id,
																							  db_session,
																							  i));	// struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
			sheet_PK_id_counter++;
		}// all_strands[i].get_size() < min_strand_size_
	}
// <end> assignment of strand into rest sheets (other than "1")


// <begin> see_whether_sheet_is_antiparallel
	utility::vector1<Size> all_distinct_sheet_ids = get_distinct_sheet_id(struct_id, db_session);
	for(Size i=1; i<all_distinct_sheet_ids.size()+1; ++i)
	{
		if (all_distinct_sheet_ids[i] == 99999) //all_strands[i].get_size() < min_strand_size_
		{
			continue;
		}
		string sheet_is_antiparallel = see_whether_sheet_is_antiparallel(
										struct_id,
										db_session,
										pose,
										all_distinct_sheet_ids[i]); //sheet id
		update_sheet_antiparallel(struct_id, db_session, all_distinct_sheet_ids[i], sheet_is_antiparallel);
	}
// <end> see_whether_sheet_is_antiparallel


// <begin> redefine sheet id
	bool sheet_id_changed =	true; //temp bool
	while (sheet_id_changed)
	{
		sheet_id_changed =	change_sheet_id_if_possible(
								struct_id,
								db_session,
								pose);
	}
// <end> redefine sheet id
	
	if (!extract_sandwich_)
	{
		return 0;
	}
	string tag = get_tag(struct_id, db_session);
	
/////////////////// <begin> assignment of sheet into sandwich_candidate_by_sheet (sw_can_by_sh)
		TR.Info << "<begin> assignment of sheet into sandwich_candidate_by_sheet (sw_can_by_sh): " << endl;
	Size sheet_j_that_will_be_used_for_pairing_with_sheet_i = 99999; // temp_value
	Size j;
	for(Size i=1; i <= all_distinct_sheet_ids.size()-1; ++i) // now I check possible combination only so i < all_distinct_sheet_ids.size()/2
	{
		bool compared = false;
		if (all_distinct_sheet_ids[i] == sheet_j_that_will_be_used_for_pairing_with_sheet_i) // useful to exclude sheet later
		{
			continue;
		}
		if (all_distinct_sheet_ids[i] == 99999) //	all_strands[i].get_size() < min_strand_size_
		{
			continue;
		}
		Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id
			
		if (num_strands_i < min_num_strands_in_sheet_)
		{
			continue;
		}

		Real lowest_avg_dis_between_sheets = 9999; //temp value
		Real avg_dis_between_sheets;

		// <begin> identify sheet_j_that_will_be_used_for_pairing_with_sheet_i to be a sandwich
		for(j=i+1; j<=all_distinct_sheet_ids.size(); ++j)
		{
			if (all_distinct_sheet_ids[j] == sheet_j_that_will_be_used_for_pairing_with_sheet_i)	// useful to exclude sheet later
			{
				continue;
			}
			if (all_distinct_sheet_ids[j] == 99999)//	all_strands[i].get_size() < min_strand_size_
			{
				continue;
			}
			Size num_strands_j = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]); // struct_id, db_session, sheet_id
			
			if (num_strands_j < min_num_strands_in_sheet_ )
			{
				continue;
			}
			compared = true;

			// now candidates of sheet pair to be sw are identified
				TR.Info << "sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
				TR.Info << "sheet_id (all_distinct_sheet_ids[j]): " << all_distinct_sheet_ids[j] << endl;

			utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
			utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[j]);

			bool these_2_sheets_are_too_close = false; // temporary 'false' designation

			// <begin> check_strand_too_closeness
			for(Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_close; ++ii)
			{
				for(Size jj=1; jj<=strands_from_sheet_j.size() && !these_2_sheets_are_too_close; ++jj)
				{
					SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
					SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
					
					bool are_strands_too_close = check_strand_too_closeness (pose, temp_strand_i, temp_strand_j);
					if (are_strands_too_close)
					{
							TR.Info << "these two sheets are too close by its strands" << endl;
						these_2_sheets_are_too_close = true;
					}
				}
			}
			// <end> check_strand_too_closeness

			if (these_2_sheets_are_too_close)
			{
				break; // break j sheet loop
			}

			// <begin> check whether strands are too distant to each other
			bool these_2_sheets_are_too_distant = false; // temporary 'false' designation
			avg_dis_between_sheets = 0;
			for(Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_close; ++ii)
			{
				vector<Real> array_avg_dis_sheet_i_j;
				for(Size jj=1; jj<=strands_from_sheet_j.size() && !these_2_sheets_are_too_close; ++jj)
				{
					SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
					SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
					
					Real avg_dis_strands = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
					array_avg_dis_sheet_i_j.push_back(avg_dis_strands);
				}

				Real sum = 0;
				for(Size k=0; k<array_avg_dis_sheet_i_j.size(); ++k)
				{
					sum += array_avg_dis_sheet_i_j[k];
				}
				avg_dis_between_sheets = sum / array_avg_dis_sheet_i_j.size();

				Real shortest_avg_dis_inter_sheet = 9999;
				
				for(Size kk=1; kk<=strands_from_sheet_j.size() && !these_2_sheets_are_too_close; ++kk)
				{
					if (shortest_avg_dis_inter_sheet > array_avg_dis_sheet_i_j[kk-1])
					{
						shortest_avg_dis_inter_sheet = array_avg_dis_sheet_i_j[kk-1];
					}
				}
					TR.Info << "shortest_avg_dis_inter_sheet: " << shortest_avg_dis_inter_sheet << endl;
				if (shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_)
				{
						TR.Info << "shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_" << endl;
					these_2_sheets_are_too_distant = true;
				}
			} // for(Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_close; ++ii)
			// <end> check whether strands are too distant to each other

				TR.Info << "avg_dis_between_sheets: " << avg_dis_between_sheets << endl;
			if (these_2_sheets_are_too_distant)
			{
				break; // break j sheet loop
			}
			if (!these_2_sheets_are_too_close && (lowest_avg_dis_between_sheets > avg_dis_between_sheets))
			{
				lowest_avg_dis_between_sheets = avg_dis_between_sheets;
				// use all_distinct_sheet_ids[j] to pair with all_distinct_sheet_ids[i]
				sheet_j_that_will_be_used_for_pairing_with_sheet_i = all_distinct_sheet_ids[j];
			}
		} // for(Size j=i+1; j<=all_distinct_sheet_ids.size(); ++j)

		// <end> identify sheet_j_that_will_be_used_for_pairing_with_sheet_i to be a sandwich

			TR.Info << "sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
			TR.Info << "sheet_id (all_distinct_sheet_ids[j-1]): " << all_distinct_sheet_ids[j-1] << endl;

			TR.Info << "sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
			TR.Info << "sheet_id (sheet_j_that_will_be_used_for_pairing_with_sheet_i): " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << endl;

		// <begin> check_sw_by_distance
		bool found_sandwich_w_these_2_sheets = false; // temporary 'false' designation
		bool chance_of_being_sandwich_w_these_2_sheets = true; // temporary 'true' designation
		utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
		utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, sheet_j_that_will_be_used_for_pairing_with_sheet_i);

			TR.Info << "<begin> check_sw_by_distance" << endl;
		while (!found_sandwich_w_these_2_sheets && chance_of_being_sandwich_w_these_2_sheets)
		{
			for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich_w_these_2_sheets; ++ii)
			{
				for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich_w_these_2_sheets; ++jj)
				{
					SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
					SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
					
					Real return_of_check_sw_by_dis_anti = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, true);
					Real return_of_check_sw_by_dis_parallel = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, false);
					
					if ( return_of_check_sw_by_dis_anti == -999 || return_of_check_sw_by_dis_parallel == -999)
					{
							TR.Info << "these sheets will not be sandwich ever because these are too close or distant to each other!" << endl;
						chance_of_being_sandwich_w_these_2_sheets = false;
					}
					
					if ( return_of_check_sw_by_dis_anti != -99 || return_of_check_sw_by_dis_parallel != -99)
					{
							TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " are in the ideal distance range" << endl;
						found_sandwich_w_these_2_sheets = true;
						chance_of_being_sandwich_w_these_2_sheets = false; // these are sandwich, but no more sheet search is needed! (this "false" false assignment is needed (confirmed! by experiment))
					}
				} //for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich_w_these_2_sheets; ++jj)
			} //for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich_w_these_2_sheets; ++ii)
			break; // no sandwich here
		} //while (!found_sandwich_w_these_2_sheets && chance_of_being_sandwich_w_these_2_sheets)
			TR.Info << "<end> check_sw_by_distance" << endl;
		// <end> check_sw_by_distance
		
		if (!found_sandwich_w_these_2_sheets)
		{
				TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " are not in the ideal distance range so these cannot constitute a sandwich" << endl;
			continue;
		}
		
		bool facing = judge_facing(struct_id, db_session, pose, all_distinct_sheet_ids[i], sheet_j_that_will_be_used_for_pairing_with_sheet_i);
		// if false, these two strand_pairs are linear to each other or do not face properly to each other
		
		if (!facing)
		{
				TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " do not face each other" << endl;
			continue;
		}
		
			TR.Info << "! writing into 'sandwich candidate by sheet' !" << endl;
		
		write_to_sw_can_by_sh (struct_id, db_session, sw_can_by_sh_PK_id_counter, tag, sw_can_by_sh_id_counter, all_distinct_sheet_ids[i], strands_from_sheet_i.size());
		sw_can_by_sh_PK_id_counter++;
		
		write_to_sw_can_by_sh (struct_id, db_session, sw_can_by_sh_PK_id_counter, tag, sw_can_by_sh_id_counter, sheet_j_that_will_be_used_for_pairing_with_sheet_i, strands_from_sheet_j.size());
		sw_can_by_sh_PK_id_counter++;
		
		sw_can_by_sh_id_counter++;
			
	} // for(Size i=1; i<all_distinct_sheet_ids.size(); ++i)
/////////////////// <end> assignment of sheet into sw_can_by_sh

			
						



/////////////////// <begin> fill a table 'sw_by_components' by secondary_structure_segments
		TR.Info << "<begin> fill a table 'sw_by_components' by secondary_structure_segments" << endl;

	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh = prepare_to_fill_sw_by_components(struct_id, db_session); 
		// It retrieves beta segments of sandwich_candidate_by_sheets, it does not make sandwich_candidate_by_components

	if (bs_of_sw_can_by_sh.size() == 0)
	{
			TR.Info << "no beta segment in sandwich_by_sheet (maybe these are too distant sheets or a beta barrel) " << endl;
			TR.Info << "<Done> for this pdb including extraction of sandwich" << endl;
		return 0;
	}

	if (write_phi_psi_)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string phi_psi_file_name = pdb_file_name + "_phi_psi_of_strand_res.txt";
		ofstream phi_psi_file;
		phi_psi_file.open(phi_psi_file_name.c_str());	
		phi_psi_file << "tag	res_num	res_AA	res_at_terminal	sheet_is_antiparallel	strand_is_at_edge	phi	psi" << endl;
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii)
		{
			if (sw_by_components_PK_id_counter > max_num_sw_per_pdb_)
			{
				break;
			}
			string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
			Size strand_is_at_edge = see_whether_strand_is_at_edge	(
										pose,
										struct_id,
										db_session,
										bs_of_sw_can_by_sh[ii].get_sheet_id(),
										sheet_antiparallel,
										bs_of_sw_can_by_sh[ii].get_start(),
										bs_of_sw_can_by_sh[ii].get_end());
			
			fill_sw_by_components (struct_id, db_session, sw_by_components_PK_id_counter, tag, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), sheet_antiparallel, bs_of_sw_can_by_sh[ii].get_strand_id(), strand_is_at_edge, bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			sw_by_components_PK_id_counter++;
			
			Size res_at_terminal;		
			for (Size res_num = bs_of_sw_can_by_sh[ii].get_start(); res_num <= bs_of_sw_can_by_sh[ii].get_end(); res_num++)
			{
				if (res_num == bs_of_sw_can_by_sh[ii].get_start() || res_num == bs_of_sw_can_by_sh[ii].get_end())
				{
					res_at_terminal = 1;	
				}
				else
				{
					res_at_terminal = 0;	
				}
				Real phi = pose.phi(res_num);
				phi_psi_file << tag << "	" << res_num << "	" << pose.residue_type(res_num).name3() << "	" << res_at_terminal << "	" <<	sheet_antiparallel << "	" << strand_is_at_edge << "	" << phi << "	";
				
				Real psi = pose.psi(res_num);
				phi_psi_file << psi << endl;
			}
		}
		phi_psi_file.close();
	} //write_phi_psi_
	
	else //!write_phi_psi_
	{
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii)
		{
			if (sw_by_components_PK_id_counter > max_num_sw_per_pdb_)
			{
				break;
			}
			string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
			Size strand_is_at_edge = see_whether_strand_is_at_edge	(
										pose,
										struct_id,
										db_session,
										bs_of_sw_can_by_sh[ii].get_sheet_id(),
										sheet_antiparallel,
										bs_of_sw_can_by_sh[ii].get_start(),
										bs_of_sw_can_by_sh[ii].get_end());
			
			fill_sw_by_components (struct_id, db_session, sw_by_components_PK_id_counter, tag, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), sheet_antiparallel, bs_of_sw_can_by_sh[ii].get_strand_id(), strand_is_at_edge, bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			sw_by_components_PK_id_counter++;			
		}
	}	//!write_phi_psi_
/////////////////// <end> fill a table 'sw_by_components' by secondary_structure_segments
		TR.Info << "<end> fill a table 'sw_by_components' by secondary_structure_segments" << endl;

	
	/////////////////// <begin> update sheet_connecting_loops (2nd judge whether each sandwich_by_sheet_id becomes sandwich_by_components)

		// get_distinct(sw_can_by_sh_id)
	utility::vector1<Size> vec_sw_can_by_sh_id =  get_vec_sw_can_by_sh_id(struct_id, db_session);

	for(Size ii=1; ii<=vec_sw_can_by_sh_id.size(); ++ii)
	{
		bool chance_of_being_canonical_sw	=	true; // not yet decided fate whether this could be canonical sandwich or not, but assumed to be true
		Size size_sw_by_components_PK_id =
		get_size_sw_by_components_PK_id(
			struct_id,
			db_session,
			vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
			);

		bool bool_no_helix_in_loop = true;
		bool bool_no_more_strand_in_loop = true;
		Size former_start_res = 0; //temporary

			// this 'jj' is used for for iteration purpose only and this 'for' loop iterates only for connecting sheets
		for(Size jj=1; (jj<=size_sw_by_components_PK_id-1) && (chance_of_being_canonical_sw) && (bool_no_helix_in_loop) && (bool_no_more_strand_in_loop); ++jj)
		{
			// get_starting_res_for_connecting_strands and its sheet_id
			std::pair<Size, Size>
			start_res_sh_id =
			get_starting_res_for_connecting_strands(
				struct_id,
				db_session,
				vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
				former_start_res);

			Size start_res = start_res_sh_id.first;
			Size sheet_id_of_start_res = start_res_sh_id.second;

			// get_next_starting_res_for_connecting_strands and its sheet_id
			std::pair<Size, Size>
			next_starting_res_for_connecting_strands =
			get_next_starting_res_for_connecting_strands (
				struct_id,
				db_session,
				vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
				start_res); //former_ending_res

			Size next_start_res = next_starting_res_for_connecting_strands.first;
			Size sheet_id_of_next_start_res = next_starting_res_for_connecting_strands.second;

			former_start_res = next_start_res;

			if (no_strand_in_loop_in_extracted_sw_)
			{
				// <begin> check whether there is other strand which is not part of current dealing two sheets
				for(Size kk=start_res+1; kk<=next_start_res-1; ++kk)
				{
					char res_ss( dssp_pose.secstruct( kk ) ) ;
					if (res_ss == 'E')
					{
						bool contains_other_strands = see_whether_other_strands_are_contained(
														struct_id,
														db_session,
														vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
														kk //res_num
														);
						if (contains_other_strands)
						{
							bool_no_more_strand_in_loop = false;
							break;
						}
					}
				}
				// <end> check whether there is other strand which is not part of current dealing two sheets
			}

			if (no_helix_in_extracted_sw_)
			{
				// <begin> check whether there is a helix as a loop in this extracted sandwich candidate
				for(Size kk=start_res+1; kk<=next_start_res-1; ++kk)
				{
					char res_ss( dssp_pose.secstruct( kk ) ) ;
					if (res_ss == 'H')
					{
						bool_no_helix_in_loop = false;
						break;
					}
				}
				// <end> check whether there is a helix as a loop in this extracted sandwich candidate
			}
			
			if (!bool_no_more_strand_in_loop || !bool_no_helix_in_loop)
			{
				delete_this_sw_can_by_sh_id(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
				chance_of_being_canonical_sw = false;
				break; // break jj 'for' loop
			}

			string LR = check_LR(dssp_pose, start_res+1, next_start_res-1);

			std::pair<string, string> PA = check_PA(dssp_pose, start_res+1, next_start_res-1);
			string PA_by_preceding_E = PA.first;
			string PA_by_following_E = PA.second;

			// check whethere it is positive, negative, away or meet
			string heading_direction = check_heading_direction(dssp_pose, start_res+1, next_start_res-1, check_N_to_C_direction_by_);

			if (heading_direction == "except")
			{
				return 0;
			}
			string parallel_EE;
			if (heading_direction == "posi" || heading_direction == "nega")
			{
				parallel_EE = "P_EE";
			}
			else
			{
				parallel_EE = "A_EE"; // to keep chacracter size be same as "parallel"
			}

			Size loop_size = (next_start_res-1) - (start_res+1) + 1;
			if (sheet_id_of_start_res == sheet_id_of_next_start_res)
				// this loop connects sheets as intra-sheet way
			{
				string cano_LR = check_canonicalness_of_LR(loop_size, true, LR); // loop_size, intra_sheet bool, LR
				string cano_PA = check_canonicalness_of_PA(loop_size, true, PA_by_preceding_E, PA_by_following_E, check_canonicalness_cutoff_); // loop_size, intra_sheet bool, 2 PAs
				string cano_parallel_EE = check_canonicalness_of_parallel_EE(loop_size, true, parallel_EE); // loop_size, intra_sheet bool, parallel_EE
				update_sheet_con(
					struct_id,
					db_session,
					sw_by_components_PK_id_counter,
					tag,
					vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
					true,
					intra_sheet_con_id_counter,
					inter_sheet_con_id_counter,
					LR,
					cano_LR,
					PA_by_preceding_E,
					PA_by_following_E,
					cano_PA,
					heading_direction,
					parallel_EE,
					cano_parallel_EE,
					loop_size,
					start_res+1, // new start_res for intra_sheet_con
					next_start_res-1 // new end_res for intra_sheet_con
					);

				sw_by_components_PK_id_counter++;
				intra_sheet_con_id_counter++;
			}
			else // this loop connects sheets as inter-sheet way
			{
				if (exclude_sandwich_that_is_linked_w_same_direction_strand_)
				{
					bool sheets_are_connected_by_same_direction_strand = check_whether_sheets_are_connected_by_same_direction_strand(struct_id,	db_session,	pose,	start_res,	next_start_res);
						TR.Info << "sheets_are_connected_by_same_direction_strand: " << sheets_are_connected_by_same_direction_strand << endl;
					if (sheets_are_connected_by_same_direction_strand)
					{
						delete_this_sw_can_by_sh_id(
							struct_id,
							db_session,
							vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
							);
						chance_of_being_canonical_sw = false;
						break; // break jj 'for' loop
					}
				}
				string cano_LR = check_canonicalness_of_LR(loop_size, false, LR);	// loop_size, intra_sheet bool, LR
				string cano_PA = check_canonicalness_of_PA(loop_size, false, PA_by_preceding_E, PA_by_following_E, check_canonicalness_cutoff_);
														  // loop_size,	intra_sheet bool, PA_ref_1, PA_ref_2, cutoff
				string cano_parallel_EE = check_canonicalness_of_parallel_EE(loop_size, false, parallel_EE);
																			// loop_size, intra_sheet bool, parallel_EE
				update_sheet_con(
					struct_id,
					db_session,
					sw_by_components_PK_id_counter,
					tag,
					vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
					false, //bool intra_sheet_con
					intra_sheet_con_id_counter,
					inter_sheet_con_id_counter,
					LR,
					cano_LR,
					PA_by_preceding_E,
					PA_by_following_E,
					cano_PA,
					heading_direction,
					parallel_EE,
					cano_parallel_EE,
					loop_size,
					start_res+1, // new start_res for inter_sheet_con
					next_start_res-1 // new end_res for inter_sheet_con
					);

				sw_by_components_PK_id_counter++;
				inter_sheet_con_id_counter++;
			}	// this loop connects sheets as inter-sheet way
		} // for(Size jj=1; (jj<=size_sw_by_components_PK_id-1) && (bool_no_helix_in_loop) && (bool_no_more_strand_in_loop); ++jj)
			TR.Info << "chance_of_being_canonical_sw: " << chance_of_being_canonical_sw << endl;
		/*
		/////////// as of 03/13/2013 temporarily suspend the usage of writing chain_B_resnum, since current 'writing chain_B_resnum' is only
		appliable when max_num_sw_per_pdb_ = 1

		// <begin> write chain_B_resNum to a file
		if (write_chain_B_resnum_ && (chance_of_being_canonical_sw))
		{
			Size tag_len = tag.length();
			string pdb_file_name = tag.substr(0, tag_len-5);
			string report_file_name = pdb_file_name + "_chain_B_resNum.txt";
			ofstream report_file;
			report_file.open(report_file_name.c_str());
			utility::vector1<SandwichFragment> chain_B_resNum = get_chain_B_resNum(struct_id, db_session);
			for(Size i=1; i<=chain_B_resNum.size(); ++i)
			{
				report_file << chain_B_resNum[i].get_resNum() << endl;
			}
			report_file.close();
		} //write_chain_B_resnum_
		// <end> write chain_B_resNum to a file
		*/ // chain_B_resnum
	}	// per each sandwich_by_sheet_id
	/////////////////// <end> update sheet_connecting_loops	(2nd judge whether each sandwich_by_sheet_id becomes sandwich_by_components)

		TR.Info << "<Done> for this pdb including extraction of sandwich" << endl;

	return 0;
} //SandwichFeatures::report_features


} //namespace strand_assembly
} //namespace features
} //namespace protocols
