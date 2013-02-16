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

//C library
#include <string>
#include <math.h> // for round

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
min_num_strands_to_deal_(5), // it should be at least 4
max_num_strands_to_deal_(13), // (in 1LD9 chain A) 16 is too many number of strands, takes too long time
min_strand_size_(2), // minimum number of residues in a strand, for edge strand definition & analysis, 4=< is recommended (in 1A8M) min_strand_size = 2, (in 1PMY) min_strand_size = 3
min_CA_CA_dis_(3.5), // (in 1A8M) 'min_CA_CA_dis_= 3.5', (in 1KIT) 'min_CA_CA_dis_= 4.0'
max_CA_CA_dis_(6.2), // (in 1A8M) 'max_CA_CA_dis_= 6.2', (in 1KIT) 'max_CA_CA_dis_= 5.7'
/*min_O_N_dis_(2.4), // 1KIT shows (renumbered residues 178-181 and residues 76-81) show that 2.5 A exist!
max_O_N_dis_(3.1),*/
min_C_O_N_angle_(120), // (in 1L0Q chain A), 138 is the smallest C_O_N_angle (C and O from one sheet, N from other sheet)
min_sheet_dis_(7.0), // 7 Angstrom may seem OK though
max_sheet_dis_(15.0), // 15 Angstrom may seem OK though
min_sheet_angle_(30.0),
max_sheet_angle_(150.0), // (in 1TEN) even 155 degree comes from same sheet!
min_sheet_torsion_cen_res_(-90.0), // with respect to central residues, one torsion angles of 1TEN is 84.9
max_sheet_torsion_cen_res_(90.0),
min_inter_sheet_dis_CA_CA_(4.0), // (in 2WYF) the distance between val and gln is 5.8 A
max_inter_sheet_dis_CA_CA_(20.0), // (in 1ASQ) the average distance between sheet 1 and 4 > 20 A (so these sheets cannot be a sandwich)
extract_sandwich_(true),
write_chain_B_resnum_(true), // if true, write chain_B_resnum file for InterfaceAnalyzer
write_phi_psi_(false), // if true, write phi_psi_file
max_num_sw_per_pdb_(1) // maximum number of sandwiches to be extracted per a pdb file
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
	
	if(option[OptionKeys::strand_assembly::max_num_sw_per_pdb].user()){
		max_num_sw_per_pdb_ = option[OptionKeys::strand_assembly::max_num_sw_per_pdb];
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
	// changed into null 'possible' because of sw_can_by_components

	Column sheet_antiparallel	("sheet_antiparallel",	new DbInteger(), true /* could be null at first, eventually it will not be null though*/, false /*no autoincrement*/); //if 0, it is parallel or combination of antiparallel and parallel

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




/****** <begin> writing sw_can_by_components (sandwich candidate by components such as strands, loops, helices) ******/

	// Columns
	// id of sw_can_by_components

	//unique
	Column sw_can_by_components_PK_id	("sw_can_by_components_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column sw_can_by_components_bs_id	("sw_can_by_components_bs_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	// could be redundant
	Column sw_can_by_components_bs_edge	("sw_can_by_components_bs_edge",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	
	// could be redundant
	Column intra_sheet_con_id	("intra_sheet_con_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	// could be redundant
	Column inter_sheet_con_id	("inter_sheet_con_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column residue_begin("residue_begin", new DbInteger(), false /*not null*/, false /*don't autoincrement*/);
	Column residue_end  ("residue_end", new DbInteger(), false /*not null*/, false /*don't autoincrement*/);



	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sw_can_by_components;
	primary_key_columns_sw_can_by_components.push_back(struct_id);
	primary_key_columns_sw_can_by_components.push_back(sw_can_by_components_PK_id);

	Schema sw_can_by_components("sw_can_by_components",  PrimaryKey(primary_key_columns_sw_can_by_components));

	// add column which is not PrimaryKey nor ForeignKey
	sw_can_by_components.add_column(tag);
	sw_can_by_components.add_column(sw_can_by_sh_id);
	sw_can_by_components.add_column(sheet_id);
	sw_can_by_components.add_column(sheet_antiparallel);
	sw_can_by_components.add_column(sw_can_by_components_bs_id);
	sw_can_by_components.add_column(sw_can_by_components_bs_edge);
	sw_can_by_components.add_column(intra_sheet_con_id);
	sw_can_by_components.add_column(inter_sheet_con_id);

	// ForeignKey
	sw_can_by_components.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
		// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	utility::vector1<Column> residue_begin_fkey_cols;
	residue_begin_fkey_cols.push_back(struct_id);
	residue_begin_fkey_cols.push_back(residue_begin);

	sw_can_by_components.add_foreign_key(ForeignKey(residue_begin_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

	utility::vector1<Column> residue_end_fkey_cols;
	residue_end_fkey_cols.push_back(struct_id);
	residue_end_fkey_cols.push_back(residue_end);

	sw_can_by_components.add_foreign_key(ForeignKey(residue_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

	sw_can_by_components.write(db_session);
/****** <end> writing sw_can_by_components ******/


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


//utility::vector1<SandwichFragment>
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
Size
SandwichFeatures::update_sheet_antiparallel(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sheet_id,
	Size antiparallel)
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
	return 0;
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
	while(res.next()){
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

//get_num_strands
Size
SandwichFeatures::get_num_strands(
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
} //get_num_strands

//prepare_to_fill_sw_can_by_components
utility::vector1<SandwichFragment>
SandwichFeatures::prepare_to_fill_sw_can_by_components(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT sw_sh.sw_can_by_sh_id AS sw_can_by_sh_id, sh.sheet_id AS sheet_id, \n"
	"	bs.segment_id AS sw_can_by_components_bs_id, bs.residue_begin AS	residue_begin, \n"
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
		Size sw_can_by_sh_id, sheet_id, sw_can_by_components_bs_id, residue_begin,	residue_end;
		res >> sw_can_by_sh_id >> sheet_id >> sw_can_by_components_bs_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(sw_can_by_sh_id, sheet_id, sw_can_by_components_bs_id, residue_begin, residue_end));
	}
	return all_strands;

} //prepare_to_fill_sw_can_by_components

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




bool
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
		
		if (vec_cen_resnum_i[ii] == -999) // I ignore a unrepresentative strand
		{
			dis_from_ii = -999;
			vec_dis_sum_from_cen_resnum_i.push_back(dis_from_ii);
			continue;
		}

		for (Size jj=0; jj<=strands_from_i.size()-1; ++jj)
		{				
			if (vec_cen_resnum_i[jj] == -999) // I ignore a unrepresentative strand
			{
				continue;
			}
			Real dis_CA_CA = pose.residue(vec_cen_resnum_i[ii]).atom("CA").xyz().distance(pose.residue(vec_cen_resnum_i[jj]).atom("CA").xyz());
			dis_from_ii = dis_from_ii + dis_CA_CA;
		}
		vec_dis_sum_from_cen_resnum_i.push_back(dis_from_ii);
	}
	// <end> get array of sum of all distances from central residues


	// <begin> get largest distance and its residue index
	Size res_index_having_the_largest_dis = -99;
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
			if (vec_cen_resnum_i[ii] == -999) // unrepresentative_strands
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
			return false;
		}
	}
	return true;
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
	if (dis_CA_CA_0 > 40)
	{
		return -999; // these sheets will not be sandwich ever, since these two sheets are too distant!
	}
	if (dis_CA_CA_0 < min_sheet_dis_ || dis_CA_CA_0 > max_sheet_dis_)
	{
		return -99;
	}

	Real dis_CA_CA_1 = pose.residue(i_resnum_1).atom("CA").xyz().distance(pose.residue(j_resnum_1).atom("CA").xyz());

	if (dis_CA_CA_1 < min_sheet_dis_ || dis_CA_CA_1 > max_sheet_dis_)
	{
		return -99;
	}

	Real dis_CA_CA_2 = pose.residue(i_resnum_2).atom("CA").xyz().distance(pose.residue(j_resnum_2).atom("CA").xyz());

	if (dis_CA_CA_2 < min_sheet_dis_ || dis_CA_CA_2 > max_sheet_dis_)
	{
		return -99;
	}

	Real dis_CA_CA_3 = pose.residue(i_resnum_3).atom("CA").xyz().distance(pose.residue(j_resnum_3).atom("CA").xyz());

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
		//	TR.Info << "let me see distances between strands " << endl;
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

			Real avg_dis_CA_CA = get_avg_dis_CA_CA(pose,	i_resnum,	i_resnum_1, i_resnum_2, i_resnum_3, j_resnum, j_resnum_1, j_resnum_2, j_resnum_3);

			if (avg_dis_CA_CA == -999)
			{
				break; // these sheets will not be sandwich ever, since these two sheets are too distant!
			}

			if (avg_dis_CA_CA == -99)
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
	return floor(x+.5);
} //round


//can this strand represent a terminal strand for inter-sheet angle calculation?
bool
SandwichFeatures::can_this_strand_represent_a_terminal(
	Pose const & pose,
	utility::vector1<SandwichFragment> strands_from_sheet_i,
	Size current_strand_id_as_i)
{

// <begin> identify the closest strand from current_strand_id_as_i
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
// <end> identify the closest strand from current_strand_id_as_i
	
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


//get_two_central_residues
std::pair<Real, Real>
SandwichFeatures::get_two_central_residues(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i)
{
	utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_i);

	// strands_from_sheet_i
	// get central residue numbers
	vector<Real> cen_res_arr_sheet_i;
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i)
	{
		if (strands_from_sheet_i[i].get_size() <= 2)
		{
			continue;
		}
		bool this_strand_represent_a_terminal = can_this_strand_represent_a_terminal (pose, strands_from_sheet_i, i); // pose. strands_from_sheet_i, current_strand_id_as_i

		if (!this_strand_represent_a_terminal)
		{
			continue;
		}

		Real to_be_rounded_i = (strands_from_sheet_i[i].get_start() + strands_from_sheet_i[i].get_end())/(2.0);
		Real cen_resnum_i = round(to_be_rounded_i);
		cen_res_arr_sheet_i.push_back(cen_resnum_i);
	}
	Size array_size = cen_res_arr_sheet_i.size();
	if (array_size == 0)
	{
		return std::make_pair(-99, -99); // this sheet is constituted with 2 residues long strands only
	}

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
	Real terminal_cen_res_pos_1 = cen_res_arr_sheet_i[index_terminal_cen_res_pos_1];
	Real terminal_cen_res_pos_2 = cen_res_arr_sheet_i[index_terminal_cen_res_pos_2];

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
	// <begin> identify four terminal central residues
	std::pair<Real, Real>
	two_central_residues =	get_two_central_residues(
													 struct_id,
													 db_session,
													 pose,
													 sheet_i);
	Real i_ter_cen_1 = two_central_residues.first;
	Real i_ter_cen_2 = two_central_residues.second;

	if (i_ter_cen_1 == -99 && i_ter_cen_2 == -99)
	{
		return false; // I can't choose two central residues since this sheet is constituted with 2 residues long strands only
	}
	two_central_residues =	get_two_central_residues(
													 struct_id,
													 db_session,
													 pose,
													 sheet_j);
	Real j_ter_cen_1 = two_central_residues.first;
	Real j_ter_cen_2 = two_central_residues.second;

	if (j_ter_cen_1 == -99 && j_ter_cen_2 == -99)
	{
		return false; // I can't choose two central residues since this sheet is constituted with 2 residues long strands only
	}
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
	Size sheet_antiparallel,
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
	
	if (sheet_antiparallel)
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
	

Size
SandwichFeatures::fill_sw_can_by_components	(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_components_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size sheet_antiparallel,
	Size sw_can_by_components_bs_id,
	Size strand_is_at_edge,
	Size residue_begin,
	Size residue_end)
{
	string insert =
	"INSERT INTO sw_can_by_components (struct_id, sw_can_by_components_PK_id, tag, sw_can_by_sh_id, sheet_id, sheet_antiparallel, sw_can_by_components_bs_id, sw_can_by_components_bs_edge, residue_begin, residue_end)  VALUES (?,?,?,?,?,?,?,?,?,?);";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	sw_can_by_components_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	sw_can_by_sh_id);
	insert_stmt.bind(5,	sheet_id);
	insert_stmt.bind(6,	sheet_antiparallel);
	insert_stmt.bind(7,	sw_can_by_components_bs_id); //bs_id
	insert_stmt.bind(8,	strand_is_at_edge);
	insert_stmt.bind(9,	residue_begin);
	insert_stmt.bind(10,	residue_end);
	basic::database::safely_write_to_database(insert_stmt);
	return 0;
} //fill_sw_can_by_components


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



//get_size_sw_can_by_components_PK_id
Size
SandwichFeatures::get_size_sw_can_by_components_PK_id(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	count(sw_can_by_components_PK_id) \n"
	"FROM\n"
	"	sw_can_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size size_of_sw_can_by_components_PK_id;
	while(res.next())
	{
		res >> size_of_sw_can_by_components_PK_id;
	}
	return size_of_sw_can_by_components_PK_id;
} //get_size_sw_can_by_components_PK_id



//get_sheet_antiparallel_info
Size
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

	Size sheet_is_antiparallel;
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
	"	sw_can_by_components \n"
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
	"	sw_can_by_components \n"
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
	"	sw_can_by_components \n"
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
	"	sw_can_by_components \n"
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



//update intra_sheet_con
Size
SandwichFeatures::update_sheet_con(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Size sw_can_by_components_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	bool intra_sheet_con, // if false, then inter_sheet_con
	Size intra_sheet_con_id,
	Size inter_sheet_con_id,
	Size start_res,
	Size end_res)
{
	if (intra_sheet_con)
	{
		string insert =
		"INSERT INTO sw_can_by_components (struct_id, sw_can_by_components_PK_id, tag, sw_can_by_sh_id, residue_begin, residue_end, intra_sheet_con_id)  VALUES (?,?,?,?,?,?,?);";
		statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
		insert_stmt.bind(1,	struct_id);
		insert_stmt.bind(2,	sw_can_by_components_PK_id_counter);
		insert_stmt.bind(3,	tag);
		insert_stmt.bind(4,	sw_can_by_sh_id);
		insert_stmt.bind(5,	start_res);
		insert_stmt.bind(6,	end_res);
		insert_stmt.bind(7,	intra_sheet_con_id);
		basic::database::safely_write_to_database(insert_stmt);
	}
	else
	{
		string insert =
		"INSERT INTO sw_can_by_components (struct_id, sw_can_by_components_PK_id, tag, sw_can_by_sh_id, residue_begin, residue_end, inter_sheet_con_id)  VALUES (?,?,?,?,?,?,?);";
		statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
		insert_stmt.bind(1,	struct_id);
		insert_stmt.bind(2,	sw_can_by_components_PK_id_counter);
		insert_stmt.bind(3,	tag);
		insert_stmt.bind(4,	sw_can_by_sh_id);
		insert_stmt.bind(5,	start_res);
		insert_stmt.bind(6,	end_res);
		insert_stmt.bind(7,	inter_sheet_con_id);
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
	"	sw_com.sw_can_by_components_PK_id AS sw_com_PK \n"
	"FROM\n"
	"	sw_can_by_components AS	sw_com\n"
	"WHERE\n"
	"	(sw_com.struct_id = ?) \n"
	"	AND (sw_com.sw_can_by_sh_id = ?) \n"
	"	AND (sw_com.sw_can_by_components_bs_id = ?);";
	
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
	"	sw_can_by_components	\n"
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

	Size sheet_PK_id_counter=1; //initial value
	Size sw_can_by_sh_PK_id_counter=1; //initial value
	Size sw_can_by_sh_id_counter=1; //initial value
	Size sw_can_by_components_PK_id_counter=1; //initial value
	Size intra_sheet_con_id_counter=1; //initial value
	Size inter_sheet_con_id_counter=1; //initial value

	utility::vector1<SandwichFragment> all_strands = get_full_strands(struct_id, db_session);
		TR.Info << "all_strands.size(): " << all_strands.size() << endl;
	
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
				write_to_sheet (struct_id, db_session, sheet_PK_id_counter, 1, get_segment_id(
																							  struct_id,
																							  db_session,
																							  i));	// struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id

				sheet_PK_id_counter++;
				write_to_sheet (struct_id, db_session, sheet_PK_id_counter, 1, get_segment_id(
																							  struct_id,
																							  db_session,
																							  j));	// struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
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
		if (all_distinct_sheet_ids[i] == 99999)
		{
			continue;
		}
		
		bool sheet_is_antiparallel = see_whether_sheet_is_antiparallel(
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
	
//<begin> assignment of sheet into sw_can_by_sh
	
	for(Size i=1; i<all_distinct_sheet_ids.size(); ++i)
	{
		if (sw_can_by_sh_id_counter > max_num_sw_per_pdb_)
		{
			break;
		}
		if (all_distinct_sheet_ids[i] == 99999)
		{
			continue;
		}
		for(Size j=i+1; j<=all_distinct_sheet_ids.size(); ++j)
		{
			if (sw_can_by_sh_id_counter > max_num_sw_per_pdb_)
			{
				break;
			}
			if (all_distinct_sheet_ids[j] == 99999)
			{
				continue;
			}
						
			Size num_strands_i = get_num_strands(struct_id, db_session, all_distinct_sheet_ids[i]);
			Size num_strands_j = get_num_strands(struct_id, db_session, all_distinct_sheet_ids[j]);
			
			if (num_strands_i < 2 || num_strands_j < 2 )
			{
				continue;
			}
			
			utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
			utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[j]);
			
			bool found_sandwich = false; // temporary 'false' designation
			bool chance_of_being_sandwich = true; // temporary 'true' designation
			
			
			// <begin> check_strand_too_closeness
			for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich; ++ii)
			{
				for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich; ++jj)
				{
					SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
					SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
					
					bool are_strands_too_close = check_strand_too_closeness (pose, temp_strand_i, temp_strand_j);
					
					if (are_strands_too_close)
					{
						chance_of_being_sandwich = false;
					}
				}
			}
			// <end> check_strand_too_closeness
			
			
			// <begin> check_strand_too_distantness
			TR.Info << "<begin> check_strand_too_distantness" << endl;
			for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich; ++ii)
			{
				vector<Real> array_avg_dis_sheet_i_j;
				for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich; ++jj)
				{
					SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
					SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
					
					Real avg_dis_strands = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
					array_avg_dis_sheet_i_j.push_back(avg_dis_strands);
				}
				
				Real shortest_avg_dis_inter_sheet = 9999;
				
				for(Size kk=1; kk<=strands_from_sheet_j.size() && chance_of_being_sandwich; ++kk)
				{
					if (shortest_avg_dis_inter_sheet > array_avg_dis_sheet_i_j[kk-1])
					{
						shortest_avg_dis_inter_sheet = array_avg_dis_sheet_i_j[kk-1];
					}
				}
				if (shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_)
				{
					chance_of_being_sandwich = false;
				}
				
			} // for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich; ++ii)
			// <end> check_strand_too_distantness
			
			
			// <begin> check_sw_by_distance
			TR.Info << "<begin> check_sw_by_distance" << endl;
			while (!found_sandwich && chance_of_being_sandwich)
			{
				for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich; ++ii)
				{
					for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich; ++jj)
					{
						SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
						SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
						
						Real return_of_check_sw_by_dis_anti = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, true);
						Real return_of_check_sw_by_dis_parallel = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, false);
						
						if ( return_of_check_sw_by_dis_anti == -999 || return_of_check_sw_by_dis_parallel == -999)
						{
								TR.Info << "these sheets will not be sandwich ever since these are too close or distant!" << endl;
							chance_of_being_sandwich = false;
							break;
						}
						
						if ( return_of_check_sw_by_dis_anti != -99 || return_of_check_sw_by_dis_parallel != -99)
						{
								TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << all_distinct_sheet_ids[j] << " are in the ideal distance range" << endl;
							found_sandwich = true;
							chance_of_being_sandwich = false; // these are sandwich, but no more sheet search is needed!
							break;
						}
					} //for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich; ++jj)
				} //for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich; ++ii)
				break; // no sandwich here
			} //while (!found_sandwich && chance_of_being_sandwich)
			// <end> check_sw_by_distance
			
			if (!found_sandwich)
			{
					TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << all_distinct_sheet_ids[j] << " are not in the ideal distance range so these cannot constitute a sandwich" << endl;
				continue;
			}
			
			bool facing = judge_facing(struct_id, db_session, pose, all_distinct_sheet_ids[i], all_distinct_sheet_ids[j]);
			// if false, these two strand_pairs are linear to each other or do not face properly to each other
			
			if (!facing)
			{
					TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << all_distinct_sheet_ids[j] << " do not face each other" << endl;
				continue;
			}
			
				TR.Info << "! writing into 'sandwich candidate by sheet' !" << endl;
			
			write_to_sw_can_by_sh (struct_id, db_session, sw_can_by_sh_PK_id_counter, tag, sw_can_by_sh_id_counter, all_distinct_sheet_ids[i], strands_from_sheet_i.size());
			sw_can_by_sh_PK_id_counter++;
			
			write_to_sw_can_by_sh (struct_id, db_session, sw_can_by_sh_PK_id_counter, tag, sw_can_by_sh_id_counter, all_distinct_sheet_ids[j], strands_from_sheet_j.size());
			sw_can_by_sh_PK_id_counter++;
			
			sw_can_by_sh_id_counter++;					
		} // for(Size j=i+1; j<=all_distinct_sheet_ids.size(); ++j)
	} // for(Size i=1; i<all_distinct_sheet_ids.size(); ++i)
//<end> assignment of sheet into sw_can_by_sh

	Size num_distinct_sheet_id = get_num_of_distinct_sheet_id(
															struct_id,
															db_session);
	if (num_distinct_sheet_id > 2)
	{
			TR.Info << "exception: current beta-sandwich maybe somehow non-canonical by definition of dssp like 1AC6 " << endl;
		return 0;
	}


// <begin> fill a table 'sw_can_by_components' by secondary_structure_segments
		TR.Info << "<begin> fill a table 'sw_can_by_components' by secondary_structure_segments" << endl;
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh = prepare_to_fill_sw_can_by_components(struct_id, db_session); // beta segments of sandwich candidate by sheets

	if (bs_of_sw_can_by_sh.size() == 0)
	{
			TR.Info << "maybe this is a beta barrel " << endl;
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
			Size sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
			Size strand_is_at_edge = see_whether_strand_is_at_edge	(
																	 pose,
																	 struct_id,
																	 db_session,
																	 bs_of_sw_can_by_sh[ii].get_sheet_id(),
																	 sheet_antiparallel,
																	 bs_of_sw_can_by_sh[ii].get_start(),
																	 bs_of_sw_can_by_sh[ii].get_end());
			
			fill_sw_can_by_components (struct_id, db_session, sw_can_by_components_PK_id_counter, tag, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), sheet_antiparallel, bs_of_sw_can_by_sh[ii].get_strand_id(), strand_is_at_edge, bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			sw_can_by_components_PK_id_counter++;
			
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
	}
	
	else //!write_phi_psi_
	{
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii)
		{
			Size sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
			Size strand_is_at_edge = see_whether_strand_is_at_edge	(
																	 pose,
																	 struct_id,
																	 db_session,
																	 bs_of_sw_can_by_sh[ii].get_sheet_id(),
																	 sheet_antiparallel,
																	 bs_of_sw_can_by_sh[ii].get_start(),
																	 bs_of_sw_can_by_sh[ii].get_end());
			
			fill_sw_can_by_components (struct_id, db_session, sw_can_by_components_PK_id_counter, tag, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), sheet_antiparallel, bs_of_sw_can_by_sh[ii].get_strand_id(), strand_is_at_edge, bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			sw_can_by_components_PK_id_counter++;			
		}
	}
// <end> fill a table 'sw_can_by_components' by secondary_structure_segments
		TR.Info << "<end> fill a table 'sw_can_by_components' by secondary_structure_segments" << endl;

	Size chance_of_being_canonical_sw = 99; // not yet decided fate whether this could be canonical sandwich or not
	
// <begin> update sheet_connecting_loops
	// get_distinct(sw_can_by_sh_id)
	utility::vector1<Size> vec_sw_can_by_sh_id =  get_vec_sw_can_by_sh_id(struct_id, db_session);

	// per each sandwich_by_sheet_id
	for(Size ii=1; ii<=vec_sw_can_by_sh_id.size(); ++ii)
	{
		Size size_sw_can_by_components_PK_id =
		get_size_sw_can_by_components_PK_id(
											struct_id,
											db_session,
											vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
											);

		bool chance_of_being_sandwich_before_checking_other_strand_check = true;
		Size former_start_res = 0; //temporary
		for(Size jj=1; (jj<=size_sw_can_by_components_PK_id-1) && (chance_of_being_sandwich_before_checking_other_strand_check); ++jj) // this jj is used for counting purpose only, this 'for' loop iterates only for connecting strands
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

			// <begin> check whether there is other strand which is not part of current dealing two sheets
			for(Size kk=start_res+1; (kk<=next_start_res-1) && (chance_of_being_sandwich_before_checking_other_strand_check); ++kk)
			{
				char res_ss( dssp_pose.secstruct( kk ) ) ;
				if (res_ss == 'E')
				{
					bool contains_other_strands = see_whether_other_strands_are_contained(struct_id,
																						  db_session,
																						  vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
																						  kk //res_num
																						  );
					if (contains_other_strands)
					{
						chance_of_being_sandwich_before_checking_other_strand_check = false;
					}
				}
			}
			// <end> check whether there is other strand which is not part of current dealing two sheets
				//TR.Info << "chance_of_being_sandwich_before_checking_other_strand_check: " << chance_of_being_sandwich_before_checking_other_strand_check << endl;

			if (!chance_of_being_sandwich_before_checking_other_strand_check)
			{
				delete_this_sw_can_by_sh_id(struct_id,
											db_session,
											vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
											);
				chance_of_being_canonical_sw = 0;
				break; // break jj 'for' loop
			}
			else
			{
				chance_of_being_canonical_sw = 1;
			}

			if (sheet_id_of_start_res == sheet_id_of_next_start_res)
				//TR.Info << "it connects as intra-sheet way" << endl;
			{
				update_sheet_con(
								 struct_id,
								 db_session,
								 sw_can_by_components_PK_id_counter,
								 tag,
								 vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
								 true,
								 intra_sheet_con_id_counter,
								 inter_sheet_con_id_counter,
								 start_res+1, // new start_res for intra_sheet_con
								 next_start_res-1 // new end_res for intra_sheet_con
								 );

				sw_can_by_components_PK_id_counter++;
				intra_sheet_con_id_counter++;
			}
			else //TR.Info << "it connects as inter-sheet way" << endl;
			{
				update_sheet_con(
								 struct_id,
								 db_session,
								 sw_can_by_components_PK_id_counter,
								 tag,
								 vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
								 false, //bool intra_sheet_con
								 intra_sheet_con_id_counter,
								 inter_sheet_con_id_counter,
								 start_res+1, // new start_res for inter_sheet_con
								 next_start_res-1 // new end_res for inter_sheet_con
								 );

				sw_can_by_components_PK_id_counter++;
				inter_sheet_con_id_counter++;
			}
		}
	}	// per each sandwich_by_sheet_id
// <end> update sheet_connecting_loops
		 
		TR.Info << "chance_of_being_canonical_sw: " << chance_of_being_canonical_sw << endl;
	// <begin> write chain_B_resNum to a file
	if (write_chain_B_resnum_ && (chance_of_being_canonical_sw == 1))
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

		TR.Info << "<Done> for this pdb including extraction of sandwich" << endl;

	return 0;
} //SandwichFeatures::report_features


} //namespace strand_assembly
} //namespace features
} //namespace protocols
