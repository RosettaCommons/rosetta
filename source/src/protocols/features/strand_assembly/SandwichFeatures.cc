
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
/// @brief Extract and analyze beta-sandwich features
/// @author Doo Nam Kim (based on Tim Jacobs' helix_assembly)
/// @overview
///		@ task 0: Determine whether we deal with given pdb file
///		@ task 1: Identify all beta-strands
///		@ task 2: Identify all beta-sheets with these strands
///		@ task 3: Identify all beta-sandwiches with these sheets
///			@ task 3-1: Merge sheets to each other if possible
///			@ task 3-2: Make beta-sandwiches with sheets that are ideal only
///				@ task 3-2-1: Exclude if this_sheet_is_surrounded_by_more_than_1_other_sheet
///				@ task 3-2-2: Exclude sheets that are too close to each other
///				@ task 3-2-3: Exclude sheets that are too distant to each other
///				@ task 3-2-4: Exclude sheets that do not face each other
///					@	task 3-2-4-1: Exclude sheets that do not face each other by an angle with two terminal residues and one central residue
///					@	task 3-2-4-2: Exclude sheets that do not face each other by an angle with four terminal residues in two edge strands
///			@ task 3-3: Write AA distribution
///			@ task 3-4: Test canonical sandwich test
///				@ task 3-4-1: Canonical sandwiches need to have low number of helix or strand residues in any loop (beta-hairpin-loop or inter-sheet-loop)
///				@ task 3-4-2: Canonical sandwiches need to not have same-direction-strands as connecting two beta-sheets
///				@ task 3-4-3: Canonical sandwiches should not be beta-barrel obviously
///		@ task 4: Write beta-sandwiches that passed canonical tests into database
///				@ task 4-1: Write hairpin_loop and inter-sheet loop
///				@ task 4-2: Write starting_loop and endng_loop
///				@ task 4-3: Write ratio_hydrophobic_philic/net_charge
///				@ task 4-4: Write total size of sandwich

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh> // for dssp application
#include <core/pose/PDBInfo.hh> // maybe for PDBInfoCOP
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh> // ScoreFunction.hh seems required for compilation of InterfaceAnalyzerMover.hh
#include <core/scoring/ScoreFunctionFactory.hh> // maybe needed for "getScoreFunction" ?

//External
#include <boost/uuid/uuid.hpp>

//Devel
#include <protocols/features/strand_assembly/SandwichFeatures.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <numeric/xyz.functions.hh> // for torsion calculations
#include <utility/vector1.hh> // for utility::vector1<Column> primary_key_columns;

//for vector
#include <numeric/xyzVector.hh>
#include <core/id/NamedAtomID.hh>

//C library
#include <math.h> // for round, floor, ceil, trunc, sqrt

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
#include <algorithm>	// for avg,min,max
#include <fstream>
#include <iostream>
#include <numeric>
//#include <stdio.h>     //for remove( ) and rename( )
#include <stdlib.h> // for abs()
#include <vector>

template <typename T, size_t N> const T* mybegin(const T (&a)[N]) { return a; }    
template <typename T, size_t N> const T* myend  (const T (&a)[N]) { return a+N; }
// reference:	http://stackoverflow.com/questions/9874802/how-can-i-get-the-max-or-min-value-in-a-vector-c


// exception handling
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.hh> // for SASA

//DSSP
#include <core/scoring/dssp/Dssp.hh>

// for parse_my_tag
#include <basic/datacache/DataMap.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


static basic::Tracer TR("protocols.features.strand_assembly.SandwichFeatures");

namespace protocols {
namespace features {
namespace strand_assembly {

// for parse_my_tag
using utility::tag::TagCOP;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

using core::id::NamedAtomID;
using numeric::xyzVector;

SandwichFeatures::SandwichFeatures()
{
	//TR << "A SandwichFeatures constructor called" << endl;
}

utility::vector1<std::string>
SandwichFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
} //features_reporter_dependencies()

void
SandwichFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

/****** <begin> writing sheet ******/
	// Columns
	// id of sheet
	//unique
	Column sheet_PK_id	("sheet_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column sheet_id	("sheet_id",	new DbInteger(), true /*could be null*/, false /*no autoincrement*/);
		// <history> changed into 'null-possible' because of sw_by_components

	Column sheet_antiparallel	("sheet_antiparallel",	new DbText(), true /* could be null at first, eventually it will not be null though*/, false /*no autoincrement*/);
		// A: antiparallel
		// P_or_mix: parallel or mix of antiparallel and parallel
	Column num_of_sheets_that_surround_this_sheet	("num_of_sheets_that_surround_this_sheet",	new DbInteger(), true /* could be null */, false /*no autoincrement*/);
	
	// unique key of original PDB file
	Column struct_id             ("struct_id",              new DbBigInt(),    false /*not null*/, false /*don't autoincrement*/);

	// ForeignKey
	Column segment_id ("segment_id",	new DbInteger(), false /*not null*/, false /*don't autoincrement*/);

	// Schema - sheet
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sh;
	primary_key_columns_sh.push_back(struct_id);
	primary_key_columns_sh.push_back(sheet_PK_id);

	Schema sheet("sheet",  PrimaryKey(primary_key_columns_sh));

	// add column which is neither PrimaryKey nor ForeignKey
	sheet.add_column(sheet_id);
	sheet.add_column(sheet_antiparallel);
	sheet.add_column(num_of_sheets_that_surround_this_sheet);

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
	//Column tag	("tag",	new DbText(), false /*not null*/, false /*no autoincrement*/);
	Column tag	("tag",	new DbText(), true /*could be null*/, false /*no autoincrement*/);

	// could be redundant
	Column sw_can_by_sh_id	("sw_can_by_sh_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);
	//Column sw_can_by_sh_id	("sw_can_by_sh_id",	new DbInteger(), true /*could be null*/, false /*no autoincrement*/);

	// could be redundant
	Column strand_num	("strand_num",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// Schema - sw_can_by_sh_id
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sw_can_by_sh;
	primary_key_columns_sw_can_by_sh.push_back(struct_id);
	primary_key_columns_sw_can_by_sh.push_back(sw_can_by_sh_PK_id);

	Schema sw_can_by_sh("sw_can_by_sh",  PrimaryKey(primary_key_columns_sw_can_by_sh));

	// add column which is neither PrimaryKey nor ForeignKey
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

	// add column which is neither PrimaryKey nor ForeignKey
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

	Column strand_edge	("strand_edge",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
		// edge strand
		// core strand

	Column num_strands_in_each_sw	("num_strands_in_each_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column num_edge_strands_in_each_sw	("num_edge_strands_in_each_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	
	Column intra_sheet_con_id	("intra_sheet_con_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column inter_sheet_con_id	("inter_sheet_con_id",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	
	Column loop_kind	("loop_kind",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
		// starting_loop
		// hairpin_loop (intra-sheet loop)
		// inter_sheet_loop
		// ending_loop
	
	Column LR	("LR",	new DbText(), true /* could be null*/, false /*no autoincrement*/);

	Column canonical_LR	("canonical_LR",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	// T, -> true, canonical chiral
	// F, -> false, non-canonical chiral
	// U, -> uncertain, this loop-size with this condition has no definite canonical chiral reference in the first place!

	Column turn_type	("turn_type",	new DbText(), true /* could be null*/, false /*no autoincrement*/);

	Column i_AA	("i_AA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column i_p1_AA	("i_p1_AA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column i_p2_AA	("i_p2_AA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column i_p3_AA	("i_p3_AA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);

	Column canonical_turn_AA	("canonical_turn_AA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);

	Column PA_by_preceding_E	("PA_by_preceding_E",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column PA_by_following_E	("PA_by_following_E",	new DbText(), true /* could be null*/, false /*no autoincrement*/);

	Column cano_PA	("cano_PA",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	// T, -> true, canonical PA
	// F, -> false, non-canonical PA
	// U, -> uncertain, this loop-size with this condition has no definite canonical PA reference in the first place!

	Column heading_direction	("heading_direction",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column parallel_EE	("parallel_EE",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column cano_parallel_EE	("cano_parallel_EE",	new DbText(), true /* could be null*/, false /*no autoincrement*/);	
	Column component_size	("component_size",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column A  ("A", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column C  ("C", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column D  ("D", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column E  ("E", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column F  ("F", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column G  ("G", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column H  ("H", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column I  ("I", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column K  ("K", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column L  ("L", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column M  ("M", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column N  ("N", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column P  ("P", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column Q  ("Q", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column R  ("R", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column S  ("S", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column T  ("T", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column V  ("V", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column W  ("W", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column Y  ("Y", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column A_core_heading  ("A_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column A_surface_heading  ("A_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column C_core_heading  ("C_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column C_surface_heading  ("C_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);


	Column D_core_heading  ("D_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column D_surface_heading  ("D_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column E_core_heading  ("E_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column E_surface_heading  ("E_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);


	Column F_core_heading  ("F_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column F_surface_heading  ("F_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column G_core_heading  ("G_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column G_surface_heading  ("G_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column H_core_heading  ("H_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column H_surface_heading  ("H_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column I_core_heading  ("I_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column I_surface_heading  ("I_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column K_core_heading  ("K_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column K_surface_heading  ("K_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column L_core_heading  ("L_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column L_surface_heading  ("L_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column M_core_heading  ("M_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column M_surface_heading  ("M_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column N_core_heading  ("N_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column N_surface_heading  ("N_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);


	Column P_core_heading  ("P_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column P_surface_heading  ("P_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column Q_core_heading  ("Q_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column Q_surface_heading  ("Q_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column R_core_heading  ("R_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column R_surface_heading  ("R_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);



	Column S_core_heading  ("S_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column S_surface_heading  ("S_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column T_core_heading  ("T_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column T_surface_heading  ("T_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);


	Column V_core_heading  ("V_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column V_surface_heading  ("V_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column W_core_heading  ("W_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column W_surface_heading  ("W_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column Y_core_heading  ("Y_core_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);
	Column Y_surface_heading  ("Y_surface_heading", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);


	Column H_percentage ("H_percentage", new DbReal(), true /*could be null*/, false /*don't autoincrement*/);
	Column E_percentage ("E_percentage", new DbReal(), true /*could be null*/, false /*don't autoincrement*/);
	Column L_percentage ("L_percentage", new DbReal(), true /*could be null*/, false /*don't autoincrement*/);

	Column number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands ("number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands ("number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands", new DbInteger(), true /*could be null*/, false /*don't autoincrement*/);

	Column residue_begin("residue_begin", new DbInteger(), false /*not null*/, false /*don't autoincrement*/);
	Column residue_end  ("residue_end", new DbInteger(), false /*not null*/, false /*don't autoincrement*/);



	Column number_of_hydrophobic_res	("number_of_hydrophobic_res",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_hydrophilic_res	("number_of_hydrophilic_res",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column number_of_CGP	("number_of_CGP",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column ratio_hydrophobic_philic_of_sw_in_percent	("ratio_hydrophobic_philic_of_sw_in_percent",	new DbReal(), true /* could be null*/, false /*no autoincrement*/);

	Column number_of_RK_in_sw	("number_of_RK_in_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_DE_in_sw	("number_of_DE_in_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column net_charge_of_sw	("net_charge_of_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column number_of_inward_pointing_W_in_sw	("number_of_inward_pointing_W_in_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_inward_pointing_L_in_core_strands_in_sw	("number_of_inward_pointing_L_in_core_strands_in_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_inward_pointing_W_in_core_strands_in_sw	("number_of_inward_pointing_W_in_core_strands_in_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_inward_pointing_Y_in_core_strands_in_sw	("number_of_inward_pointing_Y_in_core_strands_in_sw",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);
	Column avg_dihedral_angle_between_core_strands_across_facing_sheets	("avg_dihedral_angle_between_core_strands_across_facing_sheets",	new DbReal(), true /* could be null*/, false /*no autoincrement*/);
	Column sw_res_size	("sw_res_size",	new DbInteger(), true /* could be null*/, false /*no autoincrement*/);

	Column multimer_is_suspected	("multimer_is_suspected",	new DbText(), true /* could be null*/, false /*no autoincrement*/);
	Column avg_b_factor_CB_at_each_component	("avg_b_factor_CB_at_each_component",	new DbReal(), true /* could be null*/, false /*no autoincrement*/);

	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sw_by_components;
	primary_key_columns_sw_by_components.push_back(struct_id);
	primary_key_columns_sw_by_components.push_back(sw_by_components_PK_id);

	// table name
	Schema sw_by_components("sw_by_components",  PrimaryKey(primary_key_columns_sw_by_components));

	// add column which is neither PrimaryKey nor ForeignKey
	sw_by_components.add_column(tag);
	sw_by_components.add_column(sw_can_by_sh_id);
	sw_by_components.add_column(sheet_id);

	sw_by_components.add_column(sheet_antiparallel);
	sw_by_components.add_column(sw_by_components_bs_id);
	sw_by_components.add_column(strand_edge);
	sw_by_components.add_column(num_strands_in_each_sw);
	sw_by_components.add_column(num_edge_strands_in_each_sw);

	sw_by_components.add_column(intra_sheet_con_id);
	sw_by_components.add_column(inter_sheet_con_id);
	sw_by_components.add_column(loop_kind); // better to be located right after intra_sheet_con_id/inter_sheet_con_id for better readability

	sw_by_components.add_column(LR);
	sw_by_components.add_column(canonical_LR);

	sw_by_components.add_column(turn_type);
	sw_by_components.add_column(i_AA);
	sw_by_components.add_column(i_p1_AA);
	sw_by_components.add_column(i_p2_AA);
	sw_by_components.add_column(i_p3_AA);

	sw_by_components.add_column(canonical_turn_AA);

	sw_by_components.add_column(PA_by_preceding_E);
	sw_by_components.add_column(PA_by_following_E);
	sw_by_components.add_column(cano_PA);
	sw_by_components.add_column(heading_direction); // nega,posi,meet ...
	sw_by_components.add_column(parallel_EE);
	sw_by_components.add_column(cano_parallel_EE);

	sw_by_components.add_column(A);
	sw_by_components.add_column(C);
	sw_by_components.add_column(D);
	sw_by_components.add_column(E);
	sw_by_components.add_column(F);
	sw_by_components.add_column(G);
	sw_by_components.add_column(H);
	sw_by_components.add_column(I);
	sw_by_components.add_column(K);
	sw_by_components.add_column(L);
	sw_by_components.add_column(M);
	sw_by_components.add_column(N);

	sw_by_components.add_column(P);

	sw_by_components.add_column(Q);

	sw_by_components.add_column(R);
	sw_by_components.add_column(S);
	sw_by_components.add_column(T);

	sw_by_components.add_column(V);

	sw_by_components.add_column(W);
	
	sw_by_components.add_column(Y);

	sw_by_components.add_column(A_core_heading);
	sw_by_components.add_column(A_surface_heading);
	sw_by_components.add_column(C_core_heading);
	sw_by_components.add_column(C_surface_heading);


	sw_by_components.add_column(D_core_heading);
	sw_by_components.add_column(D_surface_heading);
	sw_by_components.add_column(E_core_heading);
	sw_by_components.add_column(E_surface_heading);


	sw_by_components.add_column(F_core_heading);
	sw_by_components.add_column(F_surface_heading);

	sw_by_components.add_column(G_core_heading);
	sw_by_components.add_column(G_surface_heading);


	sw_by_components.add_column(H_core_heading);
	sw_by_components.add_column(H_surface_heading);
	sw_by_components.add_column(I_core_heading);
	sw_by_components.add_column(I_surface_heading);

	sw_by_components.add_column(K_core_heading);
	sw_by_components.add_column(K_surface_heading);

	sw_by_components.add_column(L_core_heading);
	sw_by_components.add_column(L_surface_heading);

	sw_by_components.add_column(M_core_heading);
	sw_by_components.add_column(M_surface_heading);

	sw_by_components.add_column(N_core_heading);
	sw_by_components.add_column(N_surface_heading);

	sw_by_components.add_column(P_core_heading);
	sw_by_components.add_column(P_surface_heading);

	sw_by_components.add_column(Q_core_heading);
	sw_by_components.add_column(Q_surface_heading);

	sw_by_components.add_column(R_core_heading);
	sw_by_components.add_column(R_surface_heading);

	sw_by_components.add_column(S_core_heading);
	sw_by_components.add_column(S_surface_heading);
	sw_by_components.add_column(T_core_heading);
	sw_by_components.add_column(T_surface_heading);

	sw_by_components.add_column(V_core_heading);
	sw_by_components.add_column(V_surface_heading);
	sw_by_components.add_column(W_core_heading);
	sw_by_components.add_column(W_surface_heading);

	sw_by_components.add_column(Y_core_heading);
	sw_by_components.add_column(Y_surface_heading);


	sw_by_components.add_column(H_percentage);
	sw_by_components.add_column(E_percentage);
	sw_by_components.add_column(L_percentage);

	sw_by_components.add_column(number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands);
	sw_by_components.add_column(number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands);

	sw_by_components.add_column(number_of_hydrophobic_res);	//	A,V,I,L,M,F,Y,W
	sw_by_components.add_column(number_of_hydrophilic_res);	//	R,H,K,D,E,S,T,N,Q
	sw_by_components.add_column(number_of_CGP);	//	C,G,P
	sw_by_components.add_column(ratio_hydrophobic_philic_of_sw_in_percent);	//	(no_hydrophobic/no_hydrophilic)*100

	sw_by_components.add_column(number_of_RK_in_sw);	//	R,K
	sw_by_components.add_column(number_of_DE_in_sw);	//	D,E
	sw_by_components.add_column(net_charge_of_sw);
	sw_by_components.add_column(number_of_inward_pointing_W_in_sw);
	sw_by_components.add_column(number_of_inward_pointing_L_in_core_strands_in_sw);
	sw_by_components.add_column(number_of_inward_pointing_W_in_core_strands_in_sw);
	sw_by_components.add_column(number_of_inward_pointing_Y_in_core_strands_in_sw);
	sw_by_components.add_column(avg_dihedral_angle_between_core_strands_across_facing_sheets);
	sw_by_components.add_column(sw_res_size);
	sw_by_components.add_column(multimer_is_suspected);
	sw_by_components.add_column(avg_b_factor_CB_at_each_component);
	sw_by_components.add_column(component_size);

	
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





/****** <begin> writing rkde_in_strands ******/

	// Columns

	// unique_primary_key
	Column rkde_in_strands_PK_id	("rkde_in_strands_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// may not be unique
	Column residue_number	("residue_number",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);
	Column residue_type	("residue_type",	new DbText(), false /*not null*/, false /*no autoincrement*/);

	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_rkde_in_strands;
	primary_key_columns_rkde_in_strands.push_back(struct_id);
	primary_key_columns_rkde_in_strands.push_back(rkde_in_strands_PK_id);

	// table name
	Schema rkde_in_strands("rkde_in_strands",  PrimaryKey(primary_key_columns_rkde_in_strands));

	// add column which is neither PrimaryKey nor ForeignKey
	rkde_in_strands.add_column(tag);
	rkde_in_strands.add_column(sw_can_by_sh_id);

	rkde_in_strands.add_column(residue_number);
	rkde_in_strands.add_column(residue_type);
	rkde_in_strands.add_column(heading_direction);

	// ForeignKey
	rkde_in_strands.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
		// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	rkde_in_strands.write(db_session);
/****** <end> rkde_in_strands ******/



/****** <begin> writing rkde ******/

	// Columns

	// unique_primary_key
	Column rkde_PK_id	("rkde_PK_id",	new DbInteger(), false /*not null*/, false /*no autoincrement*/);

	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_rkde;
	primary_key_columns_rkde.push_back(struct_id);
	primary_key_columns_rkde.push_back(rkde_PK_id);

	// table name
	Schema rkde("rkde",  PrimaryKey(primary_key_columns_rkde));

	// add column which is neither PrimaryKey nor ForeignKey
	rkde.add_column(tag);

	rkde.add_column(residue_number);
	rkde.add_column(residue_type);

	// ForeignKey
	rkde.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
		// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	rkde.write(db_session);
/****** <end> rkde ******/

}


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
SandwichFeatures::get_full_strands(
   StructureID struct_id,
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
	while(res.next())
	{
		Size strand_id,     residue_begin,   residue_end;
		res >> strand_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(residue_begin, residue_end));
	}
	return all_strands;
}


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
SandwichFeatures::get_full_strands_from_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	string select_string =
	"SELECT\n"
	"	sss.residue_begin,\n"
	"	sss.residue_end\n"
	"FROM\n"
	"	secondary_structure_segments as sss, \n"
	"	sheet as sh \n"
	"WHERE\n"
	"	sss.struct_id = sh.struct_id \n"
	"	AND sss.dssp = 'E'\n"
	"	AND sss.struct_id = ? \n"
	"	AND sh.segment_id = sss.segment_id \n"
	"	AND sh.sheet_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,	db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sheet_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next())
	{
		Size residue_begin,   residue_end;
		res >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(residue_begin, residue_end));
	}
	return all_strands;
}


void
SandwichFeatures::report_number_of_electrostatic_interactions_of_residues(
	string	tag,
	StructureID	struct_id,
	sessionOP	db_session,
	Pose const & pose,
	string	dssp_code,
	string	heading_direction)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);

	string ElectroStatic_file_name;
	if	(dssp_code	==	"all_dssp")
	{
		ElectroStatic_file_name = pdb_file_name + "_electrostatic_interactions_of_all_residues.txt";
	}
	else // dssp_code = "E"
	{
		if (heading_direction == "surface")
		{
			ElectroStatic_file_name = pdb_file_name + "_electrostatic_interactions_of_surface_residues_in_a_strand.txt";
		}
		else // (heading_direction == "all_direction")
		{
			ElectroStatic_file_name = pdb_file_name + "_electrostatic_interactions_of_all_residues_in_a_strand.txt";
		}
	}

	ofstream ElectroStatic_file;
		
	ElectroStatic_file.open(ElectroStatic_file_name.c_str());

	string number_of_attractions_title;
	string head_attrac = "attrac_by_centroid_w_";
	number_of_attractions_title.append(head_attrac);

	string number_of_repulsions_title;
	string head_repul = "repul_by_centroid_w_";
	number_of_repulsions_title.append(head_repul);

	std::ostringstream cutoff;
	cutoff << distance_cutoff_for_electrostatic_interactions_;
	std::string cutoff_str = cutoff.str();
	number_of_attractions_title.append(cutoff_str);
	number_of_repulsions_title.append(cutoff_str);

	string tail	=	"_A";
	number_of_attractions_title.append(tail);
	number_of_repulsions_title.append(tail);

	ElectroStatic_file << "resNum	type	"	<<	number_of_attractions_title	<<	"	"	<<	number_of_repulsions_title << "	attrac_minus_repul	salt_bridge	CC_bridge	NO_bridge	sum_of_salt_CC_NO_bridges	longer_range_ion_pair"	<<	endl;

	utility::vector1<Size>	vector_of_unique_distinct_sw_ids	=	get_distinct_sw_id_from_sw_by_components_table	(struct_id,	db_session);

	std::vector<int> vec_number_of_attractions_by_centroid;
	std::vector<int> vec_number_of_repulsions_by_centroid;
	std::vector<int> vec_net_attrac_by_centroid;
	std::vector<int> vec_number_of_salt_bridges;
	std::vector<int> vec_number_of_CC_bridges;
	std::vector<int> vec_number_of_NO_bridges;
	std::vector<int> vec_sum_of_salt_CC_NO_bridges;
	std::vector<int> vec_number_of_longer_range_ion_pair;

	for(Size sw_ii=1; sw_ii<=vector_of_unique_distinct_sw_ids.size(); sw_ii++) // per each beta-sandwich
	{
		utility::vector1<Size>	vector_of_residue_num_of_rkde	=
			retrieve_residue_num_of_rkde(
				struct_id,
				db_session,
				vector_of_unique_distinct_sw_ids[sw_ii],	//sw_can_by_sh_id
				dssp_code,
				heading_direction
				);
		for(Size residue_i=1; residue_i<=vector_of_residue_num_of_rkde.size(); residue_i++)
		{
			Size	residue_num	=	vector_of_residue_num_of_rkde[residue_i];
			ElectroStatic_file	<< residue_num	<<	"	"	<<	pose.residue_type(residue_num).name3();

			numeric::xyzVector< core::Real > xyz_of_centroid_of_RKDE;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_R;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_R;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_of_K;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_DE;
			numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_DE;

			if (pose.residue_type(residue_num).name3() == "ARG")
			{
				// <begin> calculate centroid position
				xyz_of_centroid_of_RKDE.x()	=
						(pose.residue(residue_num).atom(" NE ").xyz().x()
					+	pose.residue(residue_num).atom(" CZ ").xyz().x()
					+	pose.residue(residue_num).atom(" NH1").xyz().x()
					+	pose.residue(residue_num).atom(" NH2").xyz().x())/4;

				xyz_of_centroid_of_RKDE.y()	=
						(pose.residue(residue_num).atom(" NE ").xyz().y()
					+	pose.residue(residue_num).atom(" CZ ").xyz().y()
					+	pose.residue(residue_num).atom(" NH1").xyz().y()
					+	pose.residue(residue_num).atom(" NH2").xyz().y())/4;

				xyz_of_centroid_of_RKDE.z()	=
						(pose.residue(residue_num).atom(" NE ").xyz().z()
					+	pose.residue(residue_num).atom(" CZ ").xyz().z()
					+	pose.residue(residue_num).atom(" NH1").xyz().z()
					+	pose.residue(residue_num).atom(" NH2").xyz().z())/4;
				// <end> calculate centroid position

				xyz_of_terminal_atom_1_of_R =	pose.residue(residue_num).atom(" NH1").xyz();
				xyz_of_terminal_atom_2_of_R =	pose.residue(residue_num).atom(" NH2").xyz();
			}
			else	if (pose.residue_type(residue_num).name3() == "LYS")
			{
				xyz_of_centroid_of_RKDE =	pose.residue(residue_num).atom(" NZ ").xyz();
				xyz_of_terminal_atom_of_K =	pose.residue(residue_num).atom(" NZ ").xyz();
			}
			else	if (pose.residue_type(residue_num).name3() == "ASP")
			{
				// <begin> calculate centroid position
				xyz_of_centroid_of_RKDE.x()	=
						(pose.residue(residue_num).atom(" CG ").xyz().x()
					+	pose.residue(residue_num).atom(" OD1").xyz().x()
					+	pose.residue(residue_num).atom(" OD2").xyz().x())/3;

				xyz_of_centroid_of_RKDE.y()	=
						(pose.residue(residue_num).atom(" CG ").xyz().y()
					+	pose.residue(residue_num).atom(" OD1").xyz().y()
					+	pose.residue(residue_num).atom(" OD2").xyz().y())/3;

				xyz_of_centroid_of_RKDE.z()	=
						(pose.residue(residue_num).atom(" CG ").xyz().z()
					+	pose.residue(residue_num).atom(" OD1").xyz().z()
					+	pose.residue(residue_num).atom(" OD2").xyz().z())/3;
				// <end> calculate centroid position

				xyz_of_terminal_atom_1_of_DE =	pose.residue(residue_num).atom(" OD1").xyz();
				xyz_of_terminal_atom_2_of_DE =	pose.residue(residue_num).atom(" OD2").xyz();
			}
			else	//if (pose.residue_type(residue_num).name3() == "GLU")
			{
				// <begin> calculate centroid position
				xyz_of_centroid_of_RKDE.x()	=
						(pose.residue(residue_num).atom(" CD ").xyz().x()
					+	pose.residue(residue_num).atom(" OE1").xyz().x()
					+	pose.residue(residue_num).atom(" OE2").xyz().x())/3;

				xyz_of_centroid_of_RKDE.y()	=
						(pose.residue(residue_num).atom(" CD ").xyz().y()
					+	pose.residue(residue_num).atom(" OE1").xyz().y()
					+	pose.residue(residue_num).atom(" OE2").xyz().y())/3;

				xyz_of_centroid_of_RKDE.z()	=
						(pose.residue(residue_num).atom(" CD ").xyz().z()
					+	pose.residue(residue_num).atom(" OE1").xyz().z()
					+	pose.residue(residue_num).atom(" OE2").xyz().z())/3;
				// <end> calculate centroid position

				xyz_of_terminal_atom_1_of_DE =	pose.residue(residue_num).atom(" OE1").xyz();
				xyz_of_terminal_atom_2_of_DE =	pose.residue(residue_num).atom(" OE2").xyz();
			}

			int	number_of_attractions_by_centroid	=	0;
			int	number_of_repulsions_by_centroid	=	0;

			Size	number_of_salt_bridge	=	0;
			Size	number_of_C_C_bridge	=	0;
			Size	number_of_N_O_bridge	=	0;
			Size	number_of_longer_range_ion_pair	=	0;

			for(Size other_residue_i=1; other_residue_i<=vector_of_residue_num_of_rkde.size(); other_residue_i++)
			{
				if	(other_residue_i	==	residue_i)
				{
					continue;
				}

				Size	other_residue_num	=	vector_of_residue_num_of_rkde[other_residue_i];

				if ((other_residue_num - residue_num) < primary_seq_distance_cutoff_for_electrostatic_interactions_)
					//	other_residue_num should be always higher than	residue_num
				{
					continue;
				}

				// <begin> check whether "other_residue" has low atom position uncertainty
				pose::PDBInfoCOP info = pose.pdb_info();
				Real B_factor_of_CB = info->temperature( other_residue_num, 5 ); // '5' atom will be 'H' for Gly
				if (B_factor_of_CB	>	CB_b_facor_cutoff_for_electrostatic_interactions_)
				{
					continue;	// this	"other_residue" has too high atom position uncertainty, so let's not use this residue when counting number of electrostatic interactions
				}
				// <end> check whether "other_residue" has low atom position uncertainty

				numeric::xyzVector< core::Real > xyz_of_centroid_of_other_RKDE;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_other_R;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_other_R;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_of_other_K;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_1_of_other_DE;
				numeric::xyzVector< core::Real > xyz_of_terminal_atom_2_of_other_DE;

				if (pose.residue_type(other_residue_num).name3() == "ARG")
				{
					// <begin> calculate centroid position
					xyz_of_centroid_of_other_RKDE.x()	=
						(pose.residue(other_residue_num).atom(" NE ").xyz().x()
					+	pose.residue(other_residue_num).atom(" CZ ").xyz().x()
					+	pose.residue(other_residue_num).atom(" NH1").xyz().x()
					+	pose.residue(other_residue_num).atom(" NH2").xyz().x())/4;

					xyz_of_centroid_of_other_RKDE.y()	=
						(pose.residue(other_residue_num).atom(" NE ").xyz().y()
					+	pose.residue(other_residue_num).atom(" CZ ").xyz().y()
					+	pose.residue(other_residue_num).atom(" NH1").xyz().y()
					+	pose.residue(other_residue_num).atom(" NH2").xyz().y())/4;

					xyz_of_centroid_of_other_RKDE.z()	=
						(pose.residue(other_residue_num).atom(" NE ").xyz().z()
					+	pose.residue(other_residue_num).atom(" CZ ").xyz().z()
					+	pose.residue(other_residue_num).atom(" NH1").xyz().z()
					+	pose.residue(other_residue_num).atom(" NH2").xyz().z())/4;
					// <end> calculate centroid position

					xyz_of_terminal_atom_1_of_other_R =	pose.residue(other_residue_num).atom(" NH1").xyz();
					xyz_of_terminal_atom_2_of_other_R =	pose.residue(other_residue_num).atom(" NH2").xyz();
				}
				else	if (pose.residue_type(other_residue_num).name3() == "LYS")
				{
					xyz_of_centroid_of_other_RKDE =	pose.residue(other_residue_num).atom(" NZ ").xyz();
					xyz_of_terminal_atom_of_other_K =	pose.residue(other_residue_num).atom(" NZ ").xyz();
				}
				else	if (pose.residue_type(other_residue_num).name3() == "ASP")
				{
					// <begin> calculate centroid position
					xyz_of_centroid_of_other_RKDE.x()	=
						(pose.residue(other_residue_num).atom(" CG ").xyz().x()
					+	pose.residue(other_residue_num).atom(" OD1").xyz().x()
					+	pose.residue(other_residue_num).atom(" OD2").xyz().x())/3;

					xyz_of_centroid_of_other_RKDE.y()	=
						(pose.residue(other_residue_num).atom(" CG ").xyz().y()
					+	pose.residue(other_residue_num).atom(" OD1").xyz().y()
					+	pose.residue(other_residue_num).atom(" OD2").xyz().y())/3;

					xyz_of_centroid_of_other_RKDE.z()	=
						(pose.residue(other_residue_num).atom(" CG ").xyz().z()
					+	pose.residue(other_residue_num).atom(" OD1").xyz().z()
					+	pose.residue(other_residue_num).atom(" OD2").xyz().z())/3;
					// <end> calculate centroid position

					xyz_of_terminal_atom_1_of_other_DE =	pose.residue(other_residue_num).atom(" OD1").xyz();
					xyz_of_terminal_atom_2_of_other_DE =	pose.residue(other_residue_num).atom(" OD2").xyz();
				}
				else	//if (pose.residue_type(other_residue_num).name3() == "GLU")
				{
					// <begin> calculate centroid position
					xyz_of_centroid_of_other_RKDE.x()	=
						(pose.residue(other_residue_num).atom(" CD ").xyz().x()
					+	pose.residue(other_residue_num).atom(" OE1").xyz().x()
					+	pose.residue(other_residue_num).atom(" OE2").xyz().x())/3;

					xyz_of_centroid_of_other_RKDE.y()	=
						(pose.residue(other_residue_num).atom(" CD ").xyz().y()
					+	pose.residue(other_residue_num).atom(" OE1").xyz().y()
					+	pose.residue(other_residue_num).atom(" OE2").xyz().y())/3;

					xyz_of_centroid_of_other_RKDE.z()	=
						(pose.residue(other_residue_num).atom(" CD ").xyz().z()
					+	pose.residue(other_residue_num).atom(" OE1").xyz().z()
					+	pose.residue(other_residue_num).atom(" OE2").xyz().z())/3;
					// <end> calculate centroid position

					xyz_of_terminal_atom_1_of_other_DE =	pose.residue(other_residue_num).atom(" OE1").xyz();
					xyz_of_terminal_atom_2_of_other_DE =	pose.residue(other_residue_num).atom(" OE2").xyz();
				}

				Real distance_between_centroid = xyz_of_centroid_of_RKDE.distance(xyz_of_centroid_of_other_RKDE);

				if	(pose.residue_type(residue_num).name3() == "ARG")
				{
					if	(pose.residue_type(other_residue_num).name3() == "ASP"	||	pose.residue_type(other_residue_num).name3() == "GLU")
					{
						if (distance_between_centroid < distance_cutoff_for_electrostatic_interactions_)
						{
							number_of_attractions_by_centroid++;
						}
						Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_1_of_R.distance(xyz_of_terminal_atom_1_of_other_DE);
						Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_1_of_R.distance(xyz_of_terminal_atom_2_of_other_DE);
						Real distance_3_between_terminal_atoms = xyz_of_terminal_atom_2_of_R.distance(xyz_of_terminal_atom_1_of_other_DE);
						Real distance_4_between_terminal_atoms = xyz_of_terminal_atom_2_of_R.distance(xyz_of_terminal_atom_2_of_other_DE);
						if	(distance_between_centroid < 7.0)
						{
							if	(distance_between_centroid < 4.0)
							{
								if	(distance_1_between_terminal_atoms	< 4.0 ||	distance_2_between_terminal_atoms	< 4.0
								||	distance_3_between_terminal_atoms	< 4.0 ||	distance_4_between_terminal_atoms	< 4.0)
								{
									number_of_salt_bridge++;
								}
								else
								{
									number_of_C_C_bridge++;
								}
							}
							else
							{
								if	(distance_1_between_terminal_atoms	< 4.0 ||	distance_2_between_terminal_atoms	< 4.0
								||	distance_3_between_terminal_atoms	< 4.0 ||	distance_4_between_terminal_atoms	< 4.0)
								{
									number_of_N_O_bridge++;
								}
								else
								{
									number_of_longer_range_ion_pair++;
								}
							}
						}
					}
					else
					{
						if (distance_between_centroid < distance_cutoff_for_electrostatic_interactions_)
						{
							number_of_repulsions_by_centroid++;
						}
					}
				}
				else	if	(pose.residue_type(residue_num).name3() == "LYS")
				{
					if	(pose.residue_type(other_residue_num).name3() == "ASP"	||	pose.residue_type(other_residue_num).name3() == "GLU")
					{
						if (distance_between_centroid < distance_cutoff_for_electrostatic_interactions_)
						{
							number_of_attractions_by_centroid++;
						}
						Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_of_K.distance(xyz_of_terminal_atom_1_of_other_DE);
						Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_of_K.distance(xyz_of_terminal_atom_2_of_other_DE);
						if	(distance_between_centroid < 7.0)
						{
							if	(distance_between_centroid < 4.0)
							{
								if	(distance_1_between_terminal_atoms < 4.0 ||	distance_2_between_terminal_atoms	< 4.0)
								{
									number_of_salt_bridge++;
								}
								else
								{
									number_of_C_C_bridge++;
								}
							}
							else
							{
								if	(distance_1_between_terminal_atoms < 4.0 ||	distance_2_between_terminal_atoms	< 4.0)
								{
									number_of_N_O_bridge++;
								}
								else
								{
									number_of_longer_range_ion_pair++;
								}
							}
						}
					}
					else
					{
						if (distance_between_centroid < distance_cutoff_for_electrostatic_interactions_)
						{
							number_of_repulsions_by_centroid++;
						}
					}
				}
				else // (pose.residue_type(residue_num).name3() == "ASP") || (pose.residue_type(residue_num).name3() == "GLU")
				{
					if	((pose.residue_type(other_residue_num).name3() == "ASP") || (pose.residue_type(other_residue_num).name3() == "GLU"))
					{
						if (distance_between_centroid < distance_cutoff_for_electrostatic_interactions_)
						{
							number_of_repulsions_by_centroid++;
						}
					}
					else // (pose.residue_type(other_residue_num).name3() == "ARG") || (pose.residue_type(other_residue_num).name3() == "LYS")
					{
						if (distance_between_centroid < distance_cutoff_for_electrostatic_interactions_)
						{
							number_of_attractions_by_centroid++;
						}
						if	(pose.residue_type(other_residue_num).name3() == "ARG")
						{
							Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_1_of_other_R.distance(xyz_of_terminal_atom_1_of_DE);
							Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_1_of_other_R.distance(xyz_of_terminal_atom_2_of_DE);
							Real distance_3_between_terminal_atoms = xyz_of_terminal_atom_2_of_other_R.distance(xyz_of_terminal_atom_1_of_DE);
							Real distance_4_between_terminal_atoms = xyz_of_terminal_atom_2_of_other_R.distance(xyz_of_terminal_atom_2_of_DE);

							if	(distance_between_centroid < 7.0)
							{
								if	(distance_between_centroid < 4.0)
								{
									if	(distance_1_between_terminal_atoms < 4.0 ||	distance_2_between_terminal_atoms	< 4.0	||	distance_3_between_terminal_atoms	< 4.0	 ||	distance_4_between_terminal_atoms	< 4.0)
									{
										number_of_salt_bridge++;
									}
									else
									{
										number_of_C_C_bridge++;
									}
								}
								else
								{
									if	(distance_1_between_terminal_atoms < 4.0 ||	distance_2_between_terminal_atoms	< 4.0	||	distance_3_between_terminal_atoms	< 4.0	 ||	distance_4_between_terminal_atoms	< 4.0)
									{
										number_of_N_O_bridge++;
									}
									else
									{
										number_of_longer_range_ion_pair++;
									}
								}
							}
						}
						else //	(pose.residue_type(other_residue_num).name3() == "LYS")
						{
							Real distance_1_between_terminal_atoms = xyz_of_terminal_atom_of_other_K.distance(xyz_of_terminal_atom_1_of_DE);
							Real distance_2_between_terminal_atoms = xyz_of_terminal_atom_of_other_K.distance(xyz_of_terminal_atom_2_of_DE);
							if	(distance_between_centroid < 7.0)
							{
								if	(distance_between_centroid < 4.0)
								{
									if	(distance_1_between_terminal_atoms < 4.0 ||	distance_2_between_terminal_atoms	< 4.0)
									{
										number_of_salt_bridge++;
									}
									else
									{
										number_of_C_C_bridge++;
									}
								}
								else
								{
									if	(distance_1_between_terminal_atoms < 4.0 ||	distance_2_between_terminal_atoms	< 4.0)
									{
										number_of_N_O_bridge++;
									}
									else
									{
										number_of_longer_range_ion_pair++;
									}
								}
							}
						}
					}
				}
			}

			vec_number_of_attractions_by_centroid.push_back(number_of_attractions_by_centroid);
			vec_number_of_repulsions_by_centroid.push_back(number_of_repulsions_by_centroid);
			vec_net_attrac_by_centroid.push_back(number_of_attractions_by_centroid-number_of_repulsions_by_centroid);
			vec_number_of_salt_bridges.push_back(number_of_salt_bridge);
			vec_number_of_CC_bridges.push_back(number_of_C_C_bridge);
			vec_number_of_NO_bridges.push_back(number_of_N_O_bridge);
			vec_sum_of_salt_CC_NO_bridges.push_back(number_of_salt_bridge+number_of_C_C_bridge+number_of_N_O_bridge);
			vec_number_of_longer_range_ion_pair.push_back(number_of_longer_range_ion_pair);

			ElectroStatic_file	<<	"	"	<< number_of_attractions_by_centroid	<< "	"	<<	number_of_repulsions_by_centroid	<<	"	"	<<	number_of_attractions_by_centroid-number_of_repulsions_by_centroid	<<	"	"	<<	number_of_salt_bridge	<<	"	"	<<	number_of_C_C_bridge	<<	"	"	<<	number_of_N_O_bridge	<<	"	"	<<	number_of_salt_bridge+number_of_C_C_bridge+number_of_N_O_bridge	<<	"	"	<<		number_of_longer_range_ion_pair	<<	endl;
		} // per each residue
	}

	// report avg (min~max)
	// warning: this avg (min~max) assumes that I deal with 1 sandwich per 1 pdb file
	float	avg_attrac	=	std::accumulate(vec_number_of_attractions_by_centroid.begin(),		vec_number_of_attractions_by_centroid.end(),	0)	/	static_cast<float>(vec_number_of_attractions_by_centroid.size());

	float	avg_repul	=	std::accumulate(vec_number_of_repulsions_by_centroid.begin(),		vec_number_of_repulsions_by_centroid.end(),	0)	/	static_cast<float>(vec_number_of_repulsions_by_centroid.size());

	float	avg_net_attrac	=	std::accumulate(vec_net_attrac_by_centroid.begin(),		vec_net_attrac_by_centroid.end(),	0)	/	static_cast<float>(vec_net_attrac_by_centroid.size());

	float	avg_salt_bridges	=	std::accumulate(vec_number_of_salt_bridges.begin(),		vec_number_of_salt_bridges.end(),	0)	/	static_cast<float>(vec_number_of_salt_bridges.size());

	float	avg_CC_bridges	=	std::accumulate(vec_number_of_CC_bridges.begin(),		vec_number_of_CC_bridges.end(),	0)	/	static_cast<float>(vec_number_of_CC_bridges.size());

	float	avg_NO_bridges	=	std::accumulate(vec_number_of_NO_bridges.begin(),		vec_number_of_NO_bridges.end(),	0)	/	static_cast<float>(vec_number_of_NO_bridges.size());

	float	avg_sum_of_salt_CC_NO_bridges	=	std::accumulate(vec_sum_of_salt_CC_NO_bridges.begin(),		vec_sum_of_salt_CC_NO_bridges.end(),	0)	/	static_cast<float>(vec_sum_of_salt_CC_NO_bridges.size());

	float	avg_number_of_longer_range_ion_pair	=	std::accumulate(vec_number_of_longer_range_ion_pair.begin(),		vec_number_of_longer_range_ion_pair.end(),	0)	/	static_cast<float>(vec_number_of_longer_range_ion_pair.size());


	int	min_attrac;
	int	max_attrac;
	if (vec_number_of_attractions_by_centroid.size()	==	0)
	{
		min_attrac	=	0;
		max_attrac	=	0;
	}
	else
	{
		min_attrac	=	*std::min_element(vec_number_of_attractions_by_centroid.begin(),	vec_number_of_attractions_by_centroid.end());
		max_attrac	=	*std::max_element(vec_number_of_attractions_by_centroid.begin(),	vec_number_of_attractions_by_centroid.end());
	}


	int	min_repul;
	int	max_repul;
	if (vec_number_of_repulsions_by_centroid.size()	==	0)
	{
		min_repul	=	0;
		max_repul	=	0;
	}
	else
	{
		min_repul	=	*std::min_element(vec_number_of_repulsions_by_centroid.begin(),	vec_number_of_repulsions_by_centroid.end());
		max_repul	=	*std::max_element(vec_number_of_repulsions_by_centroid.begin(),	vec_number_of_repulsions_by_centroid.end());
	}

	int	min_net_attrac;
	int	max_net_attrac;
	if	(vec_net_attrac_by_centroid.size()	==	0)
	{
		min_net_attrac	=	0;
		max_net_attrac	=	0;
	}
	else
	{
		min_net_attrac	=	*std::min_element(vec_net_attrac_by_centroid.begin(),	vec_net_attrac_by_centroid.end());
		max_net_attrac	=	*std::max_element(vec_net_attrac_by_centroid.begin(),	vec_net_attrac_by_centroid.end());
	}

	int	min_salt_bridge;
	int	max_salt_bridge;
	if (vec_number_of_salt_bridges.size()	==	0)
	{
		min_salt_bridge	=	0;
		max_salt_bridge	=	0;
	}
	else
	{
		min_salt_bridge	=	*std::min_element(vec_number_of_salt_bridges.begin(),	vec_number_of_salt_bridges.end());
		max_salt_bridge	=	*std::max_element(vec_number_of_salt_bridges.begin(),	vec_number_of_salt_bridges.end());
	}

	int	min_CC_bridge;
	int	max_CC_bridge;
	if (vec_number_of_CC_bridges.size()	==	0)
	{
		min_CC_bridge	=	0;
		max_CC_bridge	=	0;
	}
	else
	{
		min_CC_bridge	=	*std::min_element(vec_number_of_CC_bridges.begin(),	vec_number_of_CC_bridges.end());
		max_CC_bridge	=	*std::max_element(vec_number_of_CC_bridges.begin(),	vec_number_of_CC_bridges.end());
	}

	int	min_NO_bridge;
	int	max_NO_bridge;
	if (vec_number_of_NO_bridges.size()	==	0)
	{
		min_NO_bridge	=	0;
		max_NO_bridge	=	0;
	}
	else
	{
		min_NO_bridge	=	*std::min_element(vec_number_of_NO_bridges.begin(),	vec_number_of_NO_bridges.end());
		max_NO_bridge	=	*std::max_element(vec_number_of_NO_bridges.begin(),	vec_number_of_NO_bridges.end());
	}

	int	min_sum_of_salt_CC_NO_bridges;
	int	max_sum_of_salt_CC_NO_bridges;
	if	(vec_sum_of_salt_CC_NO_bridges.size()	==	0)
	{
		min_sum_of_salt_CC_NO_bridges	=	0;
		max_sum_of_salt_CC_NO_bridges	=	0;
	}
	else
	{
		min_sum_of_salt_CC_NO_bridges	=	*std::min_element(vec_sum_of_salt_CC_NO_bridges.begin(),	vec_sum_of_salt_CC_NO_bridges.end());
		max_sum_of_salt_CC_NO_bridges	=	*std::max_element(vec_sum_of_salt_CC_NO_bridges.begin(),	vec_sum_of_salt_CC_NO_bridges.end());
	}

	int	min_number_of_longer_range_ion_pair;
	int	max_number_of_longer_range_ion_pair;
	if (vec_number_of_longer_range_ion_pair.size()	==	0)
	{
		min_number_of_longer_range_ion_pair	=	0;
		max_number_of_longer_range_ion_pair	=	0;
	}
	else
	{
		min_number_of_longer_range_ion_pair	=	*std::min_element(vec_number_of_longer_range_ion_pair.begin(),	vec_number_of_longer_range_ion_pair.end());
		max_number_of_longer_range_ion_pair	=	*std::max_element(vec_number_of_longer_range_ion_pair.begin(),	vec_number_of_longer_range_ion_pair.end());
	}

	ElectroStatic_file	<<	"avg	(min~max)	";
	ElectroStatic_file	<< round_to_float(avg_attrac)	<<	" ("	<<	min_attrac	<<	"~"	<<	max_attrac	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_repul)	<<	" ("	<<	min_repul	<<	"~"	<<	max_repul	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_net_attrac)	<<	" ("	<<	min_net_attrac	<<	"~"	<<	max_net_attrac	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_salt_bridges)	<<	" ("	<<	min_salt_bridge	<<	"~"	<<	max_salt_bridge	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_CC_bridges)	<<	" ("	<<	min_CC_bridge	<<	"~"	<<	max_CC_bridge	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_NO_bridges)	<<	" ("	<<	min_NO_bridge	<<	"~"	<<	max_NO_bridge	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_sum_of_salt_CC_NO_bridges)	<<	" ("	<<	min_sum_of_salt_CC_NO_bridges	<<	"~"	<<	max_sum_of_salt_CC_NO_bridges	<<	")	";
	ElectroStatic_file	<< round_to_float(avg_number_of_longer_range_ion_pair)	<<	" ("	<<	min_number_of_longer_range_ion_pair	<<	"~"	<<	max_number_of_longer_range_ion_pair	<<	")	"	<<	endl;

	ElectroStatic_file.close();
}


//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
bool
SandwichFeatures::check_whether_strand_i_is_in_sheet(
	StructureID struct_id,
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
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	sh.sheet_id,\n"
	"	sh.segment_id,\n"
	"	sss.residue_begin,\n"
	"	sss.residue_end \n"
	"FROM\n"
	"	sheet as sh,\n"
	"	secondary_structure_segments AS sss\n"
	"WHERE\n"
	"	sh.segment_id = sss.segment_id \n"
	"	AND sss.dssp = 'E' \n"
	"	AND sh.struct_id = sss.struct_id \n"
	"	AND sh.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next())
	{
		Size sheet_id, segment_id,	residue_begin,	residue_end;
		res >> sheet_id >> segment_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(sheet_id, residue_begin, residue_end));
	}
	return all_strands;
} //get_current_strands_in_sheet


utility::vector1<Size>	
SandwichFeatures::get_distinct_sheet_id_from_sheet_table(
	StructureID struct_id,
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
} //get_distinct_sheet_id_from_sheet_table


utility::vector1<Size>	
SandwichFeatures::get_distinct_sw_id_from_sw_by_components_table(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	distinct sw_can_by_sh_id\n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> all_distinct_sw_ids;
	while(res.next())
	{
		Size distinct_sw_id;
		res >> distinct_sw_id;
		all_distinct_sw_ids.push_back(distinct_sw_id);
	}
	return all_distinct_sw_ids;
} //get_distinct_sw_id_from_sw_by_components_table




//get_max_sheet_id
Size
SandwichFeatures::get_max_sheet_id(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	max(sheet_id) \n"
	"FROM\n"
	"	sheet \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sheet_id != 99999);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size max_sheet_id;
	while(res.next())
	{
		res >> max_sheet_id;
	}
	return max_sheet_id;
} //get_max_sheet_id


//update_sheet_id
Size
SandwichFeatures::update_sheet_id(
	StructureID struct_id,
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
	StructureID struct_id,
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


//update_num_of_sheets_that_surround_this_sheet
void	
SandwichFeatures::update_num_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id,
	Size num_of_sheets_that_surround_this_sheet)
{
	string select_string =
	"UPDATE sheet set num_of_sheets_that_surround_this_sheet = ?	"
	"WHERE\n"
	"	(sheet_id = ?) \n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));

	select_statement.bind(1,	num_of_sheets_that_surround_this_sheet);
	select_statement.bind(2,	sheet_id);
	select_statement.bind(3,	struct_id);

	basic::database::safely_write_to_database(select_statement);
} //update_num_of_sheets_that_surround_this_sheet


//update_rkde_in_strands
void	
SandwichFeatures::update_rkde_in_strands(
	StructureID struct_id,
	sessionOP db_session,
	Size	rkde_in_strands_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size residue_number,
	string	residue_type,
	string	heading_direction)
{
	string insert =	"INSERT INTO rkde_in_strands (struct_id, rkde_in_strands_PK_id, tag, sw_can_by_sh_id, residue_number, residue_type,	heading_direction)  VALUES (?,?,?,?,	?,?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	rkde_in_strands_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	sw_can_by_sh_id);
	insert_stmt.bind(5,	residue_number);
	insert_stmt.bind(6,	residue_type);
	insert_stmt.bind(7,	heading_direction);	//	surface or core
	basic::database::safely_write_to_database(insert_stmt);

} //update_rkde_in_strands


//update_rkde
void	
SandwichFeatures::update_rkde(
	StructureID struct_id,
	sessionOP db_session,
	Size	rkde_PK_id_counter,
	string tag,
	Size residue_number,
	string	residue_type)
{
	string insert =	"INSERT INTO rkde (struct_id, rkde_PK_id, tag, residue_number, residue_type)  VALUES (?,?,?,	?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	rkde_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	residue_number);
	insert_stmt.bind(5,	residue_type);
	basic::database::safely_write_to_database(insert_stmt);

} //update_rkde


//get_num_of_sheets_that_surround_this_sheet
Size	
SandwichFeatures::get_num_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	string select_string =
	"SELECT\n"
	"	num_of_sheets_that_surround_this_sheet \n"
	"FROM\n"
	"	sheet \n"
	"WHERE\n"
	"	(sheet_id=?) \n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	sheet_id);
	select_statement.bind(2,	struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_of_sheets_that_surround_this_sheet;
	while(res.next())
	{
		res >> num_of_sheets_that_surround_this_sheet;
	}
	return num_of_sheets_that_surround_this_sheet;
} //get_num_of_sheets_that_surround_this_sheet



//get_chain_B_resNum
utility::vector1<Size>
SandwichFeatures::get_chain_B_resNum(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	// <begin> get first sheet_id
	utility::vector1<Size> vec_sheet_id =  get_vec_distinct_sheet_id(struct_id, db_session,	sw_can_by_sh_id);
	Size first_sh_id	=	vec_sheet_id[1];
	// <end> get first sheet_id

	string select_string =
	"SELECT\n"
	"	r.resNum \n"
	"FROM\n"
	"	sw_by_components AS sw, \n"
	"	residues AS r \n"
	"WHERE\n"
	"	(sw.sw_can_by_sh_id=?) \n"
	"	AND	(sw.sheet_id=?) \n"
	"	AND (r.resNum >= sw.residue_begin AND r.resNum <= sw.residue_end) \n"
	"	AND (sw.struct_id = r.struct_id) \n"
	"	AND (sw.struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	sw_can_by_sh_id);
	select_statement.bind(2,	first_sh_id);
	select_statement.bind(3,	struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> vec_chain_B_resNum;
	while(res.next())
	{
		Size chain_B_resNum;
		res >> chain_B_resNum;
		vec_chain_B_resNum.push_back(chain_B_resNum);
	}
	return vec_chain_B_resNum;
} //get_chain_B_resNum


//get_tag
string
SandwichFeatures::get_tag(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	tag \n"
	"FROM\n"
	"	structures \n"
	"WHERE\n"
	"	struct_id = ?;";

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
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
//		TR << "get_num_strands_in_this_sheet" << endl;
//	time_t start_time = time(NULL);
	string select_string =
	"SELECT\n"
	"	count(*) \n"
	"FROM\n"
	"	sheet \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sheet_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sheet_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_strands;
	while(res.next())
	{
		res >> num_strands;
	}
//	time_t end_time = time(NULL);
//		TR.Info << "Finished in " << (end_time - start_time) << " seconds." << endl;
	return num_strands;
} //get_num_strands_in_this_sheet


//	prepare_to_fill_sw_by_components
//	<role>	It retrieves all beta-strands of sandwich_candidate_by_sheets, it does not make sandwich_by_components
utility::vector1<SandwichFragment>
SandwichFeatures::prepare_to_fill_sw_by_components(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT \n"
	"	sw_sh.sw_can_by_sh_id AS sw_can_by_sh_id, \n"
	"	sh.sheet_id AS sheet_id, \n"
	"	sss.segment_id AS sw_by_components_bs_id, \n"
	"	sss.residue_begin AS	residue_begin, \n"
	"	sss.residue_end AS	residue_end \n"
	"FROM  \n"
	"	secondary_structure_segments AS sss, \n"
	"	sheet AS sh, \n"
	"	sw_can_by_sh AS sw_sh\n"
	"WHERE \n"
	"	(sss.struct_id = ?) \n"
	"	AND (sss.struct_id = sh.struct_id) \n"
	"	AND (sss.struct_id = sw_sh.struct_id) \n"
	"	AND (sss.dssp = 'E') \n"
	"	AND (sw_sh.sheet_id = sh.sheet_id) \n"
	"	AND (sh.segment_id = sss.segment_id);";


	/* (for sqlite3 try)
	SELECT
		sss.struct_id AS struct_id,
		sw_sh.sw_can_by_sh_id AS sw_can_by_sh_id, 
		sh.sheet_id AS sheet_id, 
		sss.segment_id AS sw_by_components_bs_id, 
		sss.residue_begin AS	residue_begin, 
		sss.residue_end AS	residue_end 
	FROM  
		secondary_structure_segments AS sss, 
		sheet AS sh, 
		sw_can_by_sh AS sw_sh
	WHERE 
		(sss.struct_id = sh.struct_id) 
		AND (sss.struct_id = sw_sh.struct_id) 
		AND (sss.dssp = 'E') 
		AND (sw_sh.sheet_id = sh.sheet_id) 
		AND (sh.segment_id = sss.segment_id);
	*/

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
SandwichFeatures::absolute_vec (numeric::xyzVector<Real> vector)
{
	Real absolute_vec = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	return absolute_vec;
}


// I tried to slightly vary the method of Nobu's recent nature paper
string
SandwichFeatures::check_LR ( // check L/R chirality of sidechain
	Pose & dssp_pose,
	Size residue_begin,
	Size residue_end)
{
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
	Real cosine_theta = dot_product_with_a_b / (absolute_vec(cross_product_u_v))*(absolute_vec(vector_a_b));
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
	Real	cosine_theta_between_v_and_preceding_E = dot_product_with_v_and_preceding_E / (absolute_vec(preceding_E_vec_a_b))*(absolute_vec(vector_v));

	Real	dot_product_with_v_and_following_E = dot_product( following_E_vec_a_b, vector_v );
	Real	cosine_theta_between_v_and_following_E = dot_product_with_v_and_following_E / (absolute_vec(following_E_vec_a_b))*(absolute_vec(vector_v));

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
	Real	cosine_theta_between_sidechains = dot_product_with_preceding_E_and_following_E / (absolute_vec(preceding_E_vec_a_b))*(absolute_vec(following_E_vec_a_b));
		// referred http://www.mvps.org/directx/articles/math/dot/index.htm
		// cosine 80 = 0.173

	if (check_N_to_C_direction_by_ == "PE") // default
	{
		Real	dot_product_with_v_and_preceding_E = dot_product( preceding_E_vec_a_b, vector_v );
		Real	cosine_theta_between_v_and_preceding_E = dot_product_with_v_and_preceding_E / (absolute_vec(preceding_E_vec_a_b))*(absolute_vec(vector_v));

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
		Real	cosine_theta_between_v_and_following_E = dot_product_with_v_and_following_E / (absolute_vec(following_E_vec_a_b))*(absolute_vec(vector_v));

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
	bool antiparalell // if 'false', find a sheet in parallel way
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
				if (antiparalell)	// find a sheet in an anti-parallel way
				{
					if (i_resnum+2 > strand_i.get_end() || j_resnum-2 < strand_j.get_start())
					{
						continue; // I want to extract strand_pairs within only valid ranges
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
					return 1; //  may have a kinkness or not, but these strands can be part of one sheet
				}
			} // strand_i.get_size() >= 3 && strand_j.get_size() >= 3)
		} // for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
	} // for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	
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
	StructureID struct_id,
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

			if (return_of_find_sheet_antiparallel == 999)	// since these two strands are too distant to each other, there is virtually no chance to be sheet!
			{
				break; // since these two strands are too distant to each other, there is virtually no chance to be sheet!
			}

			if (!return_of_find_sheet_antiparallel)
			{
				return_of_find_sheet_parallel = find_sheet (pose, temp_strand_i, temp_strand_j, false);
			}

			if (return_of_find_sheet_parallel == 999)	// since these two strands are too distant to each other, there is virtually no chance to be sheet!
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
	StructureID struct_id,
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


// see_whether_sheet_is_antiparallel
// <Definition> In order to be an antiparallel sheet, all strands in the sheet should be anti-parallel to each other
// <Overall plot>
	// 1. Get central residues of each edge strands (save -999 as resnum if the strand is at core)
	// 2. Get array of sum of all distances from central residues
	// 3. Get the largest distance and its residue index
	// 4. Search the nearest strand from the current strand
string	
SandwichFeatures::see_whether_sheet_is_antiparallel(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size i_sheet)
{
	utility::vector1<SandwichFragment> strands_from_i = get_full_strands_from_sheet(struct_id, db_session, i_sheet); // struct_id, db_session, sheet_id


	// <begin> Get central residues of edge strands
	vector<Real> vector_of_central_resnum; // array of central residues
	vector<Size> vector_of_not_representative_edge_strands; // both core and very_short_edge strands
	vector<Size> vector_of_short_edge_strands;

	for(Size i=1; i<=strands_from_i.size(); ++i)
	{
		string strand_is_at_edge = is_this_strand_at_edge	(
										pose,
										struct_id,
										db_session,
										i_sheet,
										strands_from_i[i].get_start(),
										strands_from_i[i].get_end());
		if (strand_is_at_edge == "core")
		{
			vector_of_central_resnum.push_back(-999); // won't be used, but needed
			vector_of_not_representative_edge_strands.push_back(i);
			continue;
		}// I ignore a core (unrepresentative here) strand

		if (strand_is_at_edge == "short_edge")
		{
			vector_of_central_resnum.push_back(-99); // won't be used, but needed
			vector_of_not_representative_edge_strands.push_back(i);
			vector_of_short_edge_strands.push_back(i);
			continue;
		}// I ignore a core (unrepresentative here) strand


		Real to_be_rounded = (strands_from_i[i].get_start() + strands_from_i[i].get_end())/(2.0);
		Size cen_resnum = round_to_Size(to_be_rounded);
		vector_of_central_resnum.push_back(cen_resnum);
	}
	// <end> Get central residues of edge strands


	// <begin> Get array of sum of all distances from central residues
	vector<Real> vec_dis_sum_from_cen_resnum; // array of sum of all distances from central residues
	for(Size ii=0; ii<=strands_from_i.size()-1; ++ii)
	{
		Real dis_from_ii = 0;
		
		if (vector_of_central_resnum[ii] == -999 || vector_of_central_resnum[ii] == -99) // I ignore a unrepresentative strand
		{
			dis_from_ii = -999;
			vec_dis_sum_from_cen_resnum.push_back(dis_from_ii);
			continue;
		}

		for (Size jj=0; jj<=strands_from_i.size()-1; ++jj)
		{				
			if (vector_of_central_resnum[jj] == -999 || vector_of_central_resnum[jj] == -99) // I ignore a unrepresentative strand
			{
				continue;
			}
			Real dis_CA_CA	=	pose.residue(static_cast<Size>(vector_of_central_resnum[ii])).atom("CA").xyz().distance(pose.residue(static_cast<Size>(vector_of_central_resnum[jj])).atom("CA").xyz());
			dis_from_ii = dis_from_ii + dis_CA_CA;
		}
		vec_dis_sum_from_cen_resnum.push_back(dis_from_ii);
	}
	// <end> get array of sum of all distances from central residues


	// <begin> Get the largest distance and its residue index (which should indicate the furthermost strand in a sheet)
	Size res_index_having_the_largest_dis = 0;
	Real largest_dis = -99;

	for(Size ii=0; ii<=strands_from_i.size()-1; ++ii)
	{
		if (vector_of_central_resnum[ii] == -999 || vector_of_central_resnum[ii] == -99) // I ignore a unrepresentative strand
		{
			continue; // unrepresentative_strands
		}
		if (largest_dis < vec_dis_sum_from_cen_resnum[ii])
		{
			largest_dis = vec_dis_sum_from_cen_resnum[ii];
			res_index_having_the_largest_dis = ii;
		}
	}
	// <end> Get the largest distance and its residue index (which should indicate the furthermost strand in a sheet)


	Size former_res_index_nearest_strand = res_index_having_the_largest_dis; //	just for the first step of 'while' loop
	Size former_res_index_having_the_largest_dis = res_index_having_the_largest_dis; //	just for the first step of 'while' loop
	
	Size maximum_possible_number_of_antiparallel_matches = strands_from_i.size() -1 -vector_of_short_edge_strands.size();
	Size count_antiparellel = 0;
	while (count_antiparellel < maximum_possible_number_of_antiparallel_matches)
	{
		SandwichFragment current_strand(strands_from_i[former_res_index_nearest_strand+1].get_start(), strands_from_i[former_res_index_nearest_strand+1].get_end());
			
		Real shortest_dis_inter_strand = 999;
		Size res_index_nearest_strand = 999;


		//<begin> search the nearest strand from the current strand
		for(Size ii=0; ii<=strands_from_i.size()-1; ++ii)
		{
			if (ii == former_res_index_nearest_strand || ii == former_res_index_having_the_largest_dis)
			{
				continue;
			}
			if (vector_of_central_resnum[ii] == -99) // I ignore a very short edge strand (but don't ignore a core strand)
			{
				continue; 
			}
			SandwichFragment nearest_strand_from_the_current_strand(strands_from_i[ii+1].get_start(), strands_from_i[ii+1].get_end());
			Real inter_strand_avg_dis = get_avg_dis_strands (pose, current_strand, nearest_strand_from_the_current_strand);
			if (inter_strand_avg_dis < shortest_dis_inter_strand)
			{
				shortest_dis_inter_strand = inter_strand_avg_dis;
				res_index_nearest_strand = ii;
			}
		}
		//<end> search the nearest strand from the current strand

		if (	res_index_nearest_strand == 999)
		{
			return "undetermined"; // this sheet may be consisted with very short strands only like sheet_id=7 in 1MSP
		}

		SandwichFragment nearest_strand_from_the_current_strand(strands_from_i[res_index_nearest_strand+1].get_start(), strands_from_i[res_index_nearest_strand+1].get_end());
		Size return_of_find_sheet = find_sheet (pose, current_strand, nearest_strand_from_the_current_strand, true);
			
		if ( return_of_find_sheet == 1 )
		{
			count_antiparellel++;
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
	return false; // OK, these two strand_pairs are distant to each other enough
} //SandwichFeatures::check_strand_too_closeness


// check whether these sheets are too close, the closeness is checked for every possible distances
// <where to use	now?> get_current_bs_id_and_closest_edge_bs_id_in_different_sheet
// <where to use	now?> in main, check whether strands are too distant to each other
// <where to use	now?> see_whether_sheet_is_antiparallel
Real
SandwichFeatures::get_avg_dis_strands(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j)
{
	if (strand_i.get_start() == strand_j.get_start()) // strand_i and strand_j are same!
	{
		return 0.0;
	}

	Real sum_dis_CA_CA = 0;
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
			// <tip> I don't need to worry about the possibility of having different distance results depending on directionality of strands since I calculate all possible combinatorial distances (--> confirmed by experiment!)
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
			sum_dis_CA_CA = sum_dis_CA_CA + dis_CA_CA;
		}
	}
	//	TR << "avg_dis_strands: " <<	sum_dis_CA_CA/(strand_i.get_size()*strand_j.get_size()) << endl;
	return sum_dis_CA_CA/(strand_i.get_size()*strand_j.get_size());
} //SandwichFeatures::get_avg_dis_strands


// check whether these sheets are too close, the closeness is checked for every possible distances
Real
SandwichFeatures::get_closest_distance_between_strands(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j)
{
	if (strand_i.get_start() == strand_j.get_start()) // strand_i and strand_j are same!
	{
		return 0.0;
	}

	Real closest_dis_CA_CA = 9999;
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
			// <tip> I don't need to worry about the possibility of having different distance results depending on directionality of strands since I calculate all possible combinatorial distances (--> confirmed by experiment!)
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
			if (closest_dis_CA_CA > dis_CA_CA)
			{
				closest_dis_CA_CA = dis_CA_CA;
			}
		}
	}
	//	TR << "closest_dis_CA_CA: " <<	closest_dis_CA_CA << endl;
	return closest_dis_CA_CA;
} //SandwichFeatures::get_closest_distance_between_strands



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
	if (dis_CA_CA_0 > 40){
		return -999; // these sheets will not be sandwich ever, since these two sheets are too distant!
	}
	if (dis_CA_CA_0 < min_sheet_dis_ || dis_CA_CA_0 > max_sheet_dis_){
		return -99;
	}

	Real dis_CA_CA_1 = pose.residue(i_resnum_1).atom("CA").xyz().distance(pose.residue(j_resnum_1).atom("CA").xyz());
	if (dis_CA_CA_1 < min_sheet_dis_ || dis_CA_CA_1 > max_sheet_dis_){
		return -99;
	}

	Real dis_CA_CA_2 = pose.residue(i_resnum_2).atom("CA").xyz().distance(pose.residue(j_resnum_2).atom("CA").xyz());
	if (dis_CA_CA_2 < min_sheet_dis_ || dis_CA_CA_2 > max_sheet_dis_){
		return -99;
	}

	Real dis_CA_CA_3 = pose.residue(i_resnum_3).atom("CA").xyz().distance(pose.residue(j_resnum_3).atom("CA").xyz());
	if (dis_CA_CA_3 < min_sheet_dis_ || dis_CA_CA_3 > max_sheet_dis_){
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

			else { // paralell
				i_resnum_1 = i_resnum+1;
				j_resnum_1 = j_resnum+1;

				i_resnum_2 = i_resnum+2;
				j_resnum_2 = j_resnum+2;

				i_resnum_3 = i_resnum+3;
				j_resnum_3 = j_resnum+3;

				if (i_resnum_3 > pose.total_residue() || j_resnum_3 > pose.total_residue()){
					continue;
				}
			}

			Real avg_dis_CA_CA = get_avg_dis_CA_CA(pose, i_resnum,	i_resnum_1, i_resnum_2, i_resnum_3, j_resnum, j_resnum_1, j_resnum_2, j_resnum_3);

			if (avg_dis_CA_CA == -999){
				break; // these sheets will not be sandwich ever, since these two sheets are too distant!
			}

			if (avg_dis_CA_CA == -99) { // dis_CA_CA_x < min_sheet_dis_ || dis_CA_CA_x > max_sheet_dis_
				continue;
			}
			return avg_dis_CA_CA;
		} //for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
	} //for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	return -99; // these sheets are not sandwich with these strands
} //SandwichFeatures::check_sw_by_dis


Size
SandwichFeatures::round_to_Size(
					Real x)
{
	Size rounded = static_cast <Size> (floor(x+.5));
	return rounded;
} //round


float
SandwichFeatures::round_to_float(
					float x)
{
	return floor((x	*	10)	+	0.5)	/	10;
} //round

Real
SandwichFeatures::round_to_Real(
					Real x)
{
	Real rounded = floor(x+.2);
	return rounded;
} //round_to_Real


// is_this_strand_at_edge
// <role> This function is used to write edge_strand info into db
// <note> If a strand is within "short_edge", the strand is not an "edge"
string	
SandwichFeatures::is_this_strand_at_edge	(
	Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id,
	Size residue_begin,
	Size residue_end)
{
//		TR << "is_this_strand_at_edge" << endl;
//	time_t start_time = time(NULL);
	//		TR << "residue_begin : " << residue_begin << endl;
//		TR << "residue_end : " << residue_end << endl;

	if (residue_end - residue_begin + 1 < 3)
	{
		return "short_edge"; // Like res_num: 78-79 in [1ten], 2 residues long strand is better to be classified as 'short edge'
	}

	// <begin> see whether this sheet is consisted with two strands only
		utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_id);

		if (strands_from_sheet_i.size() < 3)
		{
			return "edge"; // Since this sheet is constituted with two strands, both are edge!
		}
	// <end> see whether this sheet is consisted with two strands only


	SandwichFragment current_strand(residue_begin, residue_end);
	vector<Real> vec_inter_strand_dis;
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i)
	{
		SandwichFragment temporary_strand(strands_from_sheet_i[i].get_start(), strands_from_sheet_i[i].get_end());
		Real closest_distance_between_strands = get_closest_distance_between_strands (pose, current_strand, temporary_strand);
		vec_inter_strand_dis.push_back(closest_distance_between_strands);
	}
	
	Size size_of_vec_inter_strand_dis = vec_inter_strand_dis.size();
	
	// <begin> exclude self-strand
		Real min_inter_strand_avg_dis = 9999;
		Size index_having_self_strand = 0;
		for(Size i=0; i<=size_of_vec_inter_strand_dis-1; ++i)
		{
			if (min_inter_strand_avg_dis > vec_inter_strand_dis[i])
			{
				min_inter_strand_avg_dis = vec_inter_strand_dis[i];
				index_having_self_strand = i+1; // index of vec_inter_strand_avg_dis starts with '0'
												// while index of strands_from_sheet_i starts with '1'
			}
		}
	// <end> exclude self-strand

	// <begin> find the closest strand from current_strand
		min_inter_strand_avg_dis = 9999;
		Size	index_having_min_dis = 0;
		for(Size i=0; i<=size_of_vec_inter_strand_dis-1; ++i)
		{
			if (
				(i != index_having_self_strand-1)	//	exclude self-strand
				&&	(min_inter_strand_avg_dis > vec_inter_strand_dis[i])
				&&	(strands_from_sheet_i[i+1].get_size() > 2)
				)
			{
				min_inter_strand_avg_dis = vec_inter_strand_dis[i];
				index_having_min_dis = i+1; // index of vec_inter_strand_avg_dis starts with 0 while index of strands_from_sheet_i starts with 1
			}
		}
	// <end> find the closest strand from current_strand

		//		TR << "residue_begin of the closest strand: " << strands_from_sheet_i[index_having_min_dis].get_start() << endl;
		//		TR << "residue_end of the closest strand: " << strands_from_sheet_i[index_having_min_dis].get_end() << endl;

	// <begin> calculate minimum distance between strands
		Real to_be_rounded_i = (strands_from_sheet_i[index_having_min_dis].get_start() + strands_from_sheet_i[index_having_min_dis].get_end())/(2.0);
		Size cen_resnum_of_the_closest_strand = round_to_Size(to_be_rounded_i);

		Real min_inter_strand_dis = 9999;
		
		for(Size strand_i_res = residue_begin;
			strand_i_res <= residue_end;
			strand_i_res++)
		{
			Real dis_CA_CA = pose.residue(strand_i_res).atom("CA").xyz().distance(pose.residue(cen_resnum_of_the_closest_strand).atom("CA").xyz());
			if (min_inter_strand_dis > dis_CA_CA)
			{
				min_inter_strand_dis = dis_CA_CA;
			}
		}
	// <end> calculate minimum distance between strands

		//TR << "min_inter_strand_dis with the closest strand: " << min_inter_strand_dis << endl;

	if (min_inter_strand_dis > min_CA_CA_dis_ && min_inter_strand_dis < max_CA_CA_dis_)
	{
		// <begin> find the 2nd closest strand from current_strand
			min_inter_strand_avg_dis = 9999;
			Size index_having_second_min_dis = 0;
			for(Size i=0; i<=size_of_vec_inter_strand_dis-1; ++i)
			{
				if (
					(i != index_having_self_strand-1)
					&&	(i != index_having_min_dis-1)
					&&	(min_inter_strand_avg_dis > vec_inter_strand_dis[i])
					)
					// && strands_from_sheet_i[i+1].get_size() > 2) disabled as of 09/25/2013 to debug a crash of 9CGT
					// so from now on, a beta-strand is 'core' if it is next to short_edge_strand and within core region and within specified distance.
					//	Still it classifies a beta-strand as "edge" if it is not within specified distance like '1TEN')
				{
					min_inter_strand_avg_dis = vec_inter_strand_dis[i];
					index_having_second_min_dis = i+1; // index of vec_inter_strand_avg_dis starts with 0 while index of strands_from_sheet_i starts with 1
				}
			}
		// <end> find the 2nd closest strand from current_strand

//			TR << "residue_begin of the 2nd closest strand: " << strands_from_sheet_i[index_having_second_min_dis].get_start() << endl;
//			TR << "residue_end of the 2nd closest strand: " << strands_from_sheet_i[index_having_second_min_dis].get_end() << endl;

		// <begin> calculate minimum distance between strands
			to_be_rounded_i = (strands_from_sheet_i[index_having_second_min_dis].get_start() + strands_from_sheet_i[index_having_second_min_dis].get_end())/(2.0);
			Size cental_resnum_of_the_2nd_closest_strand = round_to_Size(to_be_rounded_i);

			min_inter_strand_dis = 9999;
			
			for(Size strand_i_res = residue_begin;
				strand_i_res <= residue_end;
				strand_i_res++)
			{
				Real dis_CA_CA = pose.residue(strand_i_res).atom("CA").xyz().distance(pose.residue(cental_resnum_of_the_2nd_closest_strand).atom("CA").xyz());
				if (min_inter_strand_dis > dis_CA_CA)
				{
					min_inter_strand_dis = dis_CA_CA;
				}
			}
		// <end> calculate minimum distance between strands

//		time_t end_time = time(NULL);
//			TR.Info << "Finished in " << (end_time - start_time) << " seconds." << endl;

		if (min_inter_strand_dis > min_CA_CA_dis_ && min_inter_strand_dis < max_CA_CA_dis_)
		{
			return "core";
		}
	}
	return "edge";
}
// <end> is_this_strand_at_edge


string
SandwichFeatures::is_this_strand_at_edge_by_looking_db(
	StructureID struct_id,
	sessionOP db_session,
	Size residue_begin)
{
	string select_string =
	"SELECT\n"
	"	strand_edge\n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"  AND residue_begin = ? ;";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,residue_begin);
	result res(basic::database::safely_read_from_database(select_statement));
	
	string edge;
	while(res.next())
	{
		res >> edge ;
	}
	return edge;
} //is_this_strand_at_edge_by_looking_db




void SandwichFeatures::process_decoy(
	Pose &dssp_pose,
	core::scoring::ScoreFunction const& scorefxn ) const
{
	scorefxn( dssp_pose );
} // process_decoy


core::scoring::ScoreFunctionOP SandwichFeatures::generate_scorefxn( bool fullatom ) {
	core::scoring::ScoreFunctionOP scorefxn( NULL );
	if ( fullatom )
	{
		scorefxn = core::scoring::getScoreFunction();
	}
	else
	{
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	}
	return scorefxn;
}



//check whether this sheet is constituted with 2 (or less) residues long only
bool
SandwichFeatures::check_whether_this_sheet_is_too_short(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_i){
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


//get_central_residues_in_each_of_two_edge_strands
std::pair<Size, Size>
SandwichFeatures::get_central_residues_in_each_of_two_edge_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i){
	utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_i);

	// get central residue numbers in edge strands
	vector<Size> vector_of_central_residues_in_sheet_i;
	for(Size i=1; i<=strands_from_sheet_i.size(); ++i){
		if (strands_from_sheet_i[i].get_size() <= 2){
			continue;
		}

		string strand_is_at_edge = is_this_strand_at_edge	(
										pose,
										struct_id,
										db_session,
										sheet_i,
										strands_from_sheet_i[i].get_start(),
										strands_from_sheet_i[i].get_end());
		if (strand_is_at_edge == "core" || strand_is_at_edge == "short_edge"){
			continue;
		}// I ignore a unrepresentative strand

		Real to_be_rounded_i = (strands_from_sheet_i[i].get_start() + strands_from_sheet_i[i].get_end())/(2.0);
		Size cen_resnum_i = round_to_Size(to_be_rounded_i);
		vector_of_central_residues_in_sheet_i.push_back(cen_resnum_i);
	}
	Size array_size = vector_of_central_residues_in_sheet_i.size();

	if (array_size == 0) // this sheet maybe beta-barrel like sheet_id = 1 in 1N8O
	{
		return std::make_pair(-99, -99);
	}
	
	// <begin> get sum of distances
	vector<Real> sum_dis_array_i;
	for(Size i=0; i<=array_size-1; ++i){
		Real sum_dis_i_and_j = 0;
		for(Size j=0; (j<=array_size-1); ++j){
			if (i == j){
				continue;
			}
			Real dis_i_and_j = pose.residue(vector_of_central_residues_in_sheet_i[i]).atom("CA").xyz().distance(pose.residue(vector_of_central_residues_in_sheet_i[j]).atom("CA").xyz());
			sum_dis_i_and_j = sum_dis_i_and_j + dis_i_and_j;
		}
		sum_dis_array_i.push_back(sum_dis_i_and_j);
	}
	// <end> get sum of distances


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

	Size terminal_cen_res_pos_1 = vector_of_central_residues_in_sheet_i[index_terminal_cen_res_pos_1];	// index of the first central_residue
	Size terminal_cen_res_pos_2 = vector_of_central_residues_in_sheet_i[index_terminal_cen_res_pos_2];	// index of the second central_residue

	return std::make_pair(terminal_cen_res_pos_1, terminal_cen_res_pos_2);
} //get_central_residues_in_each_of_two_edge_strands


Real
SandwichFeatures::get_shortest_among_4_vals(
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
} //SandwichFeatures::get_shortest_among_4_vals (simple one with just four parameters)


//Select all strand segments reported by the secondary_structure_segments and save them in a vector
utility::vector1<SandwichFragment>
SandwichFeatures::get_all_strands_in_sheet_i(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	string select_string =
	"SELECT\n"
	"	sh.sheet_id,\n"
	"	sh.segment_id,\n"
	"	sss.residue_begin,\n"
	"	sss.residue_end \n"
	"FROM\n"
	"	sheet as sh,\n"
	"	secondary_structure_segments AS sss\n"
	"WHERE\n"
	"	sh.segment_id = sss.segment_id \n"
	"	AND sss.dssp = 'E' \n" // just sanity check
	"	AND sh.struct_id = sss.struct_id \n"
	"	AND sh.sheet_id=? \n"
	"	AND sh.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,sheet_id);
	select_statement.bind(2,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<SandwichFragment> all_strands;
	while(res.next())
	{
		Size sheet_id, segment_id,	residue_begin,	residue_end;
		res >> sheet_id >> segment_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(sheet_id, segment_id,	residue_begin, residue_end));
	}
	return all_strands;
} //get_all_strands_in_sheet_i


// get_start_end_res_num_in_the_longest_strand which is used for judge_facing
utility::vector1<SandwichFragment>
SandwichFeatures::get_start_end_res_num_in_the_longest_strand(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
	// <begin> Identify the longest strand in this sheet
	utility::vector1<SandwichFragment> all_strands_in_sheet_i	=	get_all_strands_in_sheet_i(struct_id,	db_session,	sheet_id);
	Size longest_size_of_strand = 0;
	Size residue_begin_of_the_longest_strand = 0;
	Size residue_end_of_the_longest_strand = 0;
	utility::vector1<SandwichFragment> start_end_res_num_in_longest_strand;
	for(Size i=1; i<=all_strands_in_sheet_i.size() ; ++i)
	{
		Size size_of_this_strand = all_strands_in_sheet_i[i].get_size();
		if (size_of_this_strand > longest_size_of_strand)
		{
			longest_size_of_strand	=	size_of_this_strand;
			residue_begin_of_the_longest_strand = all_strands_in_sheet_i[i].get_start();
			residue_end_of_the_longest_strand = all_strands_in_sheet_i[i].get_end();
		}
	}
	// <end> Identify the longest strand in this sheet

	start_end_res_num_in_longest_strand.push_back(SandwichFragment(residue_begin_of_the_longest_strand, residue_end_of_the_longest_strand));
	
	return start_end_res_num_in_longest_strand;
} //get_start_end_res_num_in_the_longest_strand


int
SandwichFeatures::judge_facing(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i,
	Size sheet_j)
{
	// <begin> check_whether_this_sheet_is_too_short
	bool this_strand_is_too_short = check_whether_this_sheet_is_too_short(
		struct_id,
		db_session,
		sheet_i);

	if (this_strand_is_too_short)
	{
		return 0; // I can't choose two central residues since this sheet is constituted with 2 residues long strands only
	}

	this_strand_is_too_short = check_whether_this_sheet_is_too_short(
		struct_id,
		db_session,
		sheet_j);

	if (this_strand_is_too_short)
	{
		return 0; // I can't choose two central residues since this sheet is constituted with 2 residues long strands only
	}
	// <end> check_whether_this_sheet_is_too_short


	// <begin> get_start_end_res_num_in_the_longest_strand in two sheets
	utility::vector1<SandwichFragment> start_end_res_num_in_the_longest_strand_in_sheet_i =
		get_start_end_res_num_in_the_longest_strand(
			struct_id,
			db_session,
			sheet_i);

	utility::vector1<SandwichFragment> start_end_res_num_in_the_longest_strand_in_sheet_j =
		get_start_end_res_num_in_the_longest_strand(
			struct_id,
			db_session,
			sheet_j);
	// <end> get_start_end_res_num_in_the_longest_strand in two sheets


	// <begin> measure inter-sheet angle to prevent non-facing sheets
	Real angle_with_cen_res;
	if (start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_size()	>	start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_size())
		// Index of SandwichFragment starts with '1' not '0' 
	{
		Real to_be_rounded_i = (start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_start() + start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_end())/(2.0);
		Size cen_resnum_of_smaller_sheet = round_to_Size(to_be_rounded_i);

		Real distance_1 = pose.residue(cen_resnum_of_smaller_sheet).atom("CA").xyz().distance(pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_start()).atom("CA").xyz());
		Real distance_2 = pose.residue(cen_resnum_of_smaller_sheet).atom("CA").xyz().distance(pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_end()).atom("CA").xyz());

		if (distance_1 < distance_2){
			Vector const& first_res_xyz    ( pose.residue(cen_resnum_of_smaller_sheet).xyz("CA") );
			Vector const& middle_res_xyz   ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_start()).xyz("CA") );
			Vector const& third_res_xyz    ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_end()).xyz("CA") );

			angle_with_cen_res = numeric::angle_degrees(first_res_xyz, middle_res_xyz, third_res_xyz);
		}
		else{
			Vector const& first_res_xyz    ( pose.residue(cen_resnum_of_smaller_sheet).xyz("CA") );
			Vector const& middle_res_xyz   ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_end()).xyz("CA") );
			Vector const& third_res_xyz    ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_start()).xyz("CA") );

			angle_with_cen_res = numeric::angle_degrees(first_res_xyz, middle_res_xyz, third_res_xyz);
		}
	}
	else
	{
		Real to_be_rounded_i = (start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_start() + start_end_res_num_in_the_longest_strand_in_sheet_i[1].get_end())/(2.0);
		Size cen_resnum_of_smaller_sheet = round_to_Size(to_be_rounded_i);

		Real distance_1 = pose.residue(cen_resnum_of_smaller_sheet).atom("CA").xyz().distance(pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_start()).atom("CA").xyz());
		Real distance_2 = pose.residue(cen_resnum_of_smaller_sheet).atom("CA").xyz().distance(pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_end()).atom("CA").xyz());

		if (distance_1 < distance_2){
			Vector const& first_res_xyz    ( pose.residue(cen_resnum_of_smaller_sheet).xyz("CA") );
			Vector const& middle_res_xyz   ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_start()).xyz("CA") );
			Vector const& third_res_xyz    ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_end()).xyz("CA") );

			angle_with_cen_res = numeric::angle_degrees(first_res_xyz, middle_res_xyz, third_res_xyz);
		}
		else{
			Vector const& first_res_xyz    ( pose.residue(cen_resnum_of_smaller_sheet).xyz("CA") );
			Vector const& middle_res_xyz   ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_end()).xyz("CA") );
			Vector const& third_res_xyz    ( pose.residue(start_end_res_num_in_the_longest_strand_in_sheet_j[1].get_start()).xyz("CA") );

			angle_with_cen_res = numeric::angle_degrees(first_res_xyz, middle_res_xyz, third_res_xyz);
		}
	}
	// <end> measure inter-sheet angle to prevent non-facing sheets



	if (angle_with_cen_res > max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_)
	{
		return 0; // these two sheets are linear or do not face to each other properly!
	}


	// <begin> identify four terminal central residues
	std::pair<Size, Size>
	two_central_residues_in_two_edge_strands =	get_central_residues_in_each_of_two_edge_strands(
		struct_id,
		db_session,
		pose,
		sheet_i);
		
	int i_ter_cen_1 = two_central_residues_in_two_edge_strands.first;
	int i_ter_cen_2 = two_central_residues_in_two_edge_strands.second;

	if (i_ter_cen_1 == -99 || i_ter_cen_2 == -99)
	{
		return -99;
	}

	two_central_residues_in_two_edge_strands =	get_central_residues_in_each_of_two_edge_strands(
		struct_id,
		db_session,
		pose,
		sheet_j);
		
	int j_ter_cen_1 = two_central_residues_in_two_edge_strands.first;
	int j_ter_cen_2 = two_central_residues_in_two_edge_strands.second;

	if (j_ter_cen_1 == -99 || j_ter_cen_2 == -99)
	{
		return -99;
	}


	Real arr_dis_inter_sheet [4];
	arr_dis_inter_sheet[0] = pose.residue(i_ter_cen_1).atom("CA").xyz().distance(pose.residue(j_ter_cen_1).atom("CA").xyz());
	arr_dis_inter_sheet[1] = pose.residue(i_ter_cen_1).atom("CA").xyz().distance(pose.residue(j_ter_cen_2).atom("CA").xyz());
	arr_dis_inter_sheet[2] = pose.residue(i_ter_cen_2).atom("CA").xyz().distance(pose.residue(j_ter_cen_1).atom("CA").xyz());
	arr_dis_inter_sheet[3] = pose.residue(i_ter_cen_2).atom("CA").xyz().distance(pose.residue(j_ter_cen_2).atom("CA").xyz());

	Real
	shortest_dis_inter_sheet = get_shortest_among_4_vals(arr_dis_inter_sheet);

		// temporary value assignment
	Real angle_1 = 999.9;
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

	if ((angle_1 > min_sheet_angle_by_four_term_cen_res_) &&
		(angle_1 < max_sheet_angle_by_four_term_cen_res_) &&
		(angle_2 > min_sheet_angle_by_four_term_cen_res_) &&
		(angle_2 < max_sheet_angle_by_four_term_cen_res_) &&
		(torsion_i_j > min_sheet_torsion_cen_res_)
		&& (torsion_i_j < max_sheet_torsion_cen_res_))
	{
		return 1; // these two strand_pairs face to each other properly, so constitute a sandwich
	}

	else
	{
		return 0; // these two strand_pairs are linear or do not face to each other properly!
	}

} //SandwichFeatures::judge_facing


Real
SandwichFeatures::calculate_dihedral_w_4_resnums(
	Pose const & pose,	Size res1_sheet_i,		Size res2_sheet_i,		Size res1_sheet_j,		Size res2_sheet_j)
{
	Real arr_dis_inter_sheet [4];
	arr_dis_inter_sheet[0] = pose.residue(res1_sheet_i).atom("CA").xyz().distance(pose.residue(res1_sheet_j).atom("CA").xyz());
	arr_dis_inter_sheet[1] = pose.residue(res1_sheet_i).atom("CA").xyz().distance(pose.residue(res2_sheet_j).atom("CA").xyz());
	arr_dis_inter_sheet[2] = pose.residue(res2_sheet_i).atom("CA").xyz().distance(pose.residue(res1_sheet_j).atom("CA").xyz());
	arr_dis_inter_sheet[3] = pose.residue(res2_sheet_i).atom("CA").xyz().distance(pose.residue(res2_sheet_j).atom("CA").xyz());

	Real	shortest_dis_inter_sheet = get_shortest_among_4_vals(arr_dis_inter_sheet);

	Real torsion_i_j;

	if (shortest_dis_inter_sheet == arr_dis_inter_sheet[0])
	{
		Vector const& first_xyz    ( pose.residue(res2_sheet_i).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(res1_sheet_i).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(res1_sheet_j).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(res2_sheet_j).xyz("CA") );

		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}

	else if (shortest_dis_inter_sheet == arr_dis_inter_sheet[1])
	{
		Vector const& first_xyz    ( pose.residue(res2_sheet_i).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(res1_sheet_i).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(res2_sheet_j).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(res1_sheet_j).xyz("CA") );

		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}

	else if (shortest_dis_inter_sheet == arr_dis_inter_sheet[2])
	{
		Vector const& first_xyz    ( pose.residue(res1_sheet_i).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(res2_sheet_i).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(res1_sheet_j).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(res2_sheet_j).xyz("CA") );

		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}

	else // (shortest_dis_inter_sheet == arr_dis_inter_sheet[3])
	{
		Vector const& first_xyz    ( pose.residue(res1_sheet_i).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(res2_sheet_i).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(res2_sheet_j).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(res1_sheet_j).xyz("CA") );

		torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
	}
	return torsion_i_j;
}// SandwichFeatures::calculate_dihedral_w_4_resnums


Size
SandwichFeatures::write_to_sheet (
	StructureID struct_id,
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
	StructureID struct_id,
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

	

	

//get_cen_res_in_other_sheet
vector<Size>
SandwichFeatures::get_cen_res_in_other_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Size sheet_id)
{
	// <begin> get other sheet_id in same sw_can_by_sh_id
		string select_string =
		"SELECT\n"
		"	sheet_id \n"
		"FROM\n"
		"	sw_can_by_sh \n"
		"WHERE\n"
		"	(struct_id = ?) \n"
		"	AND (sw_can_by_sh_id	=	?) \n"
		"	AND (sheet_id != ?) ;";

		statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
		select_statement.bind(1,struct_id);
		select_statement.bind(2,sw_can_by_sh_id);
		select_statement.bind(3,sheet_id);
		result res(basic::database::safely_read_from_database(select_statement));

		vector<Size> other_sheet_id;
		other_sheet_id.clear();	// Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		while(res.next())
		{
			Size sheet_id;
			res >> sheet_id;
			other_sheet_id.push_back(sheet_id);
		}
	// <end> get other sheet_id in same sw_can_by_sh_id


	// <begin> get residue_begin
		select_string =
		"SELECT\n"
		"	residue_begin \n"
		"FROM\n"
		"	sheet AS sh, \n"
		"	secondary_structure_segments AS sss \n"
		"WHERE\n"
		"	(sh.struct_id = ?) \n"
		"	AND (sh.struct_id	=	sss.struct_id) \n"
		"	AND (sh.segment_id	=	sss.segment_id) \n"
		"	AND (sh.sheet_id = ?) ;";

		statement select_statement_1(basic::database::safely_prepare_statement(select_string,db_session));
		select_statement_1.bind(1,struct_id);
		select_statement_1.bind(2,other_sheet_id[other_sheet_id.size()-1]);
		result res_begin(basic::database::safely_read_from_database(select_statement_1));

		vector<Size> vector_of_residue_begin;
		vector_of_residue_begin.clear();	// Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		while(res_begin.next())
		{
			Size residue_begin;
			res_begin >> residue_begin;
			vector_of_residue_begin.push_back(residue_begin);
		}
	// <end> get residue_begin


	// <begin> get residue_end
		string select_string_2 =
		"SELECT\n"
		"	residue_end \n"
		"FROM\n"
		"	sheet AS sh, \n"
		"	secondary_structure_segments AS sss \n"
		"WHERE\n"
		"	(sh.struct_id = ?) \n"
		"	AND (sh.struct_id	=	sss.struct_id) \n"
		"	AND (sh.segment_id	=	sss.segment_id) \n"
		"	AND (sh.sheet_id = ?) ;";

		statement select_statement_2(basic::database::safely_prepare_statement(select_string_2,	db_session));
		select_statement_2.bind(1,struct_id);
		select_statement_2.bind(2,other_sheet_id[other_sheet_id.size()-1]);
		result result_end(basic::database::safely_read_from_database(select_statement_2));

		vector<Size> vector_of_residue_end;
		vector_of_residue_end.clear();	//	Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		while(result_end.next())
		{
			Size residue_end;
			result_end >> residue_end;
			vector_of_residue_end.push_back(residue_end);
		}
	// <end> get residue_end


	// <begin> get central residues
	vector<Size> vector_of_cen_residues;
	for(Size i=0; i<vector_of_residue_begin.size(); i++ )
	{
		Real to_be_rounded_i = (vector_of_residue_begin[i] + vector_of_residue_end[i])/(2.0);
		Size cen_resnum_i = round_to_Size(to_be_rounded_i);
		
		vector_of_cen_residues.push_back(cen_resnum_i);
	}
	return vector_of_cen_residues;
	// <end> get central residues
} //get_cen_res_in_other_sheet


//get_cen_residues_in_this_sheet
vector<Size>
SandwichFeatures::get_cen_residues_in_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
//		TR << "get_cen_residues_in_this_sheet" << endl;
//	time_t start_time = time(NULL);
	string select_string =
	"SELECT\n"
	"	residue_begin \n"
	"FROM\n"
	"	sheet AS sh, \n"
	"	secondary_structure_segments AS sss \n"
	"WHERE\n"
	"	(sh.struct_id = ?) \n"
	"	AND (sh.struct_id	=	sss.struct_id) \n"
	"	AND (sh.segment_id	=	sss.segment_id) \n"
	"	AND (sh.sheet_id = ?) ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sheet_id);
	result res(basic::database::safely_read_from_database(select_statement));

	vector<Size> vector_of_residue_begin;
	vector_of_residue_begin.clear();	// Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
	while(res.next())
	{
		Size residue_begin;
		res >> residue_begin;
		vector_of_residue_begin.push_back(residue_begin);
	}

	string select_string_2 =
	"SELECT\n"
	"	residue_end \n"
	"FROM\n"
	"	sheet AS sh, \n"
	"	secondary_structure_segments AS sss \n"
	"WHERE\n"
	"	(sh.struct_id = ?) \n"
	"	AND (sh.struct_id	=	sss.struct_id) \n"
	"	AND (sh.segment_id	=	sss.segment_id) \n"
	"	AND (sh.sheet_id = ?) ;";

	statement select_statement_2(basic::database::safely_prepare_statement(select_string_2,	db_session));
	select_statement_2.bind(1,struct_id);
	select_statement_2.bind(2,sheet_id);
	result result_end(basic::database::safely_read_from_database(select_statement_2));

	vector<Size> vector_of_residue_end;
	vector_of_residue_end.clear();	//	Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
	while(result_end.next())
	{
		Size residue_end;
		result_end >> residue_end;
		vector_of_residue_end.push_back(residue_end);
	}

	vector<Size> vector_of_cen_residues;
	for(Size i=0; i<vector_of_residue_begin.size(); i++ )
	{
		Real to_be_rounded_i = (vector_of_residue_begin[i] + vector_of_residue_end[i])/(2.0);
		Size cen_resnum_i = round_to_Size(to_be_rounded_i);
		
		vector_of_cen_residues.push_back(cen_resnum_i);
	}
//	time_t end_time = time(NULL);
//		TR.Info << "Finished in " << (end_time - start_time) << " seconds." << endl;
	return vector_of_cen_residues;
} //get_cen_residues_in_this_sheet

//count_AA
vector<Size>
SandwichFeatures::count_AA(
	Pose const & pose,
	Size residue_begin,
	Size residue_end)
{	
	// count AA without direction
	Size arr[] = {0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0};
	vector<Size> AA_wo_direction (arr, arr+sizeof(arr)/sizeof(arr[0]));
	for (Size ii = residue_begin; ii <= residue_end; ii++ )
	{
		if (pose.residue_type(ii).name3() == "ARG")		{			AA_wo_direction[0] = AA_wo_direction[0] + 1;		}
		else if (pose.residue_type(ii).name3() == "HIS")		{	AA_wo_direction[1] = AA_wo_direction[1] + 1;		}
		else if (pose.residue_type(ii).name3() == "LYS")		{	AA_wo_direction[2] = AA_wo_direction[2] + 1;		}
		else if (pose.residue_type(ii).name3() == "ASP")		{	AA_wo_direction[3] = AA_wo_direction[3] + 1;		}
		else if (pose.residue_type(ii).name3() == "GLU")		{	AA_wo_direction[4] = AA_wo_direction[4] + 1;		}
		else if (pose.residue_type(ii).name3() == "SER")		{	AA_wo_direction[5] = AA_wo_direction[5] + 1;		}
		else if (pose.residue_type(ii).name3() == "THR")		{	AA_wo_direction[6] = AA_wo_direction[6] + 1;		}
		else if (pose.residue_type(ii).name3() == "ASN")		{	AA_wo_direction[7] = AA_wo_direction[7] + 1;		}
		else if (pose.residue_type(ii).name3() == "GLN")		{	AA_wo_direction[8] = AA_wo_direction[8] + 1;		}
		else if (pose.residue_type(ii).name3() == "CYS")		{	AA_wo_direction[9] = AA_wo_direction[9] + 1;		}
		else if (pose.residue_type(ii).name3() == "GLY")		{	AA_wo_direction[10] = AA_wo_direction[10] + 1;		}
		else if (pose.residue_type(ii).name3() == "PRO")		{	AA_wo_direction[11] = AA_wo_direction[11] + 1;		}
		else if (pose.residue_type(ii).name3() == "ALA")		{	AA_wo_direction[12] = AA_wo_direction[12] + 1;		}
		else if (pose.residue_type(ii).name3() == "VAL")		{	AA_wo_direction[13] = AA_wo_direction[13] + 1;		}
		else if (pose.residue_type(ii).name3() == "ILE")		{	AA_wo_direction[14] = AA_wo_direction[14] + 1;		}
		else if (pose.residue_type(ii).name3() == "LEU")		{	AA_wo_direction[15] = AA_wo_direction[15] + 1;		}
		else if (pose.residue_type(ii).name3() == "MET")		{	AA_wo_direction[16] = AA_wo_direction[16] + 1;		}
		else if (pose.residue_type(ii).name3() == "PHE")		{	AA_wo_direction[17] = AA_wo_direction[17] + 1;		}
		else if (pose.residue_type(ii).name3() == "TYR")		{	AA_wo_direction[18] = AA_wo_direction[18] + 1;		}
		else if (pose.residue_type(ii).name3() == "TRP")		{	AA_wo_direction[19] = AA_wo_direction[19] + 1;		}
	}
	return AA_wo_direction;
} //count_AA

//count_AA_w_direction
vector<Size>
SandwichFeatures::count_AA_w_direction(
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	Pose const & pose_w_center_000,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end)
{
	//	TR << "count_AA_w_direction" << endl;
		
//	time_t start_time = time(NULL);

	Size arr[] = {0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0,	0,0,0,0,0};
	vector<Size> AA_w_direction (arr, arr+sizeof(arr)/sizeof(arr[0]));

//		TR << "residue_begin: " << residue_begin << endl;
//		TR << "residue_end: " << residue_end << endl;

	for (Size ii = residue_begin; ii <= residue_end; ii++ )
	{
//			TR << "resnum: " << ii << endl;

		string	heading = determine_core_heading_surface_heading_by_distance(pose_w_center_000,	ii);

		bool core_heading = "initialization_just_to_avoid_warning";

		if	(heading == "core")
		{
			core_heading = true; 			// core heading
		}
		else if	(heading == "surface")
		{
			core_heading = false;			// surface heading
		}
		else // uncertain
		{
			string	heading	=	determine_heading_direction_by_vector	(struct_id,	db_session,	pose,	sw_can_by_sh_id,	sheet_id,	residue_begin,	residue_end,	ii);
			if	(heading == "core")
			{
				core_heading = true; 			// core heading
			}
			else if	(heading == "surface")
			{
				core_heading = false;			// surface heading
			}
		}
		
		if (pose.residue_type(ii).name3() == "ARG"){
			if (core_heading)	{	AA_w_direction[0] = AA_w_direction[0] + 1;				}	//R_heading_core_num++;
			else				{	AA_w_direction[1] = AA_w_direction[1] + 1; 				}	//R_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "HIS"){
			if (core_heading)	{	AA_w_direction[2] = AA_w_direction[2] + 1;				}	//H_heading_core_num++;
			else				{	AA_w_direction[3] = AA_w_direction[3] + 1; 				}	//H_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "LYS"){
			if (core_heading)	{	AA_w_direction[4] = AA_w_direction[4] + 1;				}	//K_heading_core_num++;
			else				{	AA_w_direction[5] = AA_w_direction[5] + 1; 				}	//K_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "ASP"){
			if (core_heading)	{	AA_w_direction[6] = AA_w_direction[6] + 1;				}	//D_heading_core_num++;
			else				{	AA_w_direction[7] = AA_w_direction[7] + 1; 				}	//D_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "GLU"){
			if (core_heading)	{	AA_w_direction[8] = AA_w_direction[8] + 1;				}	//E_heading_core_num++;
			else				{	AA_w_direction[9] = AA_w_direction[9] + 1; 				}	//E_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "SER"){
			if (core_heading)	{	AA_w_direction[10] = AA_w_direction[10] + 1;				}	//S_heading_core_num++;
			else				{	AA_w_direction[11] = AA_w_direction[11] + 1; 				}	//S_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "THR"){
			if (core_heading)	{	AA_w_direction[12] = AA_w_direction[12] + 1;				}	//T_heading_core_num++;
			else				{	AA_w_direction[13] = AA_w_direction[13] + 1; 				}	//T_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "ASN"){
			if (core_heading)	{	AA_w_direction[14] = AA_w_direction[14] + 1;				}	//N_heading_core_num++;
			else				{	AA_w_direction[15] = AA_w_direction[15] + 1; 				}	//N_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "GLN"){
			if (core_heading)	{	AA_w_direction[16] = AA_w_direction[16] + 1;				}	//Q_heading_core_num++;
			else				{	AA_w_direction[17] = AA_w_direction[17] + 1; 				}	//Q_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "CYS"){
			if (core_heading)	{	AA_w_direction[18] = AA_w_direction[18] + 1;				}	//C_heading_core_num++;
			else				{	AA_w_direction[19] = AA_w_direction[19] + 1; 				}	//C_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "GLY"){
			if (core_heading)	{	AA_w_direction[20] = AA_w_direction[20] + 1;				}	//G_heading_core_num++;
			else				{	AA_w_direction[21] = AA_w_direction[21] + 1; 				}	//G_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "PRO"){
			if (core_heading)	{	AA_w_direction[22] = AA_w_direction[22] + 1;				}	//P_heading_core_num++;
			else				{	AA_w_direction[23] = AA_w_direction[23] + 1; 				}	//P_heading_surface_num++;
		}
		else	if (pose.residue_type(ii).name3() == "ALA"){
			if (core_heading)	{	AA_w_direction[24] = AA_w_direction[24] + 1;				}
			else				{	AA_w_direction[25] = AA_w_direction[25] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "VAL"){
			if (core_heading)	{	AA_w_direction[26] = AA_w_direction[26] + 1;				}
			else				{	AA_w_direction[27] = AA_w_direction[27] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "ILE"){
			if (core_heading)	{	AA_w_direction[28] = AA_w_direction[28] + 1;				}
			else				{	AA_w_direction[29] = AA_w_direction[29] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "LEU"){
			if (core_heading)	{	AA_w_direction[30] = AA_w_direction[30] + 1;				}
			else				{	AA_w_direction[31] = AA_w_direction[31] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "MET"){
			if (core_heading)	{	AA_w_direction[32] = AA_w_direction[32] + 1;				}
			else				{	AA_w_direction[33] = AA_w_direction[33] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "PHE"){
			if (core_heading)	{	AA_w_direction[34] = AA_w_direction[34] + 1;				}
			else				{	AA_w_direction[35] = AA_w_direction[35] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "TYR"){
			if (core_heading)	{	AA_w_direction[36] = AA_w_direction[36] + 1;				}
			else				{	AA_w_direction[37] = AA_w_direction[37] + 1; 				}
		}
		else	if (pose.residue_type(ii).name3() == "TRP"){
			if (core_heading)	{	AA_w_direction[38] = AA_w_direction[38] + 1;				}
			else				{	AA_w_direction[39] = AA_w_direction[39] + 1; 				}
		}
	}
//	time_t end_time = time(NULL);
//		TR.Info << "Finished in " << (end_time - start_time) << " seconds." << endl;

	return AA_w_direction;
} //count_AA_w_direction





Size
SandwichFeatures::fill_sw_by_components	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_by_components_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size sheet_id,
	string sheet_antiparallel,
	Size sw_by_components_bs_id,
	string strand_is_at_edge,
	Size component_size,
	Size residue_begin,
	Size residue_end)
{
	string insert =
	"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, sheet_id, sheet_antiparallel, sw_by_components_bs_id, strand_edge, component_size, residue_begin, residue_end, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,	?,?,?,?,?,	?,	?,?,?,	?,?,	?,?,?,?,	?,?,?,	?,?,?,?,?,?,?,?);";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	sw_by_components_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	sw_can_by_sh_id);
	insert_stmt.bind(5,	sheet_id);
	insert_stmt.bind(6,	sheet_antiparallel);
	insert_stmt.bind(7,	sw_by_components_bs_id); //bs_id
	insert_stmt.bind(8,	strand_is_at_edge);
	insert_stmt.bind(9,	component_size);
	insert_stmt.bind(10,	residue_begin);
	insert_stmt.bind(11,	residue_end);
	vector<Size>	AA_vector = count_AA(pose,	 residue_begin,	residue_end);
	insert_stmt.bind(12,	AA_vector[0]); //R_num
	insert_stmt.bind(13,	AA_vector[1]); //H_num
	insert_stmt.bind(14,	AA_vector[2]); //K_num
	insert_stmt.bind(15,	AA_vector[3]); //D
	insert_stmt.bind(16,	AA_vector[4]); //E

	insert_stmt.bind(17,	AA_vector[5]); //S
	insert_stmt.bind(18,	AA_vector[6]); //T
	insert_stmt.bind(19,	AA_vector[7]); //N
	insert_stmt.bind(20,	AA_vector[8]); //Q
	insert_stmt.bind(21,	AA_vector[9]); //C
	insert_stmt.bind(22,	AA_vector[10]); //G
	insert_stmt.bind(23,	AA_vector[11]); //P

	insert_stmt.bind(24,	AA_vector[12]); //A
	insert_stmt.bind(25,	AA_vector[13]); //V
	insert_stmt.bind(26,	AA_vector[14]); //I
	insert_stmt.bind(27,	AA_vector[15]); //L
	insert_stmt.bind(28,	AA_vector[16]); //M
	insert_stmt.bind(29,	AA_vector[17]); //F
	insert_stmt.bind(30,	AA_vector[18]); //Y
	insert_stmt.bind(31,	AA_vector[19]); //W

	basic::database::safely_write_to_database(insert_stmt);
	return 0;
} //fill_sw_by_components


string
SandwichFeatures::determine_core_heading_surface_heading_by_distance
(
	Pose const & pose_w_center_000,
	Size	ii // residue_number
)
{
	string heading;

	// <begin> determine core_heading/surface_heading by a comparison between a distance between CA and 0,0,0 and a distance between CB and 0,0,0
	//			pose::Pose pose_w_center_000 = pose;
	//			pose_w_center_000.center();

		xyzVector< core::Real > center_point(0,0,0);

		Real distance_between_CA_and_center;
		Real distance_between_CB_and_center;

		if (pose_w_center_000.residue_type(ii).name3() == "GLY")
		{
			distance_between_CA_and_center = pose_w_center_000.residue(ii).atom("CA").xyz().distance(center_point);
			distance_between_CB_and_center = pose_w_center_000.residue(ii).atom("2HA").xyz().distance(center_point);
		}
		else if (pose_w_center_000.residue_type(ii).name3() == "ALA" || pose_w_center_000.residue_type(ii).name3() == "VAL" || pose_w_center_000.residue_type(ii).name3() == "ILE" || pose_w_center_000.residue_type(ii).name3() == "SER" || pose_w_center_000.residue_type(ii).name3() == "THR" || pose_w_center_000.residue_type(ii).name3() == "CYS")
		{
			distance_between_CA_and_center = pose_w_center_000.residue(ii).atom("CA").xyz().distance(center_point);
			distance_between_CB_and_center = pose_w_center_000.residue(ii).atom("CB").xyz().distance(center_point);
		}
		else
		{
			distance_between_CA_and_center = pose_w_center_000.residue(ii).atom("CA").xyz().distance(center_point);
			distance_between_CB_and_center = pose_w_center_000.residue(ii).atom("CG").xyz().distance(center_point);
		}
		
		//			TR << "A distance between CA and center of pose: " << distance_between_CA_and_center << endl;
		//			TR << "A distance between CG (or CB) and center of pose: " << distance_between_CB_and_center << endl;
	// <end> determine core_heading/surface_heading by a comparison between a distance between CA and 0,0,0 and a distance between CB and 0,0,0

	if (distance_between_CA_and_center - distance_between_CB_and_center > 0.9)
	{
		heading = "core"; 			// core heading
	}
	else if(distance_between_CA_and_center - distance_between_CB_and_center < -0.9)
	{
		heading = "surface";			// surface heading
	}
	else
	{
		heading = "uncertain";			// surface heading
	}
	return heading;
}


string
SandwichFeatures::determine_heading_direction_by_vector
(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end,
	Size	ii // residue_number
)
{
	string heading;

	// <begin> determine core_heading/surface_heading by a vector between CA-CB of a residue and CA of the closest residue of the other sheet
		xyzVector<Real> vector_sidechain;

		if (pose.residue_type(ii).name3() != "GLY")
		{
			vector_sidechain	=	pose.xyz(NamedAtomID("CB", ii)) - pose.xyz(NamedAtomID("CA", ii));
		}
		else
		{
			vector_sidechain	=	pose.xyz(NamedAtomID("2HA", ii)) - pose.xyz(NamedAtomID("CA", ii));
		}

		Real to_be_rounded_ii = (residue_begin + residue_end)/(2.0);
		Size cen_resnum_ii = round_to_Size(to_be_rounded_ii);

		vector<Size>	vector_of_cen_residues;
		vector_of_cen_residues.clear();	// Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		vector_of_cen_residues	=	get_cen_res_in_other_sheet(struct_id, db_session, sw_can_by_sh_id,	sheet_id);

		Real shortest_dis_between_AA_and_other_sheet = 9999;
		Size jj_w_shorest_dis =	0 ; // initial value=0 just to avoid build warning at rosetta trunk
		for (Size jj = 0;	jj	<vector_of_cen_residues.size();	jj++)
		{
			Real distance = pose.residue(cen_resnum_ii).atom("CA").xyz().distance(pose.residue(vector_of_cen_residues[jj]).atom("CA").xyz());
//						
			if (distance < shortest_dis_between_AA_and_other_sheet)
			{
				shortest_dis_between_AA_and_other_sheet = distance;
				jj_w_shorest_dis = jj;
			}
		}
//
		xyzVector<Real> vector_between_AA_and_other_sheet	=	pose.xyz(NamedAtomID("CA", vector_of_cen_residues[jj_w_shorest_dis])) - pose.xyz(NamedAtomID("CA", cen_resnum_ii));

		Real	dot_product_of_vectors = dot_product( vector_sidechain, vector_between_AA_and_other_sheet );
		Real	cosine_theta = dot_product_of_vectors / (absolute_vec(vector_sidechain))*(absolute_vec(vector_between_AA_and_other_sheet));
//
//					//		TR << "cosine_theta: " << cosine_theta << endl;
//
		if (cosine_theta > 0)
		{
			heading = "core"; 			// core heading
		}
		else
		{
			heading = "surface";			// surface heading
		}
	//// <end> determine core_heading/surface_heading by a vector between CA-CB of a residue and CA of the closest residue of the other sheet


	return heading;
} //determine_heading_direction_by_vector



string
SandwichFeatures::report_heading_directions_of_all_AA_in_a_strand	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end)
{
	string heading_directions	=	"";
	for (Size	ii	=	residue_begin;	ii	<=	residue_end; ii++)
	{
		string	heading	=	determine_heading_direction_by_vector	(struct_id,	db_session,	pose,	sw_can_by_sh_id,	sheet_id,	residue_begin,	residue_end,	ii);
		heading_directions	= heading_directions	+	"	"	+	heading;
	}

	return heading_directions;
} //SandwichFeatures::report_heading_directions_of_all_AA_in_a_strand




Size
SandwichFeatures::update_sw_by_components_by_AA_w_direction	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Pose const & pose_w_center_000,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end)
{
	string insert =
	"UPDATE sw_by_components set \n"
	"R_core_heading = ? , \n"
	"R_surface_heading	=	? , \n"
	"H_core_heading = ? , \n"
	"H_surface_heading	=	? , \n"
	"K_core_heading = ? , \n"
	"K_surface_heading	=	? , \n"
	"D_core_heading = ? , \n"
	"D_surface_heading	=	? , \n"
	"E_core_heading = ? , \n"
	"E_surface_heading	=	? , \n"
	"S_core_heading = ? , \n"
	"S_surface_heading	=	? , \n"
	"T_core_heading = ? , \n"
	"T_surface_heading	=	? , \n"
	"N_core_heading = ? , \n"
	"N_surface_heading	=	? , \n"
	"Q_core_heading = ? , \n"
	"Q_surface_heading	=	? , \n"
	"C_core_heading = ? , \n"
	"C_surface_heading	=	? , \n"
	"G_core_heading = ? , \n"
	"G_surface_heading	=	? , \n"
	"P_core_heading = ? , \n"
	"P_surface_heading	=	? , \n"
	"A_core_heading = ? , \n"
	"A_surface_heading	=	? , \n"
	"V_core_heading	= ? , \n"
	"V_surface_heading	=	? ,	\n"
	"I_core_heading	= ? , \n"
	"I_surface_heading	=	? , \n"
	"L_core_heading	= ? , \n"
	"L_surface_heading	=	? , \n"
	"M_core_heading = ? , \n"
	"M_surface_heading	=	? , \n"
	"F_core_heading	= ? , \n"
	"F_surface_heading	=	? ,	\n"
	"Y_core_heading	= ? , \n"
	"Y_surface_heading	=	? , \n"
	"W_core_heading	= ? , \n"
	"W_surface_heading	=	? \n"
	"WHERE\n"
	"	residue_begin = ?	\n"
	"	AND	struct_id = ?;";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));

	vector<Size>	AA_vector = count_AA_w_direction(struct_id,	db_session,	pose,	pose_w_center_000,	sw_can_by_sh_id,	sheet_id,	residue_begin,	residue_end);
	
	insert_stmt.bind(1,	AA_vector[0]); //	R_core_heading
	insert_stmt.bind(2,	AA_vector[1]); //	R_surface_heading
	insert_stmt.bind(3,	AA_vector[2]); //	H_core_heading
	insert_stmt.bind(4,	AA_vector[3]); //
	insert_stmt.bind(5,	AA_vector[4]); //	K_core_heading
	insert_stmt.bind(6,	AA_vector[5]); //
	insert_stmt.bind(7,	AA_vector[6]); //	D_core_heading
	insert_stmt.bind(8,	AA_vector[7]); //
	insert_stmt.bind(9,	AA_vector[8]); //	E_core_heading
	insert_stmt.bind(10,	AA_vector[9]); //

	insert_stmt.bind(11,	AA_vector[10]); //	S_core_heading
	insert_stmt.bind(12,	AA_vector[11]); //	S_surface_heading
	insert_stmt.bind(13,	AA_vector[12]); //	T_core_heading
	insert_stmt.bind(14,	AA_vector[13]); //	T_surface_heading
	insert_stmt.bind(15,	AA_vector[14]); //	N_core_heading
	insert_stmt.bind(16,	AA_vector[15]); //
	insert_stmt.bind(17,	AA_vector[16]); //	Q_core_heading
	insert_stmt.bind(18,	AA_vector[17]); //

	insert_stmt.bind(19,	AA_vector[18]); //	C_core_heading
	insert_stmt.bind(20,	AA_vector[19]); //	C_surface_heading
	insert_stmt.bind(21,	AA_vector[20]); //	G_core_heading
	insert_stmt.bind(22,	AA_vector[21]); //
	insert_stmt.bind(23,	AA_vector[22]); //	P_core_heading
	insert_stmt.bind(24,	AA_vector[23]); //

	insert_stmt.bind(25,	AA_vector[24]); //	A_core_heading
	insert_stmt.bind(26,	AA_vector[25]); //	A_surface_heading
	insert_stmt.bind(27,	AA_vector[26]); //	V_core_heading
	insert_stmt.bind(28,	AA_vector[27]); //	V_surface_heading
	insert_stmt.bind(29,	AA_vector[28]); //	I_core_heading
	insert_stmt.bind(30,	AA_vector[29]); //
	insert_stmt.bind(31,	AA_vector[30]); //	L_core_heading
	insert_stmt.bind(32,	AA_vector[31]); //

	insert_stmt.bind(33,	AA_vector[32]); //	M_core_heading
	insert_stmt.bind(34,	AA_vector[33]); //	M_surface_heading
	insert_stmt.bind(35,	AA_vector[34]); //	F_core_heading
	insert_stmt.bind(36,	AA_vector[35]); //	F_surface_heading
	insert_stmt.bind(37,	AA_vector[36]); //	Y_core_heading
	insert_stmt.bind(38,	AA_vector[37]); //
	insert_stmt.bind(39,	AA_vector[38]); //	W_core_heading
	insert_stmt.bind(40,	AA_vector[39]); //

	insert_stmt.bind(41,	residue_begin);
	insert_stmt.bind(42,	struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	return 0;
} //SandwichFeatures::update_sw_by_components_by_AA_w_direction


std::pair<Size, Size>
SandwichFeatures::get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end)
{
	// <begin> retrieve current sw_by_components_bs_id
	string select_string =
	"SELECT\n"
	"	sw_by_components_bs_id \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (residue_begin = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	select_statement.bind(3,	residue_begin);
	result res(basic::database::safely_read_from_database(select_statement));

	Size current_sw_by_components_bs_id;
	while(res.next())
	{
		res >> current_sw_by_components_bs_id;
	}
	// <end> retrieve current sw_by_components_bs_id


	// <begin> retrieve other edge_strands in different sheet
	string select_string_2 =
	"SELECT\n"
	"	residue_begin, \n"
	"	residue_end \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (sheet_id != ?) \n"
	"	AND (strand_edge = \'edge\') \n"
	"	AND (residue_begin != ?);";

	statement select_statement_2(basic::database::safely_prepare_statement(select_string_2,db_session));
	select_statement_2.bind(1,	struct_id);
	select_statement_2.bind(2,	sw_can_by_sh_id);
	select_statement_2.bind(3,	sheet_id);
	select_statement_2.bind(4,	residue_begin);
	result res_2(basic::database::safely_read_from_database(select_statement_2));

	utility::vector1<SandwichFragment> other_edge_strands;
	while(res_2.next())
	{
		Size residue_begin,   residue_end;
		res_2 >> residue_begin >> residue_end;
		other_edge_strands.push_back(SandwichFragment(residue_begin, residue_end));
	}
	// <end> retrieve other edge_strands in different sheet


	// <begin> see which other edge_strand in different sheet is closest to a current strand
	SandwichFragment temp_strand_i(residue_begin, residue_end);
	Real temp_shortest = 999;
	Size residue_begin_of_nearest_strand = 0; // just initial value to avoid warning
	for(Size i=1; i<=other_edge_strands.size(); i++)
	{
		SandwichFragment temp_strand_j(other_edge_strands[i].get_start(), other_edge_strands[i].get_end());
		Real inter_strand_avg_dis = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
		if (temp_shortest > inter_strand_avg_dis)
		{
			temp_shortest = inter_strand_avg_dis;
			residue_begin_of_nearest_strand = other_edge_strands[i].get_start();
		}
	}
		// <begin> retrieve the closest sw_by_components_bs_id
		string select_string_3 =
		"SELECT\n"
		"	sw_by_components_bs_id \n"
		"FROM\n"
		"	sw_by_components \n"
		"WHERE\n"
		"	(struct_id = ?) \n"
		"	AND (sw_can_by_sh_id = ?) \n"
		"	AND (residue_begin = ?);";

		statement select_statement_3(basic::database::safely_prepare_statement(select_string_3,db_session));
		select_statement_3.bind(1,	struct_id);
		select_statement_3.bind(2,	sw_can_by_sh_id);
		select_statement_3.bind(3,	residue_begin_of_nearest_strand);
		result res_3(basic::database::safely_read_from_database(select_statement_3));

		Size closest_sw_by_components_bs_id;
		while(res_3.next())
		{
			res_3 >> closest_sw_by_components_bs_id;
		}
		// <end> retrieve the closest sw_by_components_bs_id

	// <end> see which other edge_strand in different sheet is closest to a current strand

	return std::make_pair(current_sw_by_components_bs_id, closest_sw_by_components_bs_id);
	
} // SandwichFeatures::get_current_bs_id_and_closest_edge_bs_id_in_different_sheet


Size
SandwichFeatures::report_number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size current_bs_id,
	Size closest_bs_id)
{
	// <begin> sum numbers of inward-pointing-AAs in current_bs_id and closest_bs_id
	string select_string =
	"SELECT\n"
	"	sum(R_core_heading + H_core_heading + K_core_heading + D_core_heading + E_core_heading ) \n"
	"FROM\n"
	"	sw_by_components\n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND ((sw_by_components_bs_id = ?) or (sw_by_components_bs_id = ?));";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	select_statement.bind(3,current_bs_id);
	select_statement.bind(4,closest_bs_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands;
	while(res.next())
	{
		res >> number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands ;
	}
	// <end> sum numbers of inward-pointing-AAs in current_bs_id and closest_bs_id

	// <begin> UPDATE sw_by_components table
	string insert =
	"UPDATE sw_by_components set \n"
	"number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands = ? \n"
	"WHERE\n"
	"	(sw_by_components_bs_id = ?)	\n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND	(struct_id = ?) ;";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));

	insert_stmt.bind(1,	number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands);
	insert_stmt.bind(2,	current_bs_id);
	insert_stmt.bind(3,	sw_can_by_sh_id);
	insert_stmt.bind(4,	struct_id);
	
	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sw_by_components table

	return 0;

} //	SandwichFeatures::report_number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands






Size
SandwichFeatures::report_dihedral_angle_between_core_strands_across_facing_sheets	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id)
{
	//// <begin> report average dihedral angles between core strands between facing sheets
	utility::vector1<Size> vec_sheet_id =  get_vec_distinct_sheet_id(struct_id, db_session,	sw_can_by_sh_id);

	utility::vector1<SandwichFragment> all_strands_in_sheet_i	=	get_all_strands_in_sheet_i(struct_id,	db_session,	vec_sheet_id[1]);
	utility::vector1<SandwichFragment> all_strands_in_sheet_j	=	get_all_strands_in_sheet_i(struct_id,	db_session,	vec_sheet_id[2]);
	Real total_dihedral_angle_between_core_strands_across_facing_sheets = 0;
	Size count_dihedral_angle_between_core_strands_across_facing_sheets = 0;

	for(Size i=1; i<=all_strands_in_sheet_i.size(); i++)
	{
		string edge = is_this_strand_at_edge_by_looking_db (struct_id,	db_session,	all_strands_in_sheet_i[i].get_start());
		if (edge == "core" && all_strands_in_sheet_i[i].get_size() > 3)
		{
			for(Size j=1; j<=all_strands_in_sheet_j.size(); j++)
			{
				edge = is_this_strand_at_edge_by_looking_db (struct_id,	db_session,	all_strands_in_sheet_j[j].get_start());
				if (edge == "core" && all_strands_in_sheet_j[j].get_size() > 3)
				{
					Real dihedral = calculate_dihedral_w_4_resnums(pose, all_strands_in_sheet_i[i].get_start(),	all_strands_in_sheet_i[i].get_end(),		all_strands_in_sheet_j[j].get_start(),		all_strands_in_sheet_j[j].get_end());
//						TR << "all_strands_in_sheet_i[i].get_start(): " << all_strands_in_sheet_i[i].get_start() << endl;
//						TR << "all_strands_in_sheet_i[i].get_end(): " << all_strands_in_sheet_i[i].get_end() << endl;
//						TR << "all_strands_in_sheet_j[j].get_start(): " << all_strands_in_sheet_j[j].get_start() << endl;
//						TR << "all_strands_in_sheet_j[j].get_end(): " << all_strands_in_sheet_j[j].get_end() << endl;
//						TR << "dihedral: " << dihedral << endl;
					total_dihedral_angle_between_core_strands_across_facing_sheets = total_dihedral_angle_between_core_strands_across_facing_sheets + dihedral;
					count_dihedral_angle_between_core_strands_across_facing_sheets++;
				}
			}
		}
	}
	if (count_dihedral_angle_between_core_strands_across_facing_sheets	== 0)
	{
		return 0;	// this sw is not representative to show avg dihedral angle
	}

	Real avg_dihedral_angle_between_core_strands_across_facing_sheets = total_dihedral_angle_between_core_strands_across_facing_sheets	/ count_dihedral_angle_between_core_strands_across_facing_sheets;
	Real rounded_dihedral = round_to_Real(avg_dihedral_angle_between_core_strands_across_facing_sheets);
	//// <end> report average dihedral angles between core strands between facing sheets

			
	// <begin> UPDATE sw_by_components table
	string insert =
	"UPDATE sw_by_components set \n"
	"avg_dihedral_angle_between_core_strands_across_facing_sheets = ? \n"
	"WHERE\n"
	"	(sw_can_by_sh_id = ?) \n"
	"	AND	(struct_id = ?) ;";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));

	insert_stmt.bind(1,	rounded_dihedral);
	insert_stmt.bind(2,	sw_can_by_sh_id);
	insert_stmt.bind(3,	struct_id);
	
	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sw_by_components table
	
	return 0;
} //	SandwichFeatures::report_dihedral_angle_between_core_strands_across_facing_sheets




Size
SandwichFeatures::report_avg_b_factor_CB_at_each_component	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id)
{
	//// <begin> retrieve residue_begin, residue_end at each component
	string select_string =
	"SELECT\n"
	"	residue_begin, residue_end\n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND	(sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> vector_of_residue_begin;
	utility::vector1<Size> vector_of_residue_end;

	while(res.next())
	{
		Size residue_begin,	residue_end;
		res >> residue_begin	>>	residue_end;
		vector_of_residue_begin.push_back(residue_begin);
		vector_of_residue_end.push_back(residue_end);
	}
	//// <end> retrieve residue_begin, residue_end at each component


	pose::PDBInfoCOP info = pose.pdb_info();
	if ( info )
	{
		for(Size i=1; i<=vector_of_residue_begin.size(); i++)
		{
			Real	sum_of_b_factor_CB_at_each_component	=	0;
			Size	count_atoms	=	0;
			for(Size resid=vector_of_residue_begin[i];	resid<=vector_of_residue_end[i];	resid++)
			{
				Real B_factor_of_CB = info->temperature( resid, 5 ); // '5' atom will be 'H' for Gly
				sum_of_b_factor_CB_at_each_component	=	sum_of_b_factor_CB_at_each_component	+	B_factor_of_CB;
				count_atoms++;
				//for ( Size ii = 1; ii <= info->natoms( resid ); ++ii )
				//{
				//	std::cout << "Temperature on " << resid << " " << ii << " " << info->temperature( resid, ii ) << std::endl;
				//}
			}
			Real	avg_b_factor_CB_at_each_component	=	sum_of_b_factor_CB_at_each_component/count_atoms;

			// <begin> UPDATE sw_by_components table
			string insert =
			"UPDATE sw_by_components set \n"
			"avg_b_factor_CB_at_each_component = ? \n"
			"WHERE\n"
			"	(struct_id = ?) \n"
			"	AND	(sw_can_by_sh_id = ?) \n"
			"	AND	(residue_begin = ?) ;";

			statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));

			Real	rounded_avg_b_factor = round(avg_b_factor_CB_at_each_component);
			insert_stmt.bind(1,	rounded_avg_b_factor);
			insert_stmt.bind(2,	struct_id);
			insert_stmt.bind(3,	sw_can_by_sh_id);
			insert_stmt.bind(4,	vector_of_residue_begin[i]);

			basic::database::safely_write_to_database(insert_stmt);
			// <end> UPDATE sw_by_components table

		}
	}
	return 0;
} //	SandwichFeatures::report_avg_b_factor_CB_at_each_component






Size
SandwichFeatures::report_number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size current_bs_id,
	Size closest_bs_id)
{
	// <begin> sum numbers of inward-pointing-aro_AAs in current_bs_id and closest_bs_id
	string select_string =
	"SELECT\n"
	"	sum(F_core_heading + Y_core_heading + W_core_heading) \n"
	"FROM\n"
	"	sw_by_components\n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND ((sw_by_components_bs_id = ?) or (sw_by_components_bs_id = ?));";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	select_statement.bind(3,current_bs_id);
	select_statement.bind(4,closest_bs_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands;
	while(res.next())
	{
		res >> number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands ;
	}
	// <end> sum numbers of inward-pointing-AAs in current_bs_id and closest_bs_id

	// <begin> UPDATE sw_by_components table
	string insert =
	"UPDATE sw_by_components set \n"
	"number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands = ? \n"
	"WHERE\n"
	"	(sw_by_components_bs_id = ?)	\n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND	(struct_id = ?) ;";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));

	insert_stmt.bind(1,	number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands);
	insert_stmt.bind(2,	current_bs_id);
	insert_stmt.bind(3,	sw_can_by_sh_id);
	insert_stmt.bind(4,	struct_id);
	
	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sw_by_components table

	return 0;

} //	SandwichFeatures::report_number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands



Size
SandwichFeatures::report_hydrophobic_ratio_net_charge	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id)
{
	
	// <begin> sum number_of_AA
	string select_string =
	"SELECT\n"
	"	sum(A+V+I+L+M+F+Y+W), \n"
	"	sum(R+H+K+D+E+S+T+N+Q), \n"
	"	sum(C+G+P), \n"
	"	sum(R+H+K), \n"
	"	sum(D+E) \n"
	"FROM\n"
	"	sw_by_components\n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_hydrophobic_res,	number_of_hydrophilic_res,	number_of_CGP,	number_of_RK_in_sw,	number_of_DE_in_sw;
	while(res.next())
	{
		res >> number_of_hydrophobic_res >> number_of_hydrophilic_res >> number_of_CGP >> number_of_RK_in_sw >> number_of_DE_in_sw;
	}
	// <end> sum number_of_AA


	// <begin> UPDATE sw_by_components table
	string insert =
	"UPDATE sw_by_components set \n"
	"	number_of_hydrophobic_res = ?	,	\n"
	"	number_of_hydrophilic_res = ?	,	\n"
	"	number_of_CGP = ?	,	\n"
	"	ratio_hydrophobic_philic_of_sw_in_percent = ?	,	\n"
	"	number_of_RK_in_sw = ?	,	\n"
	"	number_of_DE_in_sw = ?	,	\n"
	"	net_charge_of_sw = ?	\n"
	"WHERE\n"
	"	(sw_can_by_sh_id = ?) \n"
	"	AND	(struct_id = ?) ;";
	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));

	insert_stmt.bind(1,	number_of_hydrophobic_res);
	insert_stmt.bind(2,	number_of_hydrophilic_res);
	insert_stmt.bind(3,	number_of_CGP);
	Real ratio_hydrophobic_philic_of_sw_in_percent = (number_of_hydrophobic_res*100)/(number_of_hydrophobic_res+number_of_hydrophilic_res);
	insert_stmt.bind(4,	ratio_hydrophobic_philic_of_sw_in_percent);
	insert_stmt.bind(5,	number_of_RK_in_sw);
	insert_stmt.bind(6,	number_of_DE_in_sw);
	int net_charge_int = number_of_RK_in_sw - number_of_DE_in_sw; // for unknown reason, Size net_charge may return like '18446744073709551612', so I don't use Size here
	//	TR << "net_charge_int: " << net_charge_int << endl;
	insert_stmt.bind(7,	net_charge_int); // Net charge of His at pH 7.4 is just '+0.11' according to http://www.bmolchem.wisc.edu/courses/spring503/503-sec1/503-1a.htm
	insert_stmt.bind(8,	sw_can_by_sh_id);
	insert_stmt.bind(9,	struct_id);
	
	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sw_by_components table

	return 0;

} //	SandwichFeatures::report_hydrophobic_ratio_net_charge




// get_vec_sw_can_by_sh_id
// role:	get_distinct(sw_can_by_sh_id) as a vector
// result:	mostly get_distinct(sw_can_by_sh_id) is just 1
utility::vector1<Size>
SandwichFeatures::get_vec_sw_can_by_sh_id(
	StructureID struct_id,
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
} //get_vec_sw_can_by_sh_id



//get_size_sw_by_components_PK_id
Size
SandwichFeatures::get_size_sw_by_components_PK_id(
	StructureID struct_id,
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
	StructureID struct_id,
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
	StructureID struct_id,
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

	bool starting_res_for_connecting_strands_retrieved = false;

	Size starting_res_for_connecting_strands;
	while(res.next())
	{
		starting_res_for_connecting_strands_retrieved = true;
		res >> starting_res_for_connecting_strands;
	}

	if (!starting_res_for_connecting_strands_retrieved)
	{
		return std::make_pair(0, 0);
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
	while(res_sh_id.next())
	{
		res_sh_id >> sheet_id;
	}
	return std::make_pair(starting_res_for_connecting_strands, sheet_id);
} //	SandwichFeatures::get_starting_res_for_connecting_strands


//get_next_starting_res_for_connecting_strands
std::pair<Size, Size> //Size
SandwichFeatures::get_next_starting_res_for_connecting_strands(
	StructureID struct_id,
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
	while(res.next())
	{
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
	while(res_sh_id.next())
	{
		res_sh_id >> sh_id_of_next_start_res;
	}
	return std::make_pair(next_starting_res_for_connecting_strands, sh_id_of_next_start_res);
} //get_next_starting_res_for_connecting_strands



//update intra/inter_sheet_connection part
Size
SandwichFeatures::update_sheet_connectivity(
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	Size sw_by_components_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	string loop_kind, 
	Size intra_sheet_con_id,
	Size inter_sheet_con_id,
	string LR,
	string canonical_LR,
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
	string insert;
	Size con_id;
	Size sheet_id = 0;

	if (loop_kind == "inter_sheet") // this loop connects by a inter_sheet way
	{
		insert =
		"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, sheet_id,	LR, canonical_LR, PA_by_preceding_E, PA_by_following_E,	cano_PA,	heading_direction, parallel_EE, cano_parallel_EE,	component_size,	residue_begin, residue_end, loop_kind, inter_sheet_con_id, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,	?,?,?,	?,?,	?,?,?,?,	?,?,?,	?,?,?,?,?,?,?,?);";
		con_id = inter_sheet_con_id;
	}

	else // // this loop connects by a intra_sheet way
	{

		sheet_id =	identify_sheet_id_by_residue_end(struct_id,
														db_session,
														start_res-1); // residue_end of preceding strand

		insert =
		"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, sheet_id,	LR, canonical_LR, PA_by_preceding_E, PA_by_following_E,	cano_PA,	heading_direction, parallel_EE, cano_parallel_EE,	component_size,	residue_begin, residue_end, loop_kind, intra_sheet_con_id, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,	?,?,?,	?,?,	?,?,?,?,	?,?,?,	?,?,?,?,?,?,?,?);";
		con_id = intra_sheet_con_id;
	}

	
	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	sw_by_components_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	sw_can_by_sh_id);
	insert_stmt.bind(5,	sheet_id);
	insert_stmt.bind(6,	LR);
	insert_stmt.bind(7,	canonical_LR);
	insert_stmt.bind(8,	PA_by_preceding_E);
	insert_stmt.bind(9,	PA_by_following_E);
	insert_stmt.bind(10,	cano_PA);
	insert_stmt.bind(11,	heading_direction);
	insert_stmt.bind(12,	parallel_EE);
	insert_stmt.bind(13,	cano_parallel_EE);
	insert_stmt.bind(14,	loop_size);
	insert_stmt.bind(15,	start_res);
	insert_stmt.bind(16,	end_res);
	insert_stmt.bind(17,	loop_kind);
	insert_stmt.bind(18,	con_id);

	vector<Size>	AA_vector = count_AA(pose, start_res,	end_res);
	insert_stmt.bind(19,	AA_vector[0]); //R_num
	insert_stmt.bind(20,	AA_vector[1]); //H_num
	insert_stmt.bind(21,	AA_vector[2]); //K_num

	insert_stmt.bind(22,	AA_vector[3]); //D
	insert_stmt.bind(23,	AA_vector[4]); //E

	insert_stmt.bind(24,	AA_vector[5]); //S
	insert_stmt.bind(25,	AA_vector[6]); //T
	insert_stmt.bind(26,	AA_vector[7]); //N
	insert_stmt.bind(27,	AA_vector[8]); //Q

	insert_stmt.bind(28,	AA_vector[9]); //C
	insert_stmt.bind(29,	AA_vector[10]); //G
	insert_stmt.bind(30,	AA_vector[11]); //P

	insert_stmt.bind(31,	AA_vector[12]); //A
	insert_stmt.bind(32,	AA_vector[13]); //V
	insert_stmt.bind(33,	AA_vector[14]); //I
	insert_stmt.bind(34,	AA_vector[15]); //L
	insert_stmt.bind(35,	AA_vector[16]); //M
	insert_stmt.bind(36,	AA_vector[17]); //F
	insert_stmt.bind(37,	AA_vector[18]); //Y
	insert_stmt.bind(38,	AA_vector[19]); //W

	basic::database::safely_write_to_database(insert_stmt);

	return 0;
} //update_sheet_connectivity


//delete_this_sw_can_by_sh_id_from_sw_by_comp
Size
SandwichFeatures::delete_this_sw_can_by_sh_id_from_sw_by_comp(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
		TR << "delete_this_sw_can_by_sh_id_from_sw_by_comp with sheet_id: " << sw_can_by_sh_id << endl;
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
} //delete_this_sw_can_by_sh_id_from_sw_by_comp


//get_segment_id
Size
SandwichFeatures::get_segment_id(
	StructureID struct_id,
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


utility::vector1<Size>
SandwichFeatures::get_vec_distinct_sheet_id(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	distinct sheet_id\n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND	sw_can_by_sh_id = ? \n"
	"	AND sheet_id != 0 \n" // sheet_id = 0 when a component is 'inter_sheet_loop'
	"	AND sheet_id != 99999 ;";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,	db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));
	
	utility::vector1<Size> vec_sheet_id;
	while(res.next())
	{
		Size sheet_id;
		res >> sheet_id ;
		vec_sheet_id.push_back(sheet_id);
	}
	return vec_sheet_id;
} //get_vec_distinct_sheet_id


Size
SandwichFeatures::get_num_of_distinct_sheet_id(
	StructureID struct_id,
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
} //check_helix_existence


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

	if (loop_size == 2)  // can be applied for both hairpin and inter-sheet loop
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "T_LR";}
			else	{return "F_LR";}
		}
		else {return "U_LR";}
	}
	if (loop_size == 3)
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "T_LR";}
			else	{return "F_LR";}
		}
		else {return "U_LR";}
	}
//	if (loop_size == 4)
//	{
//		if (intra_sheet)	{return "U_LR";}
//		else
//		{
//			if (LR=="L" || LR=="BL" )	{return "T_LR";}
//			else	{return "F_LR";}
//		}
//	}
	if (loop_size == 5)
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "F_LR";}
			else	{return "T_LR";}
		}
		else {return "U_LR";}
	}
//	if (loop_size == 6)
//	{
//		if (intra_sheet)
//		{
//			if (LR=="L" || LR=="BL" ) {return "F_LR";}
//			else	{return "T_LR";}
//		}
//		else	{return "U_LR";}
//	}
//	if (loop_size == 11)
//	{
//		if (intra_sheet)	{return "U_LR";}
//		else
//		{
//			if (LR=="L" || LR=="BL" )	{return "T_LR";}
//			else	{return "F_LR";}
//		}
//	}
	else // all else loop sizes
	{
		return "U_LR";
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


//check_whether_same_direction_strands_connect_two_sheets_or_a_loop
bool
SandwichFeatures::check_whether_same_direction_strands_connect_two_sheets_or_a_loop(
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	Size start_res,
	Size next_start_res)
{
	//get other terminus of start_res
	string	select_string =
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

	//get other terminus of next_start_res
	select_string =
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

	/* as of 05/03/2013, I don't check size of strands, since I think that even with very short strands (like 2 residues long), it should not be more than 90 degree apart with following strand
	Size size_of_preceding_strand = start_res - other_end_of_start_res + 1;
	Size size_of_following_strand = other_end_of_next_start_res - next_start_res + 1;

	if (size_of_preceding_strand < 3 || size_of_following_strand < 3)
	{
		return false; // use this sandwich, since it may have very short edge strand like in [1A1N] chain A
	}*/

	///////////////////////
	// in 1QAC chain A
	// where real strand_1 starts with 4MET and ends with 7SER
	// and real strand_2 starts with 10SER and ends with 13VAL

	// other_end_of_start_res is 4MET
	// start_res is 7SER
	// next_start_res is 10SER
	// angle_start_res_being_middle is 126.5 (by both pymol and SandwichFeatures)
	
	// other_end_of_next_start_res is 13VAL
	// angle_next_start_res_being_middle is 106.1 (by both pymol and SandwichFeatures)

	// torsion_between_strands is -128.5 (by both pymol and SandwichFeatures)

	///////////////////////

	//////////////// DO NOT ERASE /////////////
	/*
	// <begin> check by angle and torsion
	// angle of other terminus of start_res, other_end_of_start_res, and next_start_res
	Vector const& first_0_xyz    ( pose.residue(other_end_of_start_res).xyz("CA") );
	Vector const& middle_0_xyz   ( pose.residue(start_res).xyz("CA") );
	Vector const& third_0_xyz    ( pose.residue(next_start_res).xyz("CA") );

	Real angle_start_res_being_middle = numeric::angle_degrees(first_0_xyz, middle_0_xyz, third_0_xyz);
		//TR.Info << "angle_start_res_being_middle: " << angle_start_res_being_middle << endl;

	// angle of start_res, next_start_res and other terminus of next_start_res
	Vector const& first_1_xyz    ( pose.residue(start_res).xyz("CA") );
	Vector const& middle_1_xyz   ( pose.residue(next_start_res).xyz("CA") );
	Vector const& third_1_xyz    ( pose.residue(other_end_of_next_start_res).xyz("CA") );

	Real angle_next_start_res_being_middle = numeric::angle_degrees(first_1_xyz, middle_1_xyz, third_1_xyz);
		TR.Info << "angle_next_start_res_being_middle: " << angle_next_start_res_being_middle << endl;

	if ((angle_start_res_being_middle > max_inter_strand_angle_to_not_be_same_direction_strands_) || (angle_next_start_res_being_middle > max_inter_strand_angle_to_not_be_same_direction_strands_))
	{
		Vector const& fourth_0_xyz   ( pose.residue(other_end_of_next_start_res).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		Real torsion_between_strands = numeric::dihedral_degrees(first_0_xyz,	middle_0_xyz, third_0_xyz, fourth_0_xyz);
			TR.Info << "torsion_between_strands: " << torsion_between_strands << endl;
			TR.Info << "abs(torsion_between_strands): " << abs(torsion_between_strands) << endl;
		if (abs(torsion_between_strands) >	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_)
		{
			return true; // don't use this sandwich, sheets are connected by same direction strand so this sandwich is not our target to extract
		}
	}
	// <end> check by angle and torsion
	*/
	//////////////// DO NOT ERASE /////////////

	// secondary check for a sandwich like 1U3J
		//TR << "start_res: " << start_res << endl;
		//TR << "other_end_of_start_res: " << other_end_of_start_res << endl;
		//TR << "other_end_of_next_start_res: " << other_end_of_next_start_res << endl;
		//TR << "next_start_res: " << next_start_res << endl;
		
	xyzVector<Real> preceding_strand	=	pose.xyz(NamedAtomID("CA", start_res)) - pose.xyz(NamedAtomID("CA", other_end_of_start_res));
	xyzVector<Real> following_strand	=	pose.xyz(NamedAtomID("CA", other_end_of_next_start_res)) - pose.xyz(NamedAtomID("CA", next_start_res));

	Real	dot_product_of_strands= dot_product( preceding_strand, following_strand );
	Real	cosine_theta = dot_product_of_strands / (absolute_vec(preceding_strand))*(absolute_vec(following_strand));
		//TR << "cosine_theta: " << cosine_theta << endl;
	if (cosine_theta >=	0)
	{
		return true; // don't use this sandwich, sheets_are_connected_by_same_direction_strand so this sandwich is not our target
	}

	return false; // use this sandwich
} //check_whether_same_direction_strands_connect_two_sheets_or_a_loop



//check_whether_hairpin_connects_short_strand
bool
SandwichFeatures::check_whether_hairpin_connects_short_strand(
	StructureID struct_id,
	sessionOP db_session,
	Size start_res,
	Size next_start_res)
{
	string	select_string =
	"SELECT\n"
	"	component_size	\n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?)\n"
	"	AND (residue_begin = ?);";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,	db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	start_res);
	result res(basic::database::safely_read_from_database(select_statement));
	
	Size component_size_1;
	while(res.next())
	{
		res >> component_size_1;
	}

	if (component_size_1 < 3)
	{
		return true; // don't use this hairpin for LR/PA
	}

	string	select_string_2 =
	"SELECT\n"
	"	component_size	\n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?)\n"
	"	AND (residue_begin = ?);";
	
	statement select_statement_2(basic::database::safely_prepare_statement(select_string_2,	db_session));
	select_statement_2.bind(1,	struct_id);
	select_statement_2.bind(2,	next_start_res);
	result res_2(basic::database::safely_read_from_database(select_statement_2));
	
	Size component_size_2;
	while(res_2.next())
	{
		res_2 >> component_size_2;
	}

	if (component_size_2 < 3)
	{
		return true; // don't use this hairpin for LR/PA
	}

	return false; // use this hairpin for LR/PA
} //check_whether_hairpin_connects_short_strand




void
SandwichFeatures::add_AA_to_terminal_loops (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size	sw_by_components_PK_id_counter,
	Size	sw_can_by_sh_id,
	string tag,
	bool starting_loop,
	Size residue_begin,
	Size residue_end)
{
	string loop_kind;
	
	if (starting_loop)
	{
		loop_kind = "starting_loop";
	}
	else // ending_loop
	{
		loop_kind = "ending_loop";
	}

	string insert =	"INSERT INTO sw_by_components (struct_id, sw_by_components_PK_id, tag, sw_can_by_sh_id, loop_kind, component_size,	residue_begin, residue_end, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,?,?,?,	?,?,?,	?,?,	?,?,?,?,	?,?,?,	?,?,?,?,?,?,?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	insert_stmt.bind(1,	struct_id);
	insert_stmt.bind(2,	sw_by_components_PK_id_counter);
	insert_stmt.bind(3,	tag);
	insert_stmt.bind(4,	sw_can_by_sh_id);
	insert_stmt.bind(5,	loop_kind);
	Size loop_size = residue_end - residue_begin + 1;
	insert_stmt.bind(6,	loop_size);
	insert_stmt.bind(7,	residue_begin);
	insert_stmt.bind(8,	residue_end);

	vector<Size>	AA_vector = count_AA(dssp_pose,	residue_begin,	residue_end);
	insert_stmt.bind(9,	AA_vector[0]); //R_num
	insert_stmt.bind(10,	AA_vector[1]); //H_num
	insert_stmt.bind(11,	AA_vector[2]); //K_num
	insert_stmt.bind(12,	AA_vector[3]); //D
	insert_stmt.bind(13,	AA_vector[4]); //E

	insert_stmt.bind(14,	AA_vector[5]); //S
	insert_stmt.bind(15,	AA_vector[6]); //T
	insert_stmt.bind(16,	AA_vector[7]); //N
	insert_stmt.bind(17,	AA_vector[8]); //Q
	insert_stmt.bind(18,	AA_vector[9]); //C
	insert_stmt.bind(19,	AA_vector[10]); //G
	insert_stmt.bind(20,	AA_vector[11]); //P

	insert_stmt.bind(21,	AA_vector[12]); //A
	insert_stmt.bind(22,	AA_vector[13]); //V
	insert_stmt.bind(23,	AA_vector[14]); //I
	insert_stmt.bind(24,	AA_vector[15]); //L
	insert_stmt.bind(25,	AA_vector[16]); //M
	insert_stmt.bind(26,	AA_vector[17]); //F
	insert_stmt.bind(27,	AA_vector[18]); //Y
	insert_stmt.bind(28,	AA_vector[19]); //W
	
	basic::database::safely_write_to_database(insert_stmt);

} //SandwichFeatures::add_AA_to_terminal_loops


Size
SandwichFeatures::add_starting_loop (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size	sw_by_components_PK_id_counter,
	Size	sw_can_by_sh_id,
	string tag)
{
	string select_string =
	"SELECT\n"
	"	min(residue_begin) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	int starting_res_of_any_strand;
	while(res.next())
	{
		res >> starting_res_of_any_strand;
	}

	Size starting_res_of_starting_loop = 0; // initial value=0 just to avoid build warning at rosetta trunk
	Size ending_res_of_starting_loop = 0 ;  // initial value=0 just to avoid build warning at rosetta trunk

	bool there_is_a_starting_loop = false;

	// I used to use Size type for i, starting_res_of_any_strand and static_cast<Size> for max_starting_loop_size_, but as of 05/02/2013, simple math like 'static_cast<Size>(starting_res_of_any_strand) - static_cast<Size>(max_starting_loop_size_)' doesn't work, so I change to use int
	for (int ii = starting_res_of_any_strand-1;
		(ii >= 1) && (ii >= (starting_res_of_any_strand) - (static_cast<int>(max_starting_loop_size_) ));
		ii-- )
	{
		char res_ss( dssp_pose.secstruct( ii ) ) ;
			
		if( res_ss == 'L')
		{
			Real dis_former_latter_AA = dssp_pose.residue(ii+1).atom("CA").xyz().distance(dssp_pose.residue(ii).atom("CA").xyz());
			if (dis_former_latter_AA < 5.0)
			{
				if (ii == starting_res_of_any_strand-1)
				{
					there_is_a_starting_loop = true;
					ending_res_of_starting_loop = ii;
				}
				starting_res_of_starting_loop = ii;
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}

	if (!there_is_a_starting_loop)
	{
		return 0;
	}

	add_AA_to_terminal_loops (struct_id,	db_session,	dssp_pose,	sw_by_components_PK_id_counter,	sw_can_by_sh_id,	tag,	true,	starting_res_of_starting_loop,	ending_res_of_starting_loop);

	return 0;
} // add_starting_loop


Size
SandwichFeatures::add_ending_loop (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size	sw_by_components_PK_id_counter,
	Size	sw_can_by_sh_id,
	string tag)
{
	string select_string =
	"SELECT\n"
	"	max(residue_end) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size ending_res_of_any_strand;
	while(res.next())
	{
		res >> ending_res_of_any_strand;
	}

	Size starting_res_of_ending_loop = 0 ;  // initial value=0 just to avoid build warning at rosetta trunk
	Size ending_res_of_ending_loop = 0;  // initial value=0 just to avoid build warning at rosetta trunk

	bool there_is_an_ending_loop = false;

	//for( Size ii = (ending_res_of_any_strand+1) ; ii <= dssp_pose.total_residue() && ii <= (ending_res_of_any_strand + max_starting_loop_size_); ii++ )
	for( Size ii = static_cast<Size>(ending_res_of_any_strand+1) ; ii <= dssp_pose.total_residue() && ii <= static_cast<Size>((static_cast<Size>(ending_res_of_any_strand) + static_cast<Size>(max_starting_loop_size_))); ii++ )
	{
		char res_ss( dssp_pose.secstruct( ii ) ) ;

		if( res_ss == 'L')
		{
			Real dis_former_latter_AA = dssp_pose.residue(ii-1).atom("CA").xyz().distance(dssp_pose.residue(ii).atom("CA").xyz());
				//TR << "dis_former_latter_AA: " << dis_former_latter_AA << endl;
			if (dis_former_latter_AA < 5.0)
			{
				if (ii == ending_res_of_any_strand+1)
				{
					there_is_an_ending_loop = true;
					starting_res_of_ending_loop = ii;
				}
				ending_res_of_ending_loop = ii;
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}

	if (!there_is_an_ending_loop)
	{
		return 0;
	}

	add_AA_to_terminal_loops (struct_id,	db_session,	dssp_pose,	sw_by_components_PK_id_counter,	sw_can_by_sh_id,	tag,	false,	starting_res_of_ending_loop,	ending_res_of_ending_loop);

	return 0;
} // add_ending_loop


Size	
SandwichFeatures::add_dssp_ratio_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sw_can_by_sh_id)
{
	Size H_num = 0;
	Size E_num = 0;
	Size L_num = 0;

	for(Size ii=1; ii<=dssp_pose.total_residue(); ii++ )
	{
		char res_ss( dssp_pose.secstruct( ii ) ) ;
		if (res_ss == 'H')		{			H_num += 1;		}
		if (res_ss == 'E')		{			E_num += 1;		}
		if (res_ss == 'L')		{			L_num += 1;		}
	}

	string update =
	"UPDATE sw_by_components set \n "
	"	H_percentage = ?	,"
	"	E_percentage = ?	,"
	"	L_percentage = ?	"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	Real rounded = round_to_Real(H_num*100/dssp_pose.total_residue());
	update_statement.bind(1,	rounded);

	rounded = round_to_Real(E_num*100/dssp_pose.total_residue());
	update_statement.bind(2,	rounded);

	rounded = round_to_Real(L_num*100/dssp_pose.total_residue());
	update_statement.bind(3,	rounded);

	update_statement.bind(4,	struct_id);
	update_statement.bind(5,	sw_can_by_sh_id);
	
	basic::database::safely_write_to_database(update_statement);

	return 0;
} // add_dssp_ratio_in_sw


void
SandwichFeatures::add_number_of_inward_pointing_W_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	sum(W_core_heading) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_inward_pointing_W_in_sw;
	while(res.next())
	{
		res >> number_of_inward_pointing_W_in_sw;
	}

	string update =
	"UPDATE sw_by_components set number_of_inward_pointing_W_in_sw = ?	"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	update_statement.bind(1,	number_of_inward_pointing_W_in_sw);
	update_statement.bind(2,	struct_id);
	update_statement.bind(3,	sw_can_by_sh_id);
	
	basic::database::safely_write_to_database(update_statement);

	//return number_of_inward_pointing_W_in_sw;
} // add_number_of_inward_pointing_W_in_sw



void
SandwichFeatures::add_number_of_inward_pointing_LWY_in_core_strands_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	sum(L_core_heading),	sum(W_core_heading),	sum(Y_core_heading) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (strand_edge = 'core') \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_inward_pointing_L_in_core_strands_in_sw;
	Size number_of_inward_pointing_W_in_core_strands_in_sw;
	Size number_of_inward_pointing_Y_in_core_strands_in_sw;
	while(res.next())
	{
		res >> number_of_inward_pointing_L_in_core_strands_in_sw	>> number_of_inward_pointing_W_in_core_strands_in_sw	>> number_of_inward_pointing_Y_in_core_strands_in_sw;
	}

	string update =
	"UPDATE sw_by_components set\n"
	"	number_of_inward_pointing_L_in_core_strands_in_sw = ?	,\n"
	"	number_of_inward_pointing_W_in_core_strands_in_sw = ?	,\n"
	" 	number_of_inward_pointing_Y_in_core_strands_in_sw = ?	\n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	update_statement.bind(1,	number_of_inward_pointing_L_in_core_strands_in_sw);
	update_statement.bind(2,	number_of_inward_pointing_W_in_core_strands_in_sw);
	update_statement.bind(3,	number_of_inward_pointing_Y_in_core_strands_in_sw);
	update_statement.bind(4,	struct_id);
	update_statement.bind(5,	sw_can_by_sh_id);
	
	basic::database::safely_write_to_database(update_statement);

	//return number_of_inward_pointing_W_in_core_strands_in_sw;
} // add_number_of_inward_pointing_LWY_in_core_strands_in_sw




Size	
SandwichFeatures::add_sw_res_size (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	min(residue_begin), max(residue_end) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size starting_res_of_sw;
	Size ending_res_of_sw;
	while(res.next())
	{
		res >> starting_res_of_sw >> ending_res_of_sw;
	}
	Size sw_res_size = ending_res_of_sw - starting_res_of_sw + 1;

	string update =
	"UPDATE sw_by_components set sw_res_size = ?	"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	update_statement.bind(1,	sw_res_size);
	update_statement.bind(2,	struct_id);
	update_statement.bind(3,	sw_can_by_sh_id);
	
	basic::database::safely_write_to_database(update_statement);

	return sw_res_size;
} // add_sw_res_size


Size	
SandwichFeatures::mark_sw_which_is_not_connected_with_continuous_atoms (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	string sw_is_not_connected_with_continuous_atoms)
{
	string update =
	"UPDATE sw_by_components set multimer_is_suspected = ?	"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	update_statement.bind(1,sw_is_not_connected_with_continuous_atoms);
	update_statement.bind(2,struct_id);
	update_statement.bind(3,sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // mark_sw_which_is_not_connected_with_continuous_atoms


Size
SandwichFeatures::add_num_strands_in_each_sw // it includes even 'short_edge_strands'
	(StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	count(*) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	strand_edge is not null \n"
	"	AND sw_can_by_sh_id = ? \n"
	"	AND struct_id = ? ;";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	sw_can_by_sh_id);
	select_statement.bind(2,	struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_strands_in_each_sw;
	while(res.next())
	{
		res >> num_strands_in_each_sw;
	}

	string update =
	"UPDATE sw_by_components set num_strands_in_each_sw = ?	"
	"WHERE\n"
	"	sw_can_by_sh_id = ? \n"
	"	AND struct_id = ?;";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	update_statement.bind(1,	num_strands_in_each_sw);
	update_statement.bind(2,	sw_can_by_sh_id);
	update_statement.bind(3,	struct_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // add_num_strands_in_each_sw


Size
SandwichFeatures::add_num_edge_strands_in_each_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	count(*) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	strand_edge=\'edge\' \n"
	"	AND sw_can_by_sh_id = ? \n"
	"	AND struct_id = ? ;";
	
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	sw_can_by_sh_id);
	select_statement.bind(2,	struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_edge_strands_in_each_sw;
	while(res.next())
	{
		res >> num_edge_strands_in_each_sw;
	}

	string update =
	"UPDATE sw_by_components set num_edge_strands_in_each_sw = ?	"
	"WHERE\n"
	"	sw_can_by_sh_id = ? \n"
	"	AND struct_id = ?;";

	statement update_statement(basic::database::safely_prepare_statement(update,	db_session));

	update_statement.bind(1,	num_edge_strands_in_each_sw);
	update_statement.bind(2,	sw_can_by_sh_id);
	update_statement.bind(3,	struct_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // add_num_edge_strands_in_each_sw


bool
SandwichFeatures::check_whether_this_pdb_should_be_excluded (
	string tag)
{
	const char* args[] = {"1W8N", "1w8n", "1W8O", "1w8o"};
		// I need to exclude these since I don't come up with how to correctly extract beta-sandwich from 1W8N

	std::vector<string> to_be_excluded (args, args+4);
	for (Size	i = 0;	i < to_be_excluded.size();	i++)
	{
		Size found = tag.find(to_be_excluded[i]);
				//TR << "found: " << found << endl;
		if (found != string::npos) // referred http://www.cplusplus.com/reference/string/string/find/
		{
			return true; // this pdb should be excluded, so don't use this pdb
		}
	}
	return false;
}

Real
SandwichFeatures::cal_min_dis_between_sheets_by_cen_res (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sheet_id_1,
	Size sheet_id_2)
{
//		TR << "cal_min_dis_between_sheets_by_cen_res" << endl;
//	time_t start_time = time(NULL);
	vector<Size>	vector_of_cen_residues_in_sheet_1;
	vector_of_cen_residues_in_sheet_1.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_cen_residues_in_sheet_1	=	get_cen_residues_in_this_sheet(struct_id, db_session,	sheet_id_1);

	vector<Size>	vector_of_cen_residues_in_sheet_2;
	vector_of_cen_residues_in_sheet_2.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_cen_residues_in_sheet_2	=	get_cen_residues_in_this_sheet(struct_id, db_session,	sheet_id_2);

	Real min_dis = 9999;

	for (Size ii=0;	ii<vector_of_cen_residues_in_sheet_1.size();	ii++){
		for (Size	jj=0;	jj<vector_of_cen_residues_in_sheet_2.size();	jj++){
			Real distance = dssp_pose.residue(vector_of_cen_residues_in_sheet_1[ii]).atom("CA").xyz().distance(dssp_pose.residue(vector_of_cen_residues_in_sheet_2[jj]).atom("CA").xyz());

			if (distance < min_dis){
				min_dis = distance;
			}
		}
	}
	//	TR << "min_dis between sheet: " << sheet_id_1 << " and sheet: " << sheet_id_2 << " is " << min_dis << endl;
//	time_t end_time = time(NULL);
//		TR.Info << "Finished in " << (end_time - start_time) << " seconds." << endl;
	return min_dis;
} //cal_min_dis_between_sheets_by_cen_res


Size
SandwichFeatures::cal_num_of_sheets_that_surround_this_sheet (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	utility::vector1<Size>	all_distinct_sheet_ids,
	Size sheet_id)
{
//		TR << "cal_num_of_sheets_that_surround_this_sheet" << endl;
//	time_t start_time = time(NULL);
	Size num_of_sheets_that_surround_sheet_id = 0;
	for(Size i=1; i <= all_distinct_sheet_ids.size(); i++)
	{ // now I check all possible combinations
		if (all_distinct_sheet_ids[i] == 99999) { //all_strands[i].get_size() < min_res_in_strand_
			continue;
		}
		if (all_distinct_sheet_ids[i] == sheet_id) {
			continue;
		}
		
		Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id
			
		if (num_strands_i < min_num_strands_in_sheet_)
		{
			continue;
		}
		
		Real min_dis_between_sheets	=	cal_min_dis_between_sheets_by_cen_res(struct_id,	db_session,	dssp_pose,	sheet_id,	all_distinct_sheet_ids[i]);
		if (min_dis_between_sheets < inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_){
			num_of_sheets_that_surround_sheet_id++;
		}
	}
	TR << "num_of_sheets_that_surround sheet_id (" << sheet_id << ") within " << inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_ << " Angstrom is " << num_of_sheets_that_surround_sheet_id << endl;
//	time_t end_time = time(NULL);
//		TR.Info << "Finished in " << (end_time - start_time) << " seconds." << endl;

	return num_of_sheets_that_surround_sheet_id;

	// if (num_of_sheets_that_surround_sheet_id > 1) // this sheet is surrounded by more than 1 other sheets!
	// if (num_of_sheets_that_surround_sheet_id == 1) // this sheet is NOT surrounded by more than 1 other sheets, so we can use these sheets to extract sandwich
} //cal_num_of_sheets_that_surround_this_sheet



//check_whether_sheets_are_connected_with_near_bb_atoms
bool
SandwichFeatures::check_whether_sheets_are_connected_with_near_bb_atoms(
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sw_can_by_sh_id)
{
	// <begin> get sheet_ids of sw_can_by_sh_id
	string select_string =
	"SELECT\n"
	"	sheet_id \n"
	"FROM\n"
	"	sw_can_by_sh \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> sheet_id_vec;
	while(res.next())
	{
		Size sheet_id;
		res >> sheet_id;
		sheet_id_vec.push_back(sheet_id);
	}
	// <end> get sheet_ids of sw_can_by_sh_id

	utility::vector1<SandwichFragment> all_strands_in_sheet_i	=	get_all_strands_in_sheet_i(struct_id,	db_session,	sheet_id_vec[1]);
		// It seems that vector1 starts with index number 1 not 0 which is typical index for typical vector
	utility::vector1<SandwichFragment> all_strands_in_sheet_j	=	get_all_strands_in_sheet_i(struct_id,	db_session,	sheet_id_vec[2]);

	// <begin> get list of residues in sheet_i, sheet_j
	utility::vector1<Size> res_num_in_sheet_i = get_list_of_residues_in_sheet_i(all_strands_in_sheet_i);
	utility::vector1<Size> res_num_in_sheet_j = get_list_of_residues_in_sheet_i(all_strands_in_sheet_j);
	// <end> get list of residues in sheet_i, sheet_j

	for(Size i=1; i<=res_num_in_sheet_i.size(); i++)
	{
		for(Size j=1; j<=res_num_in_sheet_j.size(); j++)
		{
			Real distance_1 = dssp_pose.residue(res_num_in_sheet_i[i]).atom("N").xyz().distance(dssp_pose.residue(res_num_in_sheet_j[j]).atom("O").xyz());
			Real angle_1 = 0; // just initial or angle for PRO involving case
			if (dssp_pose.residue_type(res_num_in_sheet_i[i]).name3() != "PRO")
			{
				Vector const& first_res_xyz    ( dssp_pose.residue(res_num_in_sheet_i[i]).xyz("N") );
				Vector const& middle_res_xyz    ( dssp_pose.residue(res_num_in_sheet_i[i]).xyz("H") );
				Vector const& third_res_xyz    ( dssp_pose.residue(res_num_in_sheet_j[j]).xyz("O") );
				angle_1 = numeric::angle_degrees(first_res_xyz, middle_res_xyz, third_res_xyz);
			}
			Real distance_2 = dssp_pose.residue(res_num_in_sheet_i[i]).atom("O").xyz().distance(dssp_pose.residue(res_num_in_sheet_j[j]).atom("N").xyz());
			Real angle_2 = 0; // just initial or angle for PRO involving case
			if (dssp_pose.residue_type(res_num_in_sheet_j[j]).name3() != "PRO")
			{
				Vector const& first_res_xyz_2    ( dssp_pose.residue(res_num_in_sheet_i[i]).xyz("O") );
				Vector const& middle_res_xyz_2    ( dssp_pose.residue(res_num_in_sheet_j[j]).xyz("H") );
				Vector const& third_res_xyz_2    ( dssp_pose.residue(res_num_in_sheet_j[j]).xyz("N") );
				angle_2 = numeric::angle_degrees(first_res_xyz_2, middle_res_xyz_2, third_res_xyz_2);
			}
			if (
				((distance_1 < min_N_O_dis_between_two_sheets_) && (angle_1 >  min_N_H_O_angle_between_two_sheets_))
			 || ((distance_2 < min_N_O_dis_between_two_sheets_) && (angle_2 >  min_N_H_O_angle_between_two_sheets_))
				)
			{
				return true; // don't consider this sw as canonical sw since sheets_are_connected_with_near_bb_atoms (like c.128)
			}
		}
	}
	return false; // consider this sw as canonical sw since sheets_are not connected_with_near_bb_atoms
} //check_whether_sheets_are_connected_with_near_bb_atoms





//check_whether_sw_is_not_connected_with_continuous_atoms
string	
SandwichFeatures::check_whether_sw_is_not_connected_with_continuous_atoms(
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sw_can_by_sh_id)
{
	// <begin> get starting_res_num/ending_res_num
	string select_string =
	"SELECT\n"
	"	min(residue_begin), max(residue_end) \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"	AND (sw_can_by_sh_id = ? );";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	
	result res(basic::database::safely_read_from_database(select_statement));

	Size starting_res_num, ending_res_num;
	while(res.next())
	{
		res >> starting_res_num >> ending_res_num;
	}
	// <end> get starting_res_num/ending_res_num

	for(Size ii=starting_res_num; ii<ending_res_num; ii++)
	{
		Real distance = dssp_pose.residue(ii).atom("CA").xyz().distance(dssp_pose.residue(ii+1).atom("CA").xyz());
		if (distance > 5.0)
		{
			return "multimer_suspected";
				// "don't consider this sw as canonical sw since this sandwich is not connected"
				// "This sw could be multimer like 1A78 or raw pdb file lack residues like 1A21"
		}
	}
	return "monomer"; // consider this sw as canonical sw since this sandwich is connected like 1TEN
} //check_whether_sw_is_not_connected_with_continuous_atoms



//	get_list_of_residues_in_sheet_i
utility::vector1<Size>	
SandwichFeatures::get_list_of_residues_in_sheet_i(
	utility::vector1<SandwichFragment>	all_strands_in_sheet_i)
{
	utility::vector1<Size>	list_of_residues_in_sheet_i;
	for(Size i=1; i<=all_strands_in_sheet_i.size(); i++)
	{
		for(Size j=all_strands_in_sheet_i[i].get_start(); j<=all_strands_in_sheet_i[i].get_end(); j++)
		{
			list_of_residues_in_sheet_i.push_back(j);
		}
	}
	return list_of_residues_in_sheet_i;
} // get_list_of_residues_in_sheet_i



utility::vector1<Size>
SandwichFeatures::retrieve_residue_num_of_rkde(
	StructureID struct_id,
	sessionOP db_session,
	Size	sw_can_by_sh_id,
	string	dssp_code,
	string	heading_direction)
{
	if	(dssp_code	==	"all_dssp")
	{
		string select_string =
		"SELECT\n"
		"	residue_number\n"
		"FROM\n"
		"	rkde \n"
		"WHERE\n"
		"	struct_id = ? ;";
		
		statement select_statement(basic::database::safely_prepare_statement(select_string,	db_session));
		select_statement.bind(1,	struct_id);
		result res(basic::database::safely_read_from_database(select_statement));
	
		utility::vector1<Size> vector_of_residue_num_of_rkde;
		Size	residue_num_of_rkde;
		while(res.next())
		{
			res >> residue_num_of_rkde;
			vector_of_residue_num_of_rkde.push_back(residue_num_of_rkde);
		}
		return vector_of_residue_num_of_rkde;
	}
	else // dssp_code = "E"
	{
		if	(heading_direction	==	"surface")
		{
			string select_string =
			"SELECT\n"
			"	residue_number\n"
			"FROM\n"
			"	rkde_in_strands \n"
			"WHERE\n"
			"	struct_id = ? \n"
			"	AND	heading_direction = ? \n"
			"	AND	sw_can_by_sh_id = ? ;";
			
			statement select_statement(basic::database::safely_prepare_statement(select_string,	db_session));
			select_statement.bind(1,	struct_id);
			select_statement.bind(2,	heading_direction);
			select_statement.bind(3,	sw_can_by_sh_id);
			result res(basic::database::safely_read_from_database(select_statement));
		
			utility::vector1<Size> vector_of_residue_num_of_rkde;
			Size	residue_num_of_rkde;
			while(res.next())
			{
				res >> residue_num_of_rkde;
				vector_of_residue_num_of_rkde.push_back(residue_num_of_rkde);
			}
			return vector_of_residue_num_of_rkde;
		}
		else	//heading_direction	=	"all")
		{
			string select_string =
			"SELECT\n"
			"	residue_number\n"
			"FROM\n"
			"	rkde_in_strands \n"
			"WHERE\n"
			"	struct_id = ? \n"
			"	AND	sw_can_by_sh_id = ? ;";
			
			statement select_statement(basic::database::safely_prepare_statement(select_string,	db_session));
			select_statement.bind(1,	struct_id);
			select_statement.bind(2,	sw_can_by_sh_id);
			result res(basic::database::safely_read_from_database(select_statement));
		
			utility::vector1<Size> vector_of_residue_num_of_rkde;
			Size	residue_num_of_rkde;
			while(res.next())
			{
				res >> residue_num_of_rkde;
				vector_of_residue_num_of_rkde.push_back(residue_num_of_rkde);
			}
			return vector_of_residue_num_of_rkde;
		}
	}
} //retrieve_residue_num_of_rkde



utility::vector1<Size>	
SandwichFeatures::get_vector_AA_distribution_w_direction (
	StructureID struct_id,
	sessionOP db_session,
	string heading_direction, // like core_heading, surface_heading
	string strand_location // like 'edge (strand), 'core (strand)'
	)
{
	string sum_string; // just initial declaration

	if (heading_direction == "core_heading")
	{
		sum_string =
		"SELECT\n"
		"	sum(A_core_heading), sum(C_core_heading), sum(D_core_heading), sum(E_core_heading), sum(F_core_heading), \n"
		"	sum(G_core_heading), sum(H_core_heading), sum(I_core_heading), sum(K_core_heading), sum(L_core_heading), \n"
		"	sum(M_core_heading), sum(N_core_heading), sum(P_core_heading), sum(Q_core_heading), sum(R_core_heading), \n"
		"	sum(S_core_heading), sum(T_core_heading), sum(V_core_heading), sum(W_core_heading), sum(Y_core_heading) \n"
		"FROM\n"
		"	sw_by_components \n"
		"WHERE\n"
		"	strand_edge = ? \n"
		"	AND struct_id = ? ;";
	}
	else //	heading_direction == "surface_heading"
	{
		sum_string =
		"SELECT\n"
		"	sum(A_surface_heading), sum(C_surface_heading), sum(D_surface_heading), sum(E_surface_heading), sum(F_surface_heading), \n"
		"	sum(G_surface_heading), sum(H_surface_heading), sum(I_surface_heading), sum(K_surface_heading), sum(L_surface_heading), \n"
		"	sum(M_surface_heading), sum(N_surface_heading), sum(P_surface_heading), sum(Q_surface_heading), sum(R_surface_heading), \n"
		"	sum(S_surface_heading), sum(T_surface_heading), sum(V_surface_heading), sum(W_surface_heading), sum(Y_surface_heading) \n"
		"FROM\n"
		"	sw_by_components \n"
		"WHERE\n"
		"	strand_edge = ? \n"
		"	AND struct_id = ? ;";
	}
	statement sum_statement(basic::database::safely_prepare_statement(sum_string, db_session));
	sum_statement.bind(1,	strand_location);
	sum_statement.bind(2,	struct_id);
	
	result res(basic::database::safely_read_from_database(sum_statement));

	Size num_A, num_C,	num_D,	num_E, num_F, num_G, num_H, num_I, num_K, num_L, num_M, num_N, num_P, num_Q, num_R, num_S, num_T, num_V, num_W, num_Y;

	utility::vector1<Size> vector_AA_distribution_w_direction;

	while(res.next())
	{
		res >> num_A >>  num_C >> 	num_D >> 	num_E >>  num_F >>  num_G >>  num_H >>  num_I >>  num_K >>  num_L >>  num_M >>  num_N >>  num_P >>  num_Q >>  num_R >>  num_S >>  num_T >>  num_V >>  num_W >>  num_Y;
		vector_AA_distribution_w_direction.push_back(num_A);
		vector_AA_distribution_w_direction.push_back(num_C);
		vector_AA_distribution_w_direction.push_back(num_D);
		vector_AA_distribution_w_direction.push_back(num_E);
		vector_AA_distribution_w_direction.push_back(num_F);

		vector_AA_distribution_w_direction.push_back(num_G);
		vector_AA_distribution_w_direction.push_back(num_H);
		vector_AA_distribution_w_direction.push_back(num_I);
		vector_AA_distribution_w_direction.push_back(num_K);
		vector_AA_distribution_w_direction.push_back(num_L);

		vector_AA_distribution_w_direction.push_back(num_M);
		vector_AA_distribution_w_direction.push_back(num_N);
		vector_AA_distribution_w_direction.push_back(num_P);
		vector_AA_distribution_w_direction.push_back(num_Q);
		vector_AA_distribution_w_direction.push_back(num_R);

		vector_AA_distribution_w_direction.push_back(num_S);
		vector_AA_distribution_w_direction.push_back(num_T);
		vector_AA_distribution_w_direction.push_back(num_V);
		vector_AA_distribution_w_direction.push_back(num_W);
		vector_AA_distribution_w_direction.push_back(num_Y);
	}

	return vector_AA_distribution_w_direction;
} // get_vector_AA_distribution_w_direction



string
SandwichFeatures::get_residue_location (
	StructureID struct_id,
	sessionOP db_session,
	Size	residue_num
	)
{
	string	sum_string =
		"SELECT\n"
		"	strand_edge \n"
		"FROM\n"
		"	sw_by_components \n"
		"WHERE\n"
		"	? between residue_begin and residue_end \n "
		"	AND struct_id = ? ;";

	statement sum_statement(basic::database::safely_prepare_statement(sum_string, db_session));
	sum_statement.bind(1,	residue_num);
	sum_statement.bind(2,	struct_id);
	result res(basic::database::safely_read_from_database(sum_statement));

	string strand_edge;

	while(res.next())
	{
		res >> strand_edge;
	}

	string	residue_location;
	if	(strand_edge == "edge")
	{
		residue_location = "edge";
	}
	else if	(strand_edge == "core")
	{
		residue_location = "core";
	}
	else
	{
		residue_location = "loop_or_short_edge";
	}
	return residue_location;
}


utility::vector1<Size>
SandwichFeatures::get_vector_AA_distribution_wo_direction (
	StructureID struct_id,
	sessionOP db_session,
	string loop_kind // like 'hairpin' or 'inter-sheet-loop'
	)
{
	string	sum_string =
		"SELECT\n"
		"	sum(A), sum(C), sum(D), sum(E), sum(F), \n"
		"	sum(G), sum(H), sum(I), sum(K), sum(L), \n"
		"	sum(M), sum(N), sum(P), sum(Q), sum(R), \n"
		"	sum(S), sum(T), sum(V), sum(W), sum(Y) \n"
		"FROM\n"
		"	sw_by_components \n"
		"WHERE\n"
		"	loop_kind = ? \n"
		"	AND struct_id = ? ;";

	statement sum_statement(basic::database::safely_prepare_statement(sum_string, db_session));
	sum_statement.bind(1,	loop_kind);
	sum_statement.bind(2,	struct_id);
	result res(basic::database::safely_read_from_database(sum_statement));

	Size num_A, num_C,	num_D,	num_E, num_F, num_G, num_H, num_I, num_K, num_L, num_M, num_N, num_P, num_Q, num_R, num_S, num_T, num_V, num_W, num_Y;

	utility::vector1<Size> vector_AA_distribution_wo_direction;

	while(res.next())
	{
		res >> num_A >>  num_C >> 	num_D >> 	num_E >>  num_F >>  num_G >>  num_H >>  num_I >>  num_K >>  num_L >>  num_M >>  num_N >>  num_P >>  num_Q >>  num_R >>  num_S >>  num_T >>  num_V >>  num_W >>  num_Y;
		vector_AA_distribution_wo_direction.push_back(num_A);
		vector_AA_distribution_wo_direction.push_back(num_C);
		vector_AA_distribution_wo_direction.push_back(num_D);
		vector_AA_distribution_wo_direction.push_back(num_E);
		vector_AA_distribution_wo_direction.push_back(num_F);

		vector_AA_distribution_wo_direction.push_back(num_G);
		vector_AA_distribution_wo_direction.push_back(num_H);
		vector_AA_distribution_wo_direction.push_back(num_I);
		vector_AA_distribution_wo_direction.push_back(num_K);
		vector_AA_distribution_wo_direction.push_back(num_L);

		vector_AA_distribution_wo_direction.push_back(num_M);
		vector_AA_distribution_wo_direction.push_back(num_N);
		vector_AA_distribution_wo_direction.push_back(num_P);
		vector_AA_distribution_wo_direction.push_back(num_Q);
		vector_AA_distribution_wo_direction.push_back(num_R);

		vector_AA_distribution_wo_direction.push_back(num_S);
		vector_AA_distribution_wo_direction.push_back(num_T);
		vector_AA_distribution_wo_direction.push_back(num_V);
		vector_AA_distribution_wo_direction.push_back(num_W);
		vector_AA_distribution_wo_direction.push_back(num_Y);
	}

	return vector_AA_distribution_wo_direction;
} // get_vector_AA_distribution_wo_direction



utility::vector1<Size>	
SandwichFeatures::get_vec_AA_kind (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{

	// <begin> sum number_of_AA
	string select_string =
	"SELECT\n"
	"	sum(R+K), \n"	//	positive
	"	sum(D+E), \n"	//	negative
	"	sum(S+T+N+Q), \n"	//	polar
	"	sum(H+F+Y+W), \n"	//	aromatic
	"	sum(C+G+P+A+V+I+L+M) \n"	//	hydrophobic
	"FROM\n"
	"	sw_by_components\n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> vec_AA_kind;

	Size pos,	neg,	polar,	aro,	pho;
	while(res.next())
	{
		res >> pos >> neg >> polar >> aro >> pho;
	}
	// <end> sum number_of_AA

	vec_AA_kind.push_back(pos);
	vec_AA_kind.push_back(neg);
	vec_AA_kind.push_back(polar);
	vec_AA_kind.push_back(aro);
	vec_AA_kind.push_back(pho);

	return vec_AA_kind;
} // get_vec_AA_kind


// check_whether_sw_by_sh_id_still_alive
bool	
SandwichFeatures::check_whether_sw_by_sh_id_still_alive(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	sw_by_components_PK_id \n"
	"FROM\n"
	"	sw_by_components\n"
	"WHERE\n"
	"	(sw_can_by_sh_id = ?)\n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,sw_can_by_sh_id);
	select_statement.bind(2,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	bool sw_by_sh_id_still_alive = false;
	while(res.next())
	{
		sw_by_sh_id_still_alive = true;
	}
	return sw_by_sh_id_still_alive;
} //check_whether_sw_by_sh_id_still_alive



// identify_sheet_id_by_residue_end
Size
SandwichFeatures::identify_sheet_id_by_residue_end(
	StructureID struct_id,
	sessionOP db_session,
	Size residue_end)
{
	string select_string =
	"SELECT\n"
	"	sheet_id \n"
	"FROM\n"
	"	sw_by_components \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (residue_end = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,residue_end);
	result res(basic::database::safely_read_from_database(select_statement));

	Size sheet_id;
	while(res.next())
	{
		res >> sheet_id;
	}
	return sheet_id;
} //identify_sheet_id_by_residue_end



string
SandwichFeatures::report_turn_type(
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size first_res,
	Size second_res,
	StructureID struct_id,
	sessionOP db_session)
{
	Real first_phi = pose.phi(first_res);
	Real first_psi = pose.psi(first_res);

	Real second_phi = pose.phi(second_res);
	Real second_psi = pose.psi(second_res);

		// I use mean dihedral values in Protein Science (1994), 3:2207-2216 "A revised set of potentials for beta-turn formation in proteins" by Hutchinson and Thornton
		// worth to be referred	http://en.wikipedia.org/wiki/Turn_(biochemistry)#Hairpins
		// I didn't use Brian's BetaTurnDetectionFeatures since I don't understand it fully
	string turn_type = "uncertain";
	
	if (	(first_phi > (-64-allowed_deviation_for_turn_type_id_) && first_phi < (-64+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (-27-allowed_deviation_for_turn_type_id_) && first_psi < (-27+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (-90-allowed_deviation_for_turn_type_id_) && second_phi < (-90+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (-7-allowed_deviation_for_turn_type_id_) && second_psi < (-7+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "I";
	}

	else if (	(first_phi > (-60-allowed_deviation_for_turn_type_id_) && first_phi < (-60+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (131-allowed_deviation_for_turn_type_id_) && first_psi < (131+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (84-allowed_deviation_for_turn_type_id_) && second_phi < (84+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (1-allowed_deviation_for_turn_type_id_) && second_psi < (1+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "II";
	}

	else if (	(first_phi > (-72-allowed_deviation_for_turn_type_id_) && first_phi < (-72+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (-33-allowed_deviation_for_turn_type_id_) && first_psi < (-33+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (-123-allowed_deviation_for_turn_type_id_) && second_phi < (-123+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (121-allowed_deviation_for_turn_type_id_) && second_psi < (121+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "VIII";
	}

	else if (	(first_phi > (55-allowed_deviation_for_turn_type_id_) && first_phi < (55+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (38-allowed_deviation_for_turn_type_id_) && first_psi < (38+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (78-allowed_deviation_for_turn_type_id_) && second_phi < (78+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (6-allowed_deviation_for_turn_type_id_) && second_psi < (6+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "I_prime";
	}

	else if (	(first_phi > (60-allowed_deviation_for_turn_type_id_) && first_phi < (60+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (-126-allowed_deviation_for_turn_type_id_) && first_psi < (-126+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (-91-allowed_deviation_for_turn_type_id_) && second_phi < (-91+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (1-allowed_deviation_for_turn_type_id_) && second_psi < (1+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "II_prime";
	}

	else if (	(first_phi > (-64-allowed_deviation_for_turn_type_id_) && first_phi < (-64+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (142-allowed_deviation_for_turn_type_id_) && first_psi < (142+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (-93-allowed_deviation_for_turn_type_id_) && second_phi < (-93+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (5-allowed_deviation_for_turn_type_id_) && second_psi < (5+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "VIa1";
	}

	else if (	(first_phi > (-132-allowed_deviation_for_turn_type_id_) && first_phi < (-132+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (139-allowed_deviation_for_turn_type_id_) && first_psi < (139+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (-80-allowed_deviation_for_turn_type_id_) && second_phi < (-80+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (-10-allowed_deviation_for_turn_type_id_) && second_psi < (-10+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "VIa2";
	}

	else if (	(first_phi > (-135-allowed_deviation_for_turn_type_id_) && first_phi < (-135+allowed_deviation_for_turn_type_id_))
		&&	(first_psi > (131-allowed_deviation_for_turn_type_id_) && first_psi < (131+allowed_deviation_for_turn_type_id_))
		&&	(second_phi > (-76-allowed_deviation_for_turn_type_id_) && second_phi < (-76+allowed_deviation_for_turn_type_id_))
		&&	(second_psi > (157-allowed_deviation_for_turn_type_id_) && second_psi < (157+allowed_deviation_for_turn_type_id_)))
	{
		turn_type = "VIa2";
	}

	else
	{
		turn_type = "IV";
	}


	string select_string =
	"UPDATE sw_by_components set turn_type = ?	"
	"WHERE\n"
	"	(sw_can_by_sh_id = ?) \n"
	"	AND	(residue_begin = ?) \n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	turn_type);
	select_statement.bind(2,	sw_can_by_sh_id);
	select_statement.bind(3,	first_res);
	select_statement.bind(4,	struct_id);
	basic::database::safely_write_to_database(select_statement);

	return turn_type;
} //report_turn_type








void
SandwichFeatures::report_turn_AA(
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size i,
	StructureID struct_id,
	sessionOP db_session,
	string turn_type)
{
	string canonical_turn_AA	=	"F_canonical_turn_AA";
	if	(turn_type == "I")
	{
		if (pose.residue_type(i).name3() == "LEU" || pose.residue_type(i).name3() == "ALA" || pose.residue_type(i).name3() == "GLY" || pose.residue_type(i).name3() == "PRO" ||
			pose.residue_type(i).name3() == "THR" || pose.residue_type(i).name3() == "SER" || pose.residue_type(i).name3() == "GLU" || pose.residue_type(i).name3() == "ASN" ||
			pose.residue_type(i).name3() == "ASP")
		{
			if (pose.residue_type(i+1).name3() == "LEU" ||	pose.residue_type(i+1).name3() == "ALA" || pose.residue_type(i+1).name3() == "PRO" ||
				pose.residue_type(i+1).name3() == "THR" ||	pose.residue_type(i+1).name3() == "SER" || pose.residue_type(i+1).name3() == "GLU" ||
				pose.residue_type(i+1).name3() == "ASP" ||	pose.residue_type(i+1).name3() == "LYS" )
			{
				if (pose.residue_type(i+2).name3() == "ALA" ||	pose.residue_type(i+2).name3() == "GLY" || pose.residue_type(i+2).name3() == "THR" ||
					pose.residue_type(i+2).name3() == "SER" ||	pose.residue_type(i+2).name3() == "GLU" || pose.residue_type(i+2).name3() == "ASN" ||
					pose.residue_type(i+2).name3() == "ASP" ||	pose.residue_type(i+2).name3() == "LYS" )
				{
					if (pose.residue_type(i+3).name3() == "VAL" ||	pose.residue_type(i+3).name3() == "LEU"	||	pose.residue_type(i+3).name3() == "ALA" ||
						pose.residue_type(i+3).name3() == "GLY" ||	pose.residue_type(i+3).name3() == "THR"	||	pose.residue_type(i+3).name3() == "SER" ||
						pose.residue_type(i+3).name3() == "GLU" ||	pose.residue_type(i+3).name3() == "ASP"	||	pose.residue_type(i+3).name3() == "LYS" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}
	else if	(turn_type == "II")
	{
		if (pose.residue_type(i).name3() == "VAL" || pose.residue_type(i).name3() == "LEU" || pose.residue_type(i).name3() == "ALA" || pose.residue_type(i).name3() == "GLY" ||
			pose.residue_type(i).name3() == "TYR" || pose.residue_type(i).name3() == "PRO" || pose.residue_type(i).name3() == "GLU" || pose.residue_type(i).name3() == "LYS")
		{
			if (pose.residue_type(i+1).name3() == "ALA" ||	pose.residue_type(i+1).name3() == "PRO" || pose.residue_type(i+1).name3() == "SER" ||
				pose.residue_type(i+1).name3() == "GLU" ||	pose.residue_type(i+1).name3() == "LYS" )
			{
				if (pose.residue_type(i+2).name3() == "GLY")
				{
					if (pose.residue_type(i+3).name3() == "VAL" ||	pose.residue_type(i+3).name3() == "ALA"	||	pose.residue_type(i+3).name3() == "SER" ||
						pose.residue_type(i+3).name3() == "GLU" ||	pose.residue_type(i+3).name3() == "LYS" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}
	else if	(turn_type == "VIII")
	{
		if (pose.residue_type(i).name3() == "GLY" || pose.residue_type(i).name3() == "PRO")
		{
			if (pose.residue_type(i+1).name3() == "PRO" ||	pose.residue_type(i+1).name3() == "ASP" )
			{
				if (pose.residue_type(i+2).name3() == "VAL"	||	pose.residue_type(i+2).name3() == "LEU"	||	pose.residue_type(i+2).name3() == "ASN"	||	pose.residue_type(i+2).name3() == "ASP")
				{
					if (pose.residue_type(i+3).name3() == "PRO" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}
	else if	(turn_type == "I_prime")
	{
		if (pose.residue_type(i).name3() == "ILE"	||	pose.residue_type(i).name3() == "VAL"	||	pose.residue_type(i).name3() == "LEU"	||
			pose.residue_type(i).name3() == "ALA"	||	pose.residue_type(i).name3() == "TYR"	||	pose.residue_type(i).name3() == "THR"	||
			pose.residue_type(i).name3() == "SER"	||	pose.residue_type(i).name3() == "ASP"	||	pose.residue_type(i).name3() == "LYS")
		{
			if (pose.residue_type(i+1).name3() == "PRO" ||	pose.residue_type(i+1).name3() == "GLY"	 ||	pose.residue_type(i+1).name3() == "HIS" ||
			 	pose.residue_type(i+1).name3() == "ASN" ||	pose.residue_type(i+1).name3() == "ASP" )
			{
				if (pose.residue_type(i+2).name3() == "GLY")
				{
					if (pose.residue_type(i+3).name3() == "VAL"	||	pose.residue_type(i+3).name3() == "GLU" || pose.residue_type(i+3).name3() == "ASN" ||
						pose.residue_type(i+3).name3() == "LYS"	||	pose.residue_type(i+3).name3() == "ARG" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}
	else if	(turn_type == "II_prime")
	{
			TR << "pose.residue_type(i).name3(): " << pose.residue_type(i).name3() << endl;
		if (pose.residue_type(i).name3() == "PHE"	||	pose.residue_type(i).name3() == "VAL"	||	pose.residue_type(i).name3() == "LEU"	||
			pose.residue_type(i).name3() == "ALA"	||	pose.residue_type(i).name3() == "GLY"	||
			pose.residue_type(i).name3() == "TYR"	||	pose.residue_type(i).name3() == "THR"	||
			pose.residue_type(i).name3() == "SER"	||	pose.residue_type(i).name3() == "HIS"	||	pose.residue_type(i).name3() == "GLU"	||
			pose.residue_type(i).name3() == "ASN"	||
			pose.residue_type(i).name3() == "GLN"	||	pose.residue_type(i).name3() == "ASP"	||	pose.residue_type(i).name3() == "ARG")
		{
				TR << "pose.residue_type(i+1).name3(): " << pose.residue_type(i+1).name3() << endl;
			if (pose.residue_type(i+1).name3() == "GLY")
			{
						TR << "pose.residue_type(i+2).name3(): " << pose.residue_type(i+2).name3() << endl;
				if (pose.residue_type(i+2).name3() == "LEU"	||	pose.residue_type(i+2).name3() == "ALA"	||	pose.residue_type(i+2).name3() == "GLY"	||
					pose.residue_type(i+2).name3() == "PRO"	||
					pose.residue_type(i+2).name3() == "SER"	||	pose.residue_type(i+2).name3() == "GLU"	||
					pose.residue_type(i+2).name3() == "ASN"	||	pose.residue_type(i+2).name3() == "ASP"	|| pose.residue_type(i+2).name3() == "LYS")
				{
							TR << "pose.residue_type(i+3).name3(): " << pose.residue_type(i+3).name3() << endl;
					if (pose.residue_type(i+3).name3() == "PHE"	||	pose.residue_type(i+3).name3() == "VAL"	||	pose.residue_type(i+3).name3() == "LEU"	||
						pose.residue_type(i+3).name3() == "ALA"	||
						pose.residue_type(i+3).name3() == "GLY" ||	pose.residue_type(i+3).name3() == "TYR"	||	pose.residue_type(i+3).name3() == "THR"	||
						pose.residue_type(i+3).name3() == "SER"	||	pose.residue_type(i+3).name3() == "GLU"	||
						pose.residue_type(i+3).name3() == "ASN"	||	pose.residue_type(i+3).name3() == "GLN"	||	pose.residue_type(i+3).name3() == "LYS"	||	pose.residue_type(i+3).name3() == "ARG" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}

	else if	(turn_type == "VIa1" ||	turn_type == "VIa2")
	{
		if (pose.residue_type(i).name3() == "PHE"	||	pose.residue_type(i).name3() == "VAL"	||	pose.residue_type(i).name3() == "THR"	||	pose.residue_type(i).name3() == "HIS"	||
			pose.residue_type(i).name3() == "ASN")
		{
			if (pose.residue_type(i+1).name3() == "ILE"	|| pose.residue_type(i+1).name3() == "SER" || pose.residue_type(i+1).name3() == "ASN")
			{
				if (pose.residue_type(i+2).name3() == "PRO")
				{
					if (pose.residue_type(i+3).name3() == "GLY" ||	pose.residue_type(i+3).name3() == "THR" ||	pose.residue_type(i+3).name3() == "HIS" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}

	else if	(turn_type == "VIb")
	{
		if (pose.residue_type(i).name3() == "PHE"	||	pose.residue_type(i).name3() == "GLY"	||	pose.residue_type(i).name3() == "THR"	||	pose.residue_type(i).name3() == "SER")
		{
			if (pose.residue_type(i+1).name3() == "LEU"	|| pose.residue_type(i+1).name3() == "TYR" || pose.residue_type(i+1).name3() == "THR"	 || pose.residue_type(i+1).name3() == "GLU")
			{
				if (pose.residue_type(i+2).name3() == "PRO")
				{
					if (pose.residue_type(i+3).name3() == "PHE"	||	pose.residue_type(i+3).name3() == "ALA" ||	pose.residue_type(i+3).name3() == "TYR"	||
						pose.residue_type(i+3).name3() == "THR"	||	pose.residue_type(i+3).name3() == "LYS" )
					{
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	}

	else	//(turn_type == 'IV')
	{
		canonical_turn_AA = "uncertain_canonical_turn_AA_since_turn_type_eq_IV";
	}


	string select_string =
	"UPDATE sw_by_components set \n"
	"i_AA = ? , \n"
	"i_p1_AA	=	? , \n"
	"i_p2_AA	=	? , \n"
	"i_p3_AA	=	? ,	\n"
	"canonical_turn_AA	=	? \n"
	"WHERE\n"
	"	(sw_can_by_sh_id = ?) \n"
	"	AND	(residue_begin = ?) \n"
	"	AND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	pose.residue_type(i).name3());
	select_statement.bind(2,	pose.residue_type(i+1).name3());
	select_statement.bind(3,	pose.residue_type(i+2).name3());
	select_statement.bind(4,	pose.residue_type(i+3).name3());
	select_statement.bind(5,	canonical_turn_AA);
	select_statement.bind(6,	sw_can_by_sh_id);
	select_statement.bind(7,	(i+1));
	select_statement.bind(8,	struct_id);
	basic::database::safely_write_to_database(select_statement);
	
} //report_turn_AA






void
SandwichFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
)
{
	min_num_strands_to_deal_ = tag->getOption<Size>("min_num_strands_to_deal", 4);
					// At least 4 strands should be in pdb file
	max_num_strands_to_deal_ = tag->getOption<Size>("max_num_strands_to_deal", 140);
					// example: (in all chains of 1FE8) There are 132 strands, it took ~ 7 cpu minutes to process
	min_res_in_strand_ = tag->getOption<Size>("min_res_in_strand", 2);
					// definition: minimum number of residues in a strand, for edge strand definition & analysis
					// example: 4=< is recommended (in 1A8M) min_res_in_strand = 2, (in 1PMY) min_res_in_strand = 3
	min_CA_CA_dis_ = tag->getOption<Real>("min_CA_CA_dis", 3.5);
					// definition: minimum CA_CA_distance between strands in same sheet
					// example: (in 1A8M) 'min_CA_CA_dis_= 3.5', (in 1KIT) 'min_CA_CA_dis_= 4.0'

	max_CA_CA_dis_ = tag->getOption<Real>("max_CA_CA_dis", 6.2);
					// example: (in 1A8M) 'max_CA_CA_dis_= 6.2', (in 1KIT) 'max_CA_CA_dis_= 5.7'

	min_N_O_dis_between_two_sheets_ = tag->getOption<Real>("min_N_O_dis_between_sheets", 3.3);
					//	definition: min distance between bb N and bb O between two sheets
					//	example: (in c.10.0, 4.7 Angstrom exists between O and N of edge strands of two sheets) this is a canonical sw
					//	example: (in c.128.0, 3.1 Angstrom exists between O and N of edge strands of two sheets) this is not a canonical sw

	min_N_H_O_angle_between_two_sheets_ = tag->getOption<Real>("min_N_H_O_angle_between_two_sheets", 154.0);
					//	definition: minimum N-H-O angle between bb N and bb O between two sheets
					//	example: (in c.30.0, 147.2 degree angle exists between O and N of edge strands of two sheets) this is a canonical sw
					//	example: (in c.26.0, 155.5 degree angle exists between O and N of edge strands of two sheets) this is not a canonical sw
					//	example: (in c.128.0, 156.4 degree angle exists between O and N of edge strands of two sheets) this is not a canonical sw

	min_C_O_N_angle_ = tag->getOption<Real>("min_C_O_N_angle", 120.0);
					// example: (in 1L0Q chain A), 138 is the smallest C_O_N_angle (C and O from one sheet, N from other sheet)
	min_sheet_dis_ = tag->getOption<Real>("min_sheet_dis", 7.0);
					// definition: minimum CA_CA_distance between strands in different sheets to constitute a sandwich
					// 7 Angstrom seems OK
	max_sheet_dis_ = tag->getOption<Real>("max_sheet_dis", 15.0);
					// definition: maximum CA_CA_distance between strands in different sheets to constitute a sandwich
					// 15 Angstrom seems OK
	max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_ = tag->getOption<Real>("max_sheet_angle_with_cen_res", 130.0);
					//	definition: Maximum angle between sheets (CA and CA) with two terminal residues and one central residue
					//	example: I need to remove non-facing sheets
	min_sheet_angle_by_four_term_cen_res_ = tag->getOption<Real>("min_sheet_angle", 30.0);
					//	definition: Minimum angle between sheets (CA and CA)
					//	usage: used in judge_facing "angle_1 > min_sheet_angle_by_four_term_cen_res_"
	max_sheet_angle_by_four_term_cen_res_ = tag->getOption<Real>("max_sheet_angle", 150.0);
					//	definition: Maximum angle between sheets (CA and CA)
					// In [1TEN] even 155 degree comes from same sheet!

	min_sheet_torsion_cen_res_ = tag->getOption<Real>("min_sheet_torsion_cen_res", -150.0);
					//	definition: "Minimum torsion between sheets (CA and CA) with respect to terminal central residues in each beta-sheet
					//	explanation: with respect to central residues, one torsion angles of 1TEN is 84.9
	max_sheet_torsion_cen_res_ = tag->getOption<Real>("max_sheet_torsion_cen_res", 150.0);
					//	definition: "maximum torsion between sheets (CA and CA) with respect to terminal central residues in each beta-sheet
					//	usage: used in judge_facing "torsion_i_j < max_sheet_torsion_cen_res_"
	min_num_strands_in_sheet_ = tag->getOption<Size>("min_num_strands_in_sheet", 3);
					//  definition: a sheet with < 3 strands will be ignored
					//	usage: if (num_strands_i < min_num_strands_in_sheet_)
	min_inter_sheet_dis_CA_CA_ = tag->getOption<Real>("min_inter_sheet_dis_CA_CA", 4.0);
					//	example:	(in 12E8) the distance between S52 and G64 is 4.2 A
	max_inter_sheet_dis_CA_CA_ = tag->getOption<Real>("max_inter_sheet_dis_CA_CA", 24.0);
					//	example:	(in 2WYF) the distance between val and gln is 5.8 A
					//	usage:	shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_
					//	example:	(in 1TEN) shortest_avg_dis_inter_sheet between sheet 1 and 2 = 11.6 A (and these sheets should be a sandwich)
					//	example:	(in 1A64_chain_A) shortest_avg_dis_inter_sheet between sheet 1 and 2 = 25 A (and these sheets should not be a sandwich)
					//	example:	(in 1A1M) the average distance between sheet 1 and 4 > 20 A (but these sheets should be a sandwich)
					//	example:	(in 1ASQ) the average distance between sheet 1 and 4 > 20 A (and these sheets should not be a sandwich)


	max_inter_strand_angle_to_not_be_same_direction_strands_ = tag->getOption<Real>("max_inter_strand_angle_to_not_be_same_direction_strands", 120.0);
					//	usage: 	if (angle_start_res_being_middle > max_inter_strand_angle_to_not_be_same_direction_strands_)
					//	example: (in 1BQB chain A) 121 is possible (but this should be excluded as same direction strand)
					//	example: (in 1A0Q chain L) 127 is possible for 3-6-9 angle (but this should be excluded as same direction strand)
	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_ = tag->getOption<Real>("max_abs_inter_strand_dihedral_to_not_be_same_direction_strands", 100.0);
					//	usage:	if (abs(torsion_between_strands) >	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_)
					//	example: (in 1U3J chain A) 105 is possible for 4-7-11-13 dihedral angle (but this should be excluded as same direction strand)
					//	example: (in 1QAC chain A) 128.5 is possible for 4-7-10-13 dihedral angle (but this should be excluded as same direction strand)
					//	example: (in 1A3R chain L) 130 is possible for 4-7-10-14 dihedral angle (but this should be excluded as same direction strand)
	max_num_sw_per_pdb_ = tag->getOption<Size>("max_num_sw_per_pdb", 100);
					//	definition: maximum number of sandwiches to be extracted per a pdb file
	check_N_to_C_direction_by_ = tag->getOption<string>("check_N_to_C_direction_by", "PE");
					//	definition: check N->C going direction by option
	check_canonicalness_cutoff_ = tag->getOption<Real>("check_canonicalness_cutoff", 80.0);
					//	definition:	cutoff to determine canonicalness of L/R, P/A and directionality


	///////// An option that takes the longest time ///////
	count_AA_with_direction_ = tag->getOption<bool>("count_AA_with_direction", false);
					//	definition:	if true, count AA considering direction too!
					//	<note> if true, it takes more time, but ~50 sandwiches this can be used within ~ minutes.



	///////// strictness options ///////
	exclude_desinated_pdbs_ = tag->getOption<bool>("exclude_desinated_pdbs", false);
					//	definition: if true, exclude certain designated pdbs
	exclude_sandwich_that_has_near_backbone_atoms_between_sheets_ = tag->getOption<bool>("exclude_sandwich_that_has_near_backbone_atoms_between_sheets", true);
					//	definition: if true, exclude sandwich_that_has_near_backbone_atoms_between_sheets

	max_starting_loop_size_ = tag->getOption<Size>("max_starting_loop_size", 6);
					//	definition: maximum starting loop size to extract
	max_ending_loop_size_ = tag->getOption<Size>("max_ending_loop_size", 6);
					//	definition: maximum ending loop size to extract
	no_helix_in_pdb_ = tag->getOption<bool>("no_helix_in_pdb", false);
					// if true, ignore any pdb that has helix
	max_E_in_extracted_sw_loop_ = tag->getOption<Size>("max_E_in_extracted_sw_loop", 10);
					//	definition: maximum allowable number of E residues in extracted sandwich loop
					//	usefulness: If used, it is useful to exclude [1LOQ] which is a beta-propeller

	max_H_in_extracted_sw_loop_ = tag->getOption<Size>("max_H_in_extracted_sw_loop", 10);
					//	definition: maximum allowable number of helix residues in extracted sandwich loop
					//	example: 0 would be ideal, but then only ~10% of sandwiches will be extracted among CATH classified sandwiches instead even when same_direction_strand linking sw is allowed!

	exclude_sandwich_that_is_linked_w_same_direction_strand_ = tag->getOption<bool>("exclude_sandwich_that_is_linked_w_same_direction_strand", true);
					//	definition: if true, exclude a sandwich that is linked with same_direction_strand
					//	Rationale of default=true (1)
						//	If true, it is useful to exclude [1QAC]_chain_A, [2v33]_chain_A which is a canonical sandwich but linked by same direction strands between sheets

	inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_ = tag->getOption<Real>("inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets", 13.0);
					//	definition: within this distance, sheets are considered to be too near each other
					//	example: (in 1LOQ) inter-sheet distances are 11.5~14.1
					//	Rationale of default value=13 Angstron
						//	it is useful to exclude [1LOQ] which is beta-propeller and [3BVT] which is a stacked sandwich
						//	but it also excludes [2V33] which has two canonical sandwiches near each other and [1W8O] which is a canonical sandwich near a single beta-sheet

	allowed_deviation_for_turn_type_id_ = tag->getOption<Real>("allowed_deviation_for_turn_type_id", 40.0);

	primary_seq_distance_cutoff_for_beta_sheet_capping_ = tag->getOption<int>("primary_seq_distance_cutoff_for_beta_sheet_capping", 2);


	///	electrostatic_interactions related
	distance_cutoff_for_electrostatic_interactions_ = tag->getOption<Real>("distance_cutoff_for_electrostatic_interactions", 7.0);
		//"N-O bridges follow only one, that is, they have at least a pair of side-chain functional-group nitrogen and oxygen atoms within 4A distance, but the side-chain functional-group centroids are > 4A apart."
		// source: 2002_Close-Range Electrostatic Interactions in Proteins

	CB_b_facor_cutoff_for_electrostatic_interactions_ = tag->getOption<Real>("CB_b_facor_cutoff_for_electrostatic_interactions", 60);
		//"Values of 60 or greater may imply disorder (for example, free movement of a side chain or alternative side-chain conformations). Values of 20 and 5 correspond to uncertainties of 0.5 and 0.25 angstroms, respectively."
		// source: http://spdbv.vital-it.ch/TheMolecularLevel/SPVTut/text/STut09aTN.html

	primary_seq_distance_cutoff_for_electrostatic_interactions_ = tag->getOption<Size>("primary_seq_distance_cutoff_for_electrostatic_interactions", 4);
		// rationale for default value: I hypothesize that electrostatic interaction between 38E and 41K of Tencon do not stabilize that much


	///////// development options ///////
	do_not_connect_sheets_by_loops_ = tag->getOption<bool>("do_not_connect_sheets_by_loops", false);
					//	definition: if true, don't connect sheets by loops
	extract_sandwich_ = tag->getOption<bool>("extract_sandwich", true);



	///////// writing options ///////
	write_all_info_files_ = tag->getOption<bool>("write_all_info_files", false);
					//	definition: if true, write all below

	write_AA_kind_files_ = tag->getOption<bool>("write_AA_kind_files", false);
					//	definition: if true, write files that have amino acid kinds
	write_AA_distribution_files_w_direction_ = tag->getOption<bool>("write_AA_distribution_files_w_direction", false);
					//	definition: if true, write files that have amino acid distributions with directions
	write_AA_distribution_files_wo_direction_ = tag->getOption<bool>("write_AA_distribution_files_wo_direction", false);
					//	definition: if true, write files that have amino acid distributions without directions
	write_chain_B_resnum_ = tag->getOption<bool>("write_chain_B_resnum", false);
			// if true, write chain_B_resnum file for InterfaceAnalyzer
	write_phi_psi_of_all_ = tag->getOption<bool>("write_phi_psi_of_all", false);
					//	definition: if true, write phi_psi_file
	write_phi_psi_of_E_ = tag->getOption<bool>("write_phi_psi_of_E", false);
					//	definition: if true, write phi_psi_file
	write_resfile_ = tag->getOption<bool>("write_resfile", false);
					//	definition: if true, write resfile automatically
	write_p_aa_pp_files_ = tag->getOption<bool>("write_p_aa_pp_files", false);
					//	definition: if true, write p_aa_pp_files
	write_rama_at_AA_to_files_ = tag->getOption<bool>("write_rama_at_AA_to_files", false);
					//	definition: if true, write write_rama_at_AA_to_files
	write_heading_directions_of_all_AA_in_a_strand_ = tag->getOption<bool>("write_heading_directions_of_all_AA_in_a_strand", false);
	write_electrostatic_interactions_of_surface_residues_in_a_strand_ = tag->getOption<bool>("write_electrostatic_interactions_of_surface_residues_in_a_strand", false);
	write_electrostatic_interactions_of_all_residues_in_a_strand_ = tag->getOption<bool>("write_electrostatic_interactions_of_all_residues_in_a_strand", false);
	write_electrostatic_interactions_of_all_residues_ = tag->getOption<bool>("write_electrostatic_interactions_of_all_residues", false);
	write_beta_sheet_capping_info_ = tag->getOption<bool>("write_beta_sheet_capping_info", false);
		// reference: 2008_beta-Sheet capping- Signals that initiate and terminate beta-sheet formation, Journal of Structural Biology FarzadFard et al.,
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////				SandwichFeatures			///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///@brief collect all the feature data for the pose
core::Size
SandwichFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
		TR.Info << "======================= <begin> report_features =========================" << endl;

	string tag = get_tag(struct_id, db_session);
	bool canonical_sw_extracted_from_this_pdb_file = false;

	// <begin> exclude_desinated_pdbs
	if (exclude_desinated_pdbs_)
	{
		bool this_pdb_should_be_excluded	=	check_whether_this_pdb_should_be_excluded(tag);
		if (this_pdb_should_be_excluded)
		{
				TR.Info << "Exit early since this pdb should be excluded " << endl;
			return 0;
		}
	}
	// <end> exclude_desinated_pdbs

//	utility::vector1< pose::PoseOP > singlechain_poses;
//	singlechain_poses = pose.split_by_chain();
//		TR << "singlechain_poses.size(): " << singlechain_poses.size() << endl;
//
//	for(Size pose_i=1; pose_i<=singlechain_poses.size(); pose_i++)
//	{
//		pose::Pose pose ( singlechain_poses[ pose_i ] );


	pose::Pose dssp_pose ( pose ); //copy of pose, since the original pose is called as 'const'
	core::scoring::dssp::Dssp dssp( dssp_pose );
	dssp.insert_ss_into_pose( dssp_pose );



	if (no_helix_in_pdb_){
		bool helix_existence = check_helix_existence(dssp_pose);
		if (helix_existence){
				TR.Info << "Exit early since this pdb has helix " << endl;
			return 0;
		}
	}

	if (write_all_info_files_)
	{
		write_AA_kind_files_ = true;
		write_AA_distribution_files_w_direction_ = true;
		write_AA_distribution_files_wo_direction_ = true;
		write_chain_B_resnum_ = true;
		write_phi_psi_of_all_ = true;
		write_phi_psi_of_E_	= true;
		write_resfile_	= true;
		write_p_aa_pp_files_	= true;
		write_rama_at_AA_to_files_	=	true;
		write_heading_directions_of_all_AA_in_a_strand_	=	true;
		write_electrostatic_interactions_of_surface_residues_in_a_strand_	=	true;
		write_electrostatic_interactions_of_all_residues_in_a_strand_	=	true;
		write_electrostatic_interactions_of_all_residues_	=	true;
		write_beta_sheet_capping_info_	=	true;
	}


	pose::Pose pose_w_center_000 = pose;
	if (count_AA_with_direction_ ||	write_resfile_)
	{
		pose_w_center_000.center();
	}


	Size sheet_PK_id_counter=1; //initial value
	Size sw_can_by_sh_PK_id_counter=1; //initial value
	Size rkde_in_strands_PK_id_counter=1; //initial value
	Size rkde_PK_id_counter=1; //initial value
	Size sw_can_by_sh_id_counter=1; //initial value
	Size sw_by_components_PK_id_counter=1; //initial value
	Size intra_sheet_con_id_counter=1; //initial value
	Size inter_sheet_con_id_counter=1; //initial value

	utility::vector1<SandwichFragment> all_strands = get_full_strands(struct_id, db_session);

	if ((static_cast<Size>(all_strands.size()) < min_num_strands_to_deal_) || (static_cast<Size>(all_strands.size()) > max_num_strands_to_deal_)){
			TR.Info << "Exit early since all_strands.size(): " << all_strands.size() << endl;
		return 0;
	}

	// <begin> assignment of strands into the sheet & define very first sheet ("1")
	bool first_sheet_assigned = false;
	for(Size i=1; i<all_strands.size() && !first_sheet_assigned; ++i)	// I don't need the last strand since this double for loops are exhaustive search for all pairs of strands
	{
		if (all_strands[i].get_size() < min_res_in_strand_)
		{
			continue;
		}
		for(Size j=i+1; j<=all_strands.size() && !first_sheet_assigned; ++j)	// I need the last strand for this second for loop
		{
			if (all_strands[j].get_size() < min_res_in_strand_)
			{
				continue;
			}
			SandwichFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
			SandwichFragment temp_strand_j(all_strands[j].get_start(), all_strands[j].get_end());

			// <begin> anti-parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand
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
			if (return_of_find_sheet_parallel == 999)	// since these two strands are too distant to each other, there is virtually no chance to be sheet!
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
	// <end> assignment of strand into sheet & define very first sheet ("1")


	// <begin> assignment of strand into rest sheets (other than "1")
	for(Size i=1; i<=all_strands.size(); ++i)
	{
		if (all_strands[i].get_size() >= min_res_in_strand_) // the length of this beta strand > min_res_in_strand_
		{
			bool strand_i_is_in_any_sheet = check_whether_strand_i_is_in_sheet(struct_id, db_session, get_segment_id(
				struct_id,
				db_session,
				i));
			if (!strand_i_is_in_any_sheet) //  this strand is not in any sheet, so this strand needs to be assigned into any sheet
			{
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
		} //all_strands[i].get_size() >= min_res_in_strand_

		else // all_strands[i].get_size() < min_res_in_strand_  // "this strand is too small, assign it into '99999' sheet"
		{
			write_to_sheet (struct_id, db_session, sheet_PK_id_counter, 99999, get_segment_id(
																							  struct_id,
																							  db_session,
																							  i));	// struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
			sheet_PK_id_counter++;
		}// all_strands[i].get_size() < min_res_in_strand_
	}
	// <end> assignment of strand into rest sheets (other than "1")


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


	// <begin> see_whether_sheet_is_antiparallel
	utility::vector1<Size> all_distinct_sheet_ids = get_distinct_sheet_id_from_sheet_table(struct_id, db_session);
	for(Size i=1; i<=all_distinct_sheet_ids.size(); i++)
	{
		if (all_distinct_sheet_ids[i] == 99999) //all_strands[i].get_size() < min_res_in_strand_
		{
			continue;
		}
		Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

		if (num_strands_i < min_num_strands_in_sheet_)
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


	if (!extract_sandwich_)
	{
			TR.Info << "Exit early since the user doesn't request to extract beta-sandwich " << endl;
		return 0;
	}
	
/////////////////// <begin> assignment of sheet into sandwich_candidate_by_sheet (sw_can_by_sh)
		TR.Info << "<begin> assignment of sheet into sandwich_candidate_by_sheet (sw_can_by_sh) " << endl;

	// <begin> assign num_of_sheets_that_surround_this_sheet only ONCE!
	for(Size i=1; i <= all_distinct_sheet_ids.size(); i++)
	{
		Size num_of_sheets_that_surround_this_sheet	=	cal_num_of_sheets_that_surround_this_sheet(struct_id,	db_session,	dssp_pose,	all_distinct_sheet_ids,	all_distinct_sheet_ids[i]);
		update_num_of_sheets_that_surround_this_sheet (struct_id,	db_session,	all_distinct_sheet_ids[i],	num_of_sheets_that_surround_this_sheet);
	}
	// <end> assign num_of_sheets_that_surround_this_sheet only ONCE!

	Size sheet_j_that_will_be_used_for_pairing_with_sheet_i = 0; // temp value
	for(Size i=1; i <= all_distinct_sheet_ids.size()-1; i++) // now I check all possible combinations
	{
		if (all_distinct_sheet_ids[i] == sheet_j_that_will_be_used_for_pairing_with_sheet_i) // useful to exclude sheet later
		{
			continue;
		}
		if (all_distinct_sheet_ids[i] == 99999) //	all_strands[i].get_size() < min_res_in_strand_
		{
			continue;
		}
		Size num_strands_in_i_sheet = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

		if (num_strands_in_i_sheet < min_num_strands_in_sheet_)
		{
			continue;
		}
		Size num_of_sheets_that_surround_this_sheet	=	get_num_of_sheets_that_surround_this_sheet(struct_id,	db_session,	all_distinct_sheet_ids[i]);

		if (num_of_sheets_that_surround_this_sheet > 1)
		{
			continue; // i
		}
		Real lowest_avg_dis_between_sheets = 9999; //temp value
		Real avg_dis_between_sheets;

		// <begin> identify sheet_j that_will_be_used_for_pairing_with_sheet_i to be a sandwich
		// I need to iterate all 'j' for loop to find the closest sheet from sheet_id (i)
		for(Size j=i+1; j<=all_distinct_sheet_ids.size(); j++)
		{
				//TR << "all_distinct_sheet_ids[j]: " << all_distinct_sheet_ids[j] << endl;
			if (all_distinct_sheet_ids[j] == sheet_j_that_will_be_used_for_pairing_with_sheet_i)	// useful to exclude sheet later
			{
				continue;
			}
			if (all_distinct_sheet_ids[j] == 99999)//	all_strands[i].get_size() < min_res_in_strand_
			{
				continue;
			}
			Size num_strands_in_j_sheet = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]); // struct_id, db_session, sheet_id
			
			if (num_strands_in_j_sheet < min_num_strands_in_sheet_ )
			{
				continue;
			}

			Size num_of_sheets_that_surround_this_sheet	=	get_num_of_sheets_that_surround_this_sheet(struct_id,	db_session,	all_distinct_sheet_ids[j]);

			if (num_of_sheets_that_surround_this_sheet > 1)
			{
				continue;  // j
			}


			/////////////////// DO NOT ERASE THESE TRACERS ///////////
			//				TR << "Now a preliminary candidate of sheet pair to be a sw is identified after checking that these sheets are not surrounded by more than 1 other sheet" << endl;
			//				TR.Info << "sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
			//				TR.Info << "sheet_id (all_distinct_sheet_ids[j]): " << all_distinct_sheet_ids[j] << endl;
			/////////////////// DO NOT ERASE THESE TRACERS ///////////


			utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
			utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[j]);

			

			// <begin> check whether strands are too distant to each other
			bool these_2_sheets_are_too_distant = false; // temporary 'false' designation
			avg_dis_between_sheets = 0;
			for(Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_distant; ++ii) // this '&& !these_2_sheets_are_too_distant' is needed for better performance
			{
				vector<Real> array_avg_dis_sheet_i_j;
				for(Size jj=1; jj<=strands_from_sheet_j.size(); ++jj)
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
				
				for(Size kk=1; kk<=strands_from_sheet_j.size(); ++kk)
				{
					if (shortest_avg_dis_inter_sheet > array_avg_dis_sheet_i_j[kk-1])
					{
						shortest_avg_dis_inter_sheet = array_avg_dis_sheet_i_j[kk-1];
					}
				}
					//TR.Info << "shortest_avg_dis_inter_sheet: " << shortest_avg_dis_inter_sheet << endl;
				if (shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_)
				{
					//	TR.Info << "shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_" << endl;
					these_2_sheets_are_too_distant = true;
				}
			} // for(Size ii=1; ii<=strands_from_sheet_i.size(); ++ii)
			// <end> check whether strands are too distant to each other



				//TR.Info << "avg_dis_between_sheets: " << avg_dis_between_sheets << endl;
			if (these_2_sheets_are_too_distant)
			{
				continue; // continue j sheet loop
			}
			
			if (lowest_avg_dis_between_sheets > avg_dis_between_sheets)
			{
				lowest_avg_dis_between_sheets = avg_dis_between_sheets;

				// <begin> check_strand_too_closeness
					bool these_2_sheets_are_too_close = false; // temporary 'false' designation

					for(Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_close; ++ii) // && !these_2_sheets_are_too_close NEEDED for better performance
					{
						for(Size jj=1; jj<=strands_from_sheet_j.size() && !these_2_sheets_are_too_close; ++jj) // && !these_2_sheets_are_too_close NEEDED for better performance
						{
							SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
							SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());
							
							bool are_strands_too_close = check_strand_too_closeness (pose, temp_strand_i, temp_strand_j);
							if (are_strands_too_close)
							{
							//		TR.Info << "these two sheets are too close when I calculate its distance by their strands" << endl;
								these_2_sheets_are_too_close = true;
							}
						}
					}
					if (these_2_sheets_are_too_close)
					{
						continue; // continue in j sheet loop
					}
				// <end> check_strand_too_closeness

				sheet_j_that_will_be_used_for_pairing_with_sheet_i = all_distinct_sheet_ids[j];	// use all_distinct_sheet_ids[j] to pair with all_distinct_sheet_ids[i]
			}

			
		} // for(Size j=i+1; j<=all_distinct_sheet_ids.size(); ++j)
		
		// <end> identify sheet_j_that_will_be_used_for_pairing_with_sheet_i to be a sandwich


		if (sheet_j_that_will_be_used_for_pairing_with_sheet_i == 0)
		{
			continue; // continue i sheet 'for' loop
		}
		
		/////////////////// DO NOT ERASE THESE TRACERS ///////////
			//			TR.Info << "Now a real pair of candidates between the closest sheets is identified " << endl;
			//			TR.Info << "candidate 1: sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
			//			TR.Info << "candidate 2: sheet_id (sheet_j_that_will_be_used_for_pairing_with_sheet_i): " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << endl;
		/////////////////// DO NOT ERASE THESE TRACERS ///////////

		// <begin> check_sw_by_distance
		bool found_sandwich_w_these_2_sheets = false; // temporary 'false' designation
		bool chance_of_being_sandwich_w_these_2_sheets = true; // temporary 'true' designation
		utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
		utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, sheet_j_that_will_be_used_for_pairing_with_sheet_i);

//			TR.Info << "<begin> check_sw_by_distance" << endl;
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
							//TR.Info << "these sheets will not be sandwich ever because these are too close or distant to each other!" << endl;
						chance_of_being_sandwich_w_these_2_sheets = false;
					}
					
					if ( return_of_check_sw_by_dis_anti != -99 || return_of_check_sw_by_dis_parallel != -99)
					{
						//	TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " are in the ideal distance range" << endl;
						found_sandwich_w_these_2_sheets = true;
						chance_of_being_sandwich_w_these_2_sheets = false; // these are sandwich, but no more sheet search is needed! (this "false" false assignment is needed (confirmed! by experiment))
					}
				} //for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich_w_these_2_sheets; ++jj)
			} //for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich_w_these_2_sheets; ++ii)
			break; // no sandwich here
		} //while (!found_sandwich_w_these_2_sheets && chance_of_being_sandwich_w_these_2_sheets)
//			TR.Info << "<end> check_sw_by_distance" << endl;
		// <end> check_sw_by_distance
		
		if (!found_sandwich_w_these_2_sheets)
		{
			continue;
		}
		
		int facing = judge_facing(struct_id, db_session, pose, all_distinct_sheet_ids[i], sheet_j_that_will_be_used_for_pairing_with_sheet_i);
		// if false, these two strand_pairs are linear to each other or do not face properly to each other
		
		if	(facing == 0)
		{
			//	TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " do not face each other" << endl;
			continue; // skip this sheet
		}
		else if (facing == -99)
		{
				TR.Info << "at least one sheet (either " << all_distinct_sheet_ids[i] << " or " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << ")  may be a beta-barrel like sheet_id = 1 in 1N8O" << endl;
			continue; // skip this sheet
		}
		else

				//	TR.Info << "! writing into 'sandwich candidate by sheet' !" << endl;
		
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
		// It retrieves all beta-strands of sandwich_candidate_by_sheets, it does not make sandwich_by_components

	if (bs_of_sw_can_by_sh.size() == 0)
	{
			TR.Info << "no beta segment in sandwich_by_sheet (maybe these are too distant sheets or a beta barrel or \"non-canonical\" like 1MSP"") " << endl;
			TR.Info << "<Exit-Done> for this pdb including extraction of sandwich" << endl;
		return 0;
	}

	if (write_phi_psi_of_E_)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string phi_psi_file_name = pdb_file_name + "_phi_psi_of_strand_res.txt";
		ofstream phi_psi_file;
		phi_psi_file.open(phi_psi_file_name.c_str());	
		phi_psi_file << "tag	res_num	res_type	res_at_terminal	sheet_is_antiparallel	strand_is_at_edge	phi	psi" << endl;
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii)
		{
			if (bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_)
			{
				break;
			}
			string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
			string strand_is_at_edge = is_this_strand_at_edge	(
										pose,
										struct_id,
										db_session,
										bs_of_sw_can_by_sh[ii].get_sheet_id(),
										bs_of_sw_can_by_sh[ii].get_start(),
										bs_of_sw_can_by_sh[ii].get_end());

			Size component_size = bs_of_sw_can_by_sh[ii].get_size();
			fill_sw_by_components (struct_id, db_session, pose, sw_by_components_PK_id_counter, tag, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), sheet_antiparallel, bs_of_sw_can_by_sh[ii].get_strand_id(), strand_is_at_edge, component_size, bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
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
				Real psi = pose.psi(res_num);
				phi_psi_file << tag << "	" << res_num << "	" << pose.residue_type(res_num).name3() << "	" << res_at_terminal << "	" <<	sheet_antiparallel << "	" << strand_is_at_edge << "	" << phi << "	" << psi << endl;
			}
		}
		phi_psi_file.close();
	} //write_phi_psi_of_E_

	else	// write_phi_psi_of_E_ = false
	{
		if (write_phi_psi_of_all_)
		{
			Size tag_len = tag.length();
			string pdb_file_name = tag.substr(0, tag_len-5);
			string phi_psi_file_name = pdb_file_name + "_phi_psi_of_all_res.txt";
			ofstream phi_psi_file;
			phi_psi_file.open(phi_psi_file_name.c_str());	
			phi_psi_file << "tag	res_num	res_type	dssp	phi	psi" << endl;
			for(core::Size ii=1; ii<=dssp_pose.total_residue(); ii++ )
			{
				char res_ss( dssp_pose.secstruct( ii ) ) ;
				Real phi = pose.phi(ii);
				Real psi = pose.psi(ii);
				phi_psi_file << tag << "	" << ii << "	" << pose.residue_type(ii).name3() << "	" << res_ss	<< "	" << phi << "	" << psi << endl;
			}
		}
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii)
		{
			if (bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_)
			{
				break;
			}
			string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());

			string strand_is_at_edge = is_this_strand_at_edge	(
										pose,
										struct_id,
										db_session,
										bs_of_sw_can_by_sh[ii].get_sheet_id(),
										bs_of_sw_can_by_sh[ii].get_start(),
										bs_of_sw_can_by_sh[ii].get_end());

			Size component_size = bs_of_sw_can_by_sh[ii].get_size();
			fill_sw_by_components (struct_id, db_session, pose, sw_by_components_PK_id_counter, tag, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), sheet_antiparallel, bs_of_sw_can_by_sh[ii].get_strand_id(), strand_is_at_edge, component_size,	bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			sw_by_components_PK_id_counter++;			
		}
	}	//!write_phi_psi_of_E_



	if (count_AA_with_direction_)
	{
		//// <begin> count AA with direction
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii) // per each beta-strand
		{
			if (bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_)
			{
				break;
			}
			update_sw_by_components_by_AA_w_direction (struct_id, db_session, pose, pose_w_center_000,	bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),	bs_of_sw_can_by_sh[ii].get_sheet_id(), bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
		}
		//// <end> count AA with direction
	}

/////////////////// <end> fill a table 'sw_by_components' by secondary_structure_segments
		TR.Info << "<end> fill a table 'sw_by_components' by secondary_structure_segments" << endl;


	if (do_not_connect_sheets_by_loops_)
	{
			TR.Info << "<Exit> do_not_connect_sheets_by_loops_:" << do_not_connect_sheets_by_loops_ << endl;
		return 0;
	}
	
	/////////////////// <begin> update beta-hairpin or inter_sheet_connecting_loops (2nd judgement whether each sandwich_by_sheet_id can become a sandwich_by_components)

	// get_distinct(sw_can_by_sh_id)
	utility::vector1<Size> vec_sw_can_by_sh_id =  get_vec_sw_can_by_sh_id(struct_id, db_session);
		
	for(Size ii=1; ii<=vec_sw_can_by_sh_id.size(); ii++)
	{ // I think that mostly vec_sw_can_by_sh_id.size() = just 1
			TR << "Can sw_candidate_by_sheets_id " << vec_sw_can_by_sh_id[ii] << " be a canonical sw?" << endl;
		bool chance_of_being_canonical_sw	=	true; // not yet decided fate whether this could be canonical sandwich or not, but assumed to be true for now
		Size size_sw_by_components_PK_id =
		get_size_sw_by_components_PK_id(
			struct_id,
			db_session,
			vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
			);

		bool bool_proper_num_helix_in_loop = true;
		bool bool_proper_num_E_in_loop = true;
		Size former_start_res = 0; //temporary

			// this 'jj' is used for 'for' iteration purpose only and this 'for' loop iterates only for connecting sheets/strands
		for(Size jj=1; jj<=size_sw_by_components_PK_id-1; ++jj) {
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

			if (start_res == 0)	// no proper retrieval of start_res
			{
				chance_of_being_canonical_sw = false;
				break; // break jj 'for' loop
			}
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

//////////// <begin> check numbers of helix and strand residues in loops
			//////////// <begin> check whether there is a helix as a loop in this extracted sandwich candidate
			Size helix_num = 0;
			for(Size kk=start_res+1; kk<=next_start_res-1; ++kk){
				char res_ss( dssp_pose.secstruct( kk ) ) ;
				if (res_ss == 'H'){
					helix_num += 1;
				}
			}
			if (helix_num > max_H_in_extracted_sw_loop_){
					//TR << "helix_num > max_H_in_extracted_sw_loop_ " << endl;
				bool_proper_num_helix_in_loop = false;
			}
			//////////// <end> check whether there is a helix as a loop in this extracted sandwich candidate


			//////////// <begin> check whether there is a strand as a loop in this extracted sandwich candidate
			Size E_num = 0;
			for(Size kk=start_res+1; kk<=next_start_res-1; ++kk){
				char res_ss( dssp_pose.secstruct( kk ) ) ;
				if (res_ss == 'E'){
					E_num += 1;
				}
			}
			if (E_num > max_E_in_extracted_sw_loop_){
					TR << "E_num > max_E_in_extracted_sw_loop_ " << endl;
				bool_proper_num_E_in_loop = false;
			}
			//////////// <end> check whether there is a strand as a loop in this extracted sandwich candidate


			if (!bool_proper_num_helix_in_loop || !bool_proper_num_E_in_loop){
				delete_this_sw_can_by_sh_id_from_sw_by_comp(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
				chance_of_being_canonical_sw = false;
				break; // break jj 'for' loop
			}
//////////// <begin> check numbers of helix and strand residues in loops


			string LR = check_LR(dssp_pose, start_res+1, next_start_res-1);

			std::pair<string, string> PA = check_PA(dssp_pose, start_res+1, next_start_res-1);
			string PA_by_preceding_E = PA.first;
			string PA_by_following_E = PA.second;

			// check whethere it is positive, negative, away or meet
			string heading_direction = check_heading_direction(dssp_pose, start_res+1, next_start_res-1, check_N_to_C_direction_by_);

			if (heading_direction == "except")
			{
					TR.Info << "Exit-Exception:: check_N_to_C_direction_by should be either PF or FE !" << endl;
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
				// this loop is a beta-hairpin loop (that connects sheets as intra-sheet way)
			{
//					TR << "start_res: " << start_res << endl;
//					TR << "next_start_res: " << next_start_res << endl;	
				bool loop_is_surrounded_by_same_direction_strands = check_whether_same_direction_strands_connect_two_sheets_or_a_loop(struct_id,	db_session,	pose,	start_res,	next_start_res);
					//TR.Info << "loop_is_surrounded_by_same_direction_strands: " << loop_is_surrounded_by_same_direction_strands << endl;

				string canonical_LR = "-";
				string cano_PA =  "-";
				string cano_parallel_EE =  "-";
				string loop_kind =  "loop_connecting_same_direction_strands_within_same_sheet";

				bool	hairpin_connects_short_strand = check_whether_hairpin_connects_short_strand(struct_id,	db_session,	start_res,	next_start_res);

				if (!loop_is_surrounded_by_same_direction_strands)
				{
					canonical_LR = check_canonicalness_of_LR(loop_size, true, LR); // loop_size, intra_sheet bool, LR
					cano_PA = check_canonicalness_of_PA(loop_size, true, PA_by_preceding_E, PA_by_following_E, check_canonicalness_cutoff_); // loop_size, intra_sheet bool, 2 PAs
					cano_parallel_EE = check_canonicalness_of_parallel_EE(loop_size, true, parallel_EE); // loop_size, intra_sheet bool, parallel_EE

					if (hairpin_connects_short_strand) // so it is not recommended for LR/PA analysis
					{
						loop_kind =  "hairpin_connecting_short_strand";
					}
					else // this hairpin-loop is ideal for LR/PA
					{
						loop_kind =  "hairpin";
					}

				}

				update_sheet_connectivity(
					struct_id,
					db_session,
					pose,
					sw_by_components_PK_id_counter,
					tag,
					vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
					loop_kind, // "hairpin" or "a_loop_that_connects_same_direction_strands_within_same_sheet"
					intra_sheet_con_id_counter,
					inter_sheet_con_id_counter,
					LR,
					canonical_LR,
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

				if (!loop_is_surrounded_by_same_direction_strands	&&	!hairpin_connects_short_strand) // so it is recommended for LR/PA analysis
				{
					if (((next_start_res-1) - (start_res+1)) == 1) // this loop is a 2-residues long hairpin
					{
						string turn_type	=	report_turn_type
												(pose,
												vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
												start_res+1,	// new start_res for intra_sheet_con
												next_start_res-1, // new end_res for intra_sheet_con
												struct_id,
												db_session
												);

						report_turn_AA
							(pose,
							vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
							start_res,	// i
							struct_id,
							db_session,
							turn_type
							);
					}
				}

				sw_by_components_PK_id_counter++;
				intra_sheet_con_id_counter++;
			}
			else // this loop connects sheets as inter-sheet way
			{
				if (exclude_sandwich_that_is_linked_w_same_direction_strand_)
				{
					bool sheets_are_connected_by_same_direction_strand = check_whether_same_direction_strands_connect_two_sheets_or_a_loop(struct_id,	db_session,	pose,	start_res,	next_start_res);
					if (sheets_are_connected_by_same_direction_strand)
					{
							TR.Info << "sheets_are_connected_by_same_direction_strand, so delete " << vec_sw_can_by_sh_id[ii] << "sw_can_by_sh_id" << endl;
						delete_this_sw_can_by_sh_id_from_sw_by_comp(
							struct_id,
							db_session,
							vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
							);
						chance_of_being_canonical_sw = false;
						break; // break jj 'for' loop
					}
				}


				string canonical_LR = check_canonicalness_of_LR(loop_size, false, LR);	// loop_size, intra_sheet bool, LR
				string cano_PA = check_canonicalness_of_PA(loop_size, false, PA_by_preceding_E, PA_by_following_E, check_canonicalness_cutoff_);
														  // loop_size,	intra_sheet bool, PA_ref_1, PA_ref_2, cutoff
				string cano_parallel_EE = check_canonicalness_of_parallel_EE(loop_size, false, parallel_EE);
														  // loop_size, intra_sheet bool, parallel_EE
				update_sheet_connectivity(
					struct_id,
					db_session,
					pose,
					sw_by_components_PK_id_counter,
					tag,
					vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
					"loop_connecting_two_sheets",
					intra_sheet_con_id_counter,
					inter_sheet_con_id_counter,
					LR,
					canonical_LR,
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

				if (((next_start_res-1) - (start_res+1)) == 1) // this loop is a 2-residues long inter-sheet_loop
				{
					string	turn_type	=	report_turn_type
												(pose,
												vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
												start_res+1,	// new start_res for intra_sheet_con
												next_start_res-1, // new end_res for intra_sheet_con
												struct_id,
												db_session
												);

					report_turn_AA
						(pose,
						vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
						start_res,	// i
						struct_id,
						db_session,
						turn_type
						);
				}

				sw_by_components_PK_id_counter++;
				inter_sheet_con_id_counter++;
			}	// this loop connects sheets as inter-sheet way

		} // for(Size jj=1; (jj<=size_sw_by_components_PK_id-1) && (bool_no_helix_in_loop) && (bool_no_more_strand_in_loop); ++jj)
	/////////////////// <end> update beta-hairpin or inter_sheet_connecting_loops (2nd judgement whether each sandwich_by_sheet_id becomes sandwich_by_components)
	

		if (exclude_sandwich_that_has_near_backbone_atoms_between_sheets_)
		{
			bool sheets_are_connected_with_near_bb_atoms = check_whether_sheets_are_connected_with_near_bb_atoms(struct_id,	db_session,	dssp_pose,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
			if (sheets_are_connected_with_near_bb_atoms)
			{
					TR.Info << "sheets are connected with near bb atoms, so delete " << vec_sw_can_by_sh_id[ii] << " sw_can_by_sh_id" << endl;
				delete_this_sw_can_by_sh_id_from_sw_by_comp(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
				chance_of_being_canonical_sw = false;
			}
		}

			TR.Info << "chance_of_being_canonical_sw: " << chance_of_being_canonical_sw << endl;

		if (chance_of_being_canonical_sw)
		{
			canonical_sw_extracted_from_this_pdb_file = true;

			add_starting_loop(struct_id,	db_session,	dssp_pose,	sw_by_components_PK_id_counter,	vec_sw_can_by_sh_id[ii],	tag);
				//	struct_id,	db_session,	dssp_pose,	sw_by_components_PK_id, sw_can_by_sh_id, tag				
			sw_by_components_PK_id_counter++;

			add_ending_loop(struct_id,	db_session,	dssp_pose,	sw_by_components_PK_id_counter,	vec_sw_can_by_sh_id[ii],	tag);
				//	struct_id,	db_session,	dssp_pose,	sw_by_components_PK_id, sw_can_by_sh_id, tag
			sw_by_components_PK_id_counter++;


			// <begin> mark beta-sandwiches that is not connected by continuous atoms like 1A78
				string sw_is_not_connected_with_continuous_atoms = check_whether_sw_is_not_connected_with_continuous_atoms(
					struct_id,
					db_session,
					dssp_pose,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
					//TR.Info << "sw_is_not_connected_with_continuous_atoms: " << sw_is_not_connected_with_continuous_atoms << endl;
				mark_sw_which_is_not_connected_with_continuous_atoms(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
					sw_is_not_connected_with_continuous_atoms);
			// <end> mark beta-sandwiches that is not connected by continuous atoms like 1A78



			add_dssp_ratio_in_sw(struct_id,	db_session,	dssp_pose,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

			if (count_AA_with_direction_)
			{
				add_number_of_inward_pointing_W_in_sw(struct_id,	db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);

				add_number_of_inward_pointing_LWY_in_core_strands_in_sw(struct_id,	db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
			}

			Size	sw_res_size	=	add_sw_res_size(struct_id,	db_session,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

			add_num_strands_in_each_sw(struct_id,	db_session,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

			add_num_edge_strands_in_each_sw(struct_id,	db_session,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);
				
			report_hydrophobic_ratio_net_charge(struct_id,	db_session,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

			report_dihedral_angle_between_core_strands_across_facing_sheets(struct_id,	db_session, pose,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);


			report_avg_b_factor_CB_at_each_component(struct_id,	db_session, pose,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

			if (write_AA_kind_files_)
			{
				// <begin> write AA_kind to a file
				Size tag_len = tag.length();
				string pdb_file_name = tag.substr(0, tag_len-5);

				string AA_kind_file_name = pdb_file_name + "_AA_kind.txt";
				ofstream AA_kind_file;
				
				AA_kind_file.open(AA_kind_file_name.c_str());
				
				utility::vector1<Size> vec_AA_kind = get_vec_AA_kind(struct_id,	db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);

					// positive, negative, polar, aromatic, hydrophobic
				AA_kind_file << "pdb_name	Pos_percent	Neg_percent	Pol_percent	Aro_percent	Hydropho_percent	" ; // as Jenny's thesis
				AA_kind_file << "Pos_raw_count	Neg_raw_count	Pol_raw_count	Aro_raw_count	Hydropho_raw_count" << endl; // to calculate net_charge_of_sw

				AA_kind_file << pdb_file_name << "	" ;

				for (Size i =1; i<=(vec_AA_kind.size()); i++)
				{
					Real percent = vec_AA_kind[i]*100/static_cast<Real>(sw_res_size);
					Size rounded_percent = round_to_Size(percent);
					AA_kind_file << rounded_percent << "	" ;
				}

				for (Size i =1; i<=(vec_AA_kind.size()); i++)
				{
					AA_kind_file << vec_AA_kind[i] << "	" ;
				}

				AA_kind_file << endl; // to marshal
				
				AA_kind_file.close();
				// <end> write AA_kind to a file
			}
	
		}
		else // chance_of_being_canonical_sw = false
		{
			delete_this_sw_can_by_sh_id_from_sw_by_comp(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
		}

		// <begin> write chain_B_resNum to a file
		if (write_chain_B_resnum_ && chance_of_being_canonical_sw)
		{
			Size tag_len = tag.length();
			string pdb_file_name = tag.substr(0, tag_len-5);
			string report_file_name = pdb_file_name + "_chain_B_resNum.txt";
			ofstream report_file;
			report_file.open(report_file_name.c_str());
			utility::vector1<SandwichFragment> chain_B_resNum	= get_chain_B_resNum
																	(
																		struct_id,
																		db_session,
																		vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
																	);
			for(Size i=1; i<=chain_B_resNum.size(); ++i)
			{
				report_file << chain_B_resNum[i].get_resNum() << endl;
			}
			report_file.close();
		}
		// <end> write chain_B_resNum to a file

	}	// per each sandwich_candidate_by_sheet_id


	/// <begin> write_heading_direction_of_all_AA_in_a_strand
	if (write_heading_directions_of_all_AA_in_a_strand_ && canonical_sw_extracted_from_this_pdb_file)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string heading_file_name = pdb_file_name + "_heading_direction_of_all_AA_in_a_strand.txt";
		ofstream heading_file;
			
		heading_file.open(heading_file_name.c_str());
		heading_file << "residue_begin	residue_end	heading_directions" << endl;

		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++) // per each beta-strand
		{
			string	heading_directions	=	report_heading_directions_of_all_AA_in_a_strand (struct_id, db_session, pose, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),	bs_of_sw_can_by_sh[ii].get_sheet_id(), bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			heading_file << bs_of_sw_can_by_sh[ii].get_start() << "	" << bs_of_sw_can_by_sh[ii].get_end() << "	"	<<	heading_directions	<< endl;
		}

		heading_file.close();
	}
	/// <end> write_heading_direction_of_all_AA_in_a_strand


	/// <begin> write_electrostatic_interactions_of_surface_residues_in_a_strand
	if (canonical_sw_extracted_from_this_pdb_file)
	{
		if (write_electrostatic_interactions_of_surface_residues_in_a_strand_	||	write_electrostatic_interactions_of_all_residues_in_a_strand_)
		{
			// <begin> store RKDE in strands to a database table
			for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++) // per each beta-strand
			{
				Size residue_begin	=	bs_of_sw_can_by_sh[ii].get_start();
				Size residue_end	=	bs_of_sw_can_by_sh[ii].get_end();
				for (Size	residue_num	=	residue_begin;	residue_num	<=	residue_end; residue_num++)
				{
					if (
						(pose.residue_type(residue_num).name3() == "ARG")
						||	(pose.residue_type(residue_num).name3() == "LYS")
						||	(pose.residue_type(residue_num).name3() == "ASP")
						||	(pose.residue_type(residue_num).name3() == "GLU")
						)
					{
						string	heading	=	determine_heading_direction_by_vector	(struct_id,	db_session,	pose,	bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),	bs_of_sw_can_by_sh[ii].get_sheet_id(),	residue_begin,	residue_end,	residue_num);
						update_rkde_in_strands(
								struct_id,
								db_session,
								rkde_in_strands_PK_id_counter,
								tag,
								bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),	//sw_can_by_sh_id
								residue_num,	//residue_number,
								pose.residue_type(residue_num).name3(),	//residue_type,
								heading // "surface" or "core"
								);
						rkde_in_strands_PK_id_counter++;
					}
				}
			}
			// <end> store RKDE in strands to a database table
		}

		if (write_electrostatic_interactions_of_surface_residues_in_a_strand_)
		{
			// <begin> report number of electrostatic_interactions_of_surface_residues_in_a_strand
			report_number_of_electrostatic_interactions_of_residues(tag,	struct_id,	db_session,	pose, "E",	"surface");
			// <end> report number of electrostatic_interactions_of_surface_residues_in_a_strand
		}

		if (write_electrostatic_interactions_of_all_residues_in_a_strand_)
		{
			// <begin> report number of electrostatic_interactions_of_all_residues_in_a_strand
			report_number_of_electrostatic_interactions_of_residues(tag,	struct_id,	db_session,	pose, "E",	"all_direction");
			// <end> report number of electrostatic_interactions_of_all_residues_in_a_strand
		}

		if (write_electrostatic_interactions_of_all_residues_)
		{
			// <begin> store RKDE to a database table
			for(Size ii=1; ii<=pose.total_residue(); ii++ )
			{
				if (
					(pose.residue_type(ii).name3() == "ARG")
					||	(pose.residue_type(ii).name3() == "LYS")
					||	(pose.residue_type(ii).name3() == "ASP")
					||	(pose.residue_type(ii).name3() == "GLU")
					)
				{
					update_rkde(
							struct_id,
							db_session,
							rkde_PK_id_counter,
							tag,
							ii,	//residue_number,
							pose.residue_type(ii).name3()	//residue_type,
							);
					rkde_PK_id_counter++;
				}
			}
			// <end> store RKDE to a database table

			// <begin> report number of electrostatic_interactions_of_all_residues
			report_number_of_electrostatic_interactions_of_residues(tag,	struct_id,	db_session,	pose, "all_dssp", "all_direction");
			// <end> report number of electrostatic_interactions_of_all_residues

		}

	}
	/// <end> write_electrostatic_interactions_of_surface_residues_in_a_strand


	// <begin> write_beta_sheet_capping_info_
	if (write_beta_sheet_capping_info_ && canonical_sw_extracted_from_this_pdb_file)
	{
		int	number_of_strands_where_capping_is_checked	=	0;
		int	number_of_strands_with_2_cappings	=	0;
		int	number_of_strands_with_1_capping	=	0;

		utility::vector1<Size> residue_begin_of_strands_with_1_capping;
		utility::vector1<Size> residue_begin_of_strands_with_0_capping;

		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ii++) // per each beta-strand
		{
			int residue_begin	=	bs_of_sw_can_by_sh[ii].get_start();
			int residue_end	=	bs_of_sw_can_by_sh[ii].get_end();
			if (residue_end	-	residue_begin	<	2)
			{
				continue;	// this strand should be too short like "short_edge"
			}
			number_of_strands_where_capping_is_checked++;

			bool	capping_near_begin	=	false;
			bool	capping_near_end	=	false;
			for(int jj	=	residue_begin	-	primary_seq_distance_cutoff_for_beta_sheet_capping_;
				jj	<=	residue_begin	+	primary_seq_distance_cutoff_for_beta_sheet_capping_;
				jj++) // per each beta-strand
			{
				if	(jj	<=	0)
				{
					capping_near_begin	=	true;	// I assume that capping_near_begin is not necessary for the first strand
					break;
					//continue;
				}
				if (
				(pose.residue_type(jj).name3() == "GLY")
				||	(pose.residue_type(jj).name3() == "ASP")
				||	(pose.residue_type(jj).name3() == "ASN")
				||	(pose.residue_type(jj).name3() == "PRO")
				)
				{
					capping_near_begin	=	true;
					break;
				}
			}
			for(int jj	=	residue_end	-	primary_seq_distance_cutoff_for_beta_sheet_capping_;
				jj	<=	residue_end	+	primary_seq_distance_cutoff_for_beta_sheet_capping_;
				jj++) // per each beta-strand
			{
				if	(jj	>	static_cast<int>(pose.total_residue()))
				{
					capping_near_end	=	true;	// I assume that capping_near_end is not necessary for the last strand
					break;
					//continue;
				}
				if (
				(pose.residue_type(jj).name3() == "GLY")
				||	(pose.residue_type(jj).name3() == "ASP")
				||	(pose.residue_type(jj).name3() == "ASN")
				||	(pose.residue_type(jj).name3() == "PRO")
				)
				{
					capping_near_end	=	true;
					break;
				}
			}
			if (capping_near_begin	&&	capping_near_end)
			{
				number_of_strands_with_2_cappings++;
			}
			else if	(capping_near_begin	||	capping_near_end)
			{
				number_of_strands_with_1_capping++;
				residue_begin_of_strands_with_1_capping.push_back(residue_begin);
			}
			else
			{
				residue_begin_of_strands_with_0_capping.push_back(residue_begin);
			}
		}
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string info_file_name = pdb_file_name + "_beta_sheet_capping_info.txt";
		ofstream info_file;

		info_file.open(info_file_name.c_str());
		info_file << "number_of_strands_where_capping_is_checked	number_of_strands_with_2_cappings	number_of_strands_with_1_capping		number_of_strands_with_0_capping	residue_begin_of_strands_with_1_capping	residue_begin_of_strands_with_0_capping" << endl;

		info_file
			<< number_of_strands_where_capping_is_checked << "	"
			<< number_of_strands_with_2_cappings	<<	"	"
			<< number_of_strands_with_1_capping	<<	"	"
			<<	number_of_strands_where_capping_is_checked	-	number_of_strands_with_2_cappings	-	number_of_strands_with_1_capping	<<	"	";

		for	(Size	kk	=	1;		kk	<=	residue_begin_of_strands_with_1_capping.size();	kk++)
		{
			if (kk	==	residue_begin_of_strands_with_1_capping.size())
			{
				info_file <<	residue_begin_of_strands_with_1_capping[kk];
			}
			else
			{
				info_file <<	residue_begin_of_strands_with_1_capping[kk]	<<	",";
			}
		}

		info_file << "	";

		for	(Size	kk	=	1;		kk	<=	residue_begin_of_strands_with_0_capping.size();	kk++)
		{
			if (kk	==	residue_begin_of_strands_with_0_capping.size())
			{
				info_file <<	residue_begin_of_strands_with_0_capping[kk];
			}
			else
			{
				info_file <<	residue_begin_of_strands_with_0_capping[kk]	<<	",";
			}
		}

		info_file << endl;
		info_file.close();
	}
	// <end> write_beta_sheet_capping_info_


	// <begin> write AA_dis to a file
	if (write_AA_distribution_files_w_direction_ && canonical_sw_extracted_from_this_pdb_file)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string AA_dis_file_name = pdb_file_name + "_AA_distribution_w_direction_sorted_alphabetically.txt";
		ofstream AA_dis_file;
		
		AA_dis_file.open(AA_dis_file_name.c_str());	
		utility::vector1<Size> vec_core_heading_at_core_strand = get_vector_AA_distribution_w_direction (struct_id,	db_session, "core_heading", "core");
																							// struct_id,	db_session, heading_direction, strand_location
		utility::vector1<Size> vec_surface_heading_at_core_strand = get_vector_AA_distribution_w_direction (struct_id,	db_session, "surface_heading", "core");
																							// struct_id,	db_session, heading_direction, strand_location
		utility::vector1<Size> vec_core_heading_at_edge_strand = get_vector_AA_distribution_w_direction (struct_id,	db_session, "core_heading", "edge");
																							// struct_id,	db_session, heading_direction, strand_location
		utility::vector1<Size> vec_surface_heading_at_edge_strand = get_vector_AA_distribution_w_direction (struct_id,	db_session, "surface_heading", "edge");
																							// struct_id,	db_session, heading_direction, strand_location
		AA_dis_file << "core_heading_at_core_strand	surface_heading_at_core_strand	core_heading_at_edge_strand	surface_heading_at_edge_strand" << endl;
		for (Size i =1; i<=(vec_core_heading_at_core_strand.size()); i++)
		{
			AA_dis_file << vec_core_heading_at_core_strand[i] << "	" << vec_surface_heading_at_core_strand[i] << "	" << vec_core_heading_at_edge_strand[i] << "	" << vec_surface_heading_at_edge_strand[i] << endl;
		}
		AA_dis_file.close();
	} 
	// <end> write AA_dis to a file



	core::scoring::ScoreFunctionOP centroid_scorefxn( generate_scorefxn( false /*fullatom*/ ) );
	core::scoring::ScoreFunctionOP fullatom_scorefxn( generate_scorefxn( true /*fullatom*/ ) );

	process_decoy( dssp_pose, pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn);

	TR.Info << "Current energy terms that we deal with:"	<< endl;
	dssp_pose.energies().show_total_headers( std::cout );
	std::cout << std::endl;

	// <begin> write p_aa_pp (Probability of amino acid at phipsi) to a file (ref. https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/score_types.html )
	if (write_p_aa_pp_files_ && canonical_sw_extracted_from_this_pdb_file)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string p_aa_pp_file_name = pdb_file_name + "_p_aa_pp_at_each_AA.txt";
		ofstream p_aa_pp_file;
		
		p_aa_pp_file.open(p_aa_pp_file_name.c_str());
		p_aa_pp_file << "residue_number	res_type	p_aa_pp" << endl;

		for(Size ii=1; ii<=dssp_pose.total_residue(); ii++ )
		{
			core::scoring::EnergyMap em1 = dssp_pose.energies().residue_total_energies(ii);
			Real resi_p_aa_pp = em1[core::scoring::p_aa_pp];

			p_aa_pp_file << ii << "	" << dssp_pose.residue_type(ii).name3()	<<	"	"	<<	resi_p_aa_pp << endl;
		}
		p_aa_pp_file.close();
	}
	// <end> write p_aa_pp (Probability of amino acid at phipsi) to a file (ref. https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/score_types.html )


	// <begin> write rama (ramachandran preferences) of residue to a file (ref. https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/score_types.html )
	if (write_rama_at_AA_to_files_ && canonical_sw_extracted_from_this_pdb_file)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string rama_file_name = pdb_file_name + "_rama_at_each_AA.txt";
		ofstream rama_file;
		
		rama_file.open(rama_file_name.c_str());
		rama_file << "residue_number		res_type	rama" << endl;

		for(Size ii=1; ii<=dssp_pose.total_residue(); ii++ )
		{
			core::scoring::EnergyMap em1 = dssp_pose.energies().residue_total_energies(ii);
			Real rama_at_this_AA = em1[core::scoring::rama];

			rama_file << ii << "	" << dssp_pose.residue_type(ii).name3()	<<	"	"	<<	rama_at_this_AA << endl;
		}
		rama_file.close();
	}
	// <end> write rama (ramachandran preferences) of residue to a file (ref. https://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/score_types.html )


	// <begin> write AA_distribution_without_direction to a file
	if (write_AA_distribution_files_wo_direction_ && canonical_sw_extracted_from_this_pdb_file)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string AA_dis_file_name = pdb_file_name + "_AA_distribution_wo_direction_sorted_alphabetically.txt";
		ofstream AA_dis_file;
		
		AA_dis_file.open(AA_dis_file_name.c_str());	
		utility::vector1<Size> vector_of_hairpin_AA = get_vector_AA_distribution_wo_direction (struct_id,	db_session, "hairpin");
		utility::vector1<Size> vector_of_inter_sheet_loop_AA = get_vector_AA_distribution_wo_direction (struct_id,	db_session, "loop_connecting_two_sheets");
		
		AA_dis_file << "hairpin_AA	inter_sheet_loop_AA" << endl;
		for (Size i =1; i<=(vector_of_hairpin_AA.size()); i++)
		{
			AA_dis_file << vector_of_hairpin_AA[i] << "	" << vector_of_inter_sheet_loop_AA[i] << endl;
		}
		AA_dis_file.close();
	} 
	// <end> write AA_distribution_without_direction to a file
	

	//// <begin> report number_of_inward_pointing_charged_AAs/aro_AAs_in_a_pair_of_edge_strands
	if (count_AA_with_direction_ && canonical_sw_extracted_from_this_pdb_file)
	{
		for(Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii)
		{
			bool sw_by_sh_id_still_alive =	check_whether_sw_by_sh_id_still_alive (struct_id, db_session,	bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id());
			if (!sw_by_sh_id_still_alive)
			{
				continue;
			}

			if (bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_)
			{
				break;
			}
			string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
			string strand_is_at_edge = is_this_strand_at_edge	(
										pose,
										struct_id,
										db_session,
										bs_of_sw_can_by_sh[ii].get_sheet_id(),
										bs_of_sw_can_by_sh[ii].get_start(),
										bs_of_sw_can_by_sh[ii].get_end());
			if (strand_is_at_edge == "edge")
			{
				std::pair<Size, Size> current_bs_id_and_closest_edge_bs_id_in_different_sheet	=	get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (struct_id, db_session, pose,	bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),	bs_of_sw_can_by_sh[ii].get_sheet_id(), bs_of_sw_can_by_sh[ii].get_start(),	bs_of_sw_can_by_sh[ii].get_end());
				Size current_bs_id = current_bs_id_and_closest_edge_bs_id_in_different_sheet.first;
				Size closest_bs_id = current_bs_id_and_closest_edge_bs_id_in_different_sheet.second;

				report_number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands (struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), current_bs_id,	closest_bs_id);
				report_number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands (struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), current_bs_id,	closest_bs_id);
			}
		}
	}
	//// <end> report number_of_inward_pointing_charged_AAs/aro_AAs_in_a_pair_of_edge_strands


	///////////// development
	// <begin> write resfile automatically
	if (write_resfile_ && canonical_sw_extracted_from_this_pdb_file)
	{
		Size tag_len = tag.length();
		string pdb_file_name = tag.substr(0, tag_len-5);
		string resfile_name = pdb_file_name + "_resfile.txt";
		ofstream resfile_stream;
		
		resfile_stream.open(resfile_name.c_str());
		
		resfile_stream << "EX 1 NOTAA C" << endl;
		resfile_stream << "USE_INPUT_SC" << endl;
		resfile_stream << "start" << endl;

		for (Size i =1; i<=(pose.total_residue()); i++)
		{
			string residue_location = get_residue_location (struct_id,	db_session,	i);
//			resfile_stream << vector_of_hairpin_AA[i] << "	" << vector_of_inter_sheet_loop_AA[i] << endl;
		}
		resfile_stream.close();
	} 
	// <end> write resfile automatically
	///////////// development


		TR.Info << "<Exit-Done> for this pdb including extraction of sandwich" << endl;
//	}
	return 0;
} //SandwichFeatures::report_features


} //namespace strand_assembly
} //namespace features
} //namespace protocols
