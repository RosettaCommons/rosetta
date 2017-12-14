// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/SandwichFeatures.cc
/// @brief Extract and analyze beta-sandwich features
/// @author Doo Nam Kim (doonam.kim@gmail.com)
/// @overview
///  @ task 0: Determine whether we deal with given pdb file
///  @ task 1: Identify all beta-strands
///  @ task 2: Identify all beta-sheets with these strands
///  @ task 3: Identify all beta-sandwiches with these sheets
///   @ task 3-1: Merge sheets to each other if possible
///   @ task 3-2: Make beta-sandwiches with sheets that are ideal only
///    @ task 3-2-1: Exclude if this_sheet_is_surrounded_by_more_than_1_other_sheet
///    @ task 3-2-2: Exclude sheets that are too close to each other
///    @ task 3-2-3: Exclude sheets that are too distant to each other
///    @ task 3-2-4: Exclude sheets that do not face each other
///     @ task 3-2-4-1: Exclude sheets that do not face each other by an angle with two terminal residues and one central residue
///     @ task 3-2-4-2: Exclude sheets that do not face each other by an angle with four terminal residues in two edge strands
///   @ task 3-4: Test canonical sandwich test
///    @ task 3-4-1: Canonical sandwiches need to have low number of helix or strand residues in any loop (beta-hairpin-loop or inter-sheet-loop)
///    @ task 3-4-2: Canonical sandwiches need to not have same-direction-strands as connecting two beta-sheets
///    @ task 3-4-3: Canonical sandwiches should not be beta-barrel obviously
///  @ task 4: Write beta-sandwiches that passed canonical tests into database
///    @ task 4-1: Write hairpin_loop and inter-sheet loop
///    @ task 4-2: Write starting_loop and endng_loop
///    @ task 4-3: Write ratio_hydrophobic_philic/net_charge
///    @ task 4-4: Write total size of sandwich
///  @ task 5: Write resfiles automatically

//Devel
#include <protocols/features/strand_assembly/CheckForSandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.hh>
#include <protocols/features/strand_assembly/WriteToFileFromSandwichFeatures.hh>

#include <basic/database/schema_generator/DbDataType.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/strand_assembly/SandwichFeaturesCreator.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

static basic::Tracer TR( "protocols.features.strand_assembly.SandwichFeatures" );

// for parse_my_tag
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
using utility::tag::TagCOP;

using namespace std;
using namespace core;
using core::pose::Pose;
using cppdb::result;
using cppdb::statement;
using utility::vector1;
using utility::sql_database::sessionOP;

//using core::id::NamedAtomID;
//using numeric::xyzVector;

SandwichFeatures::SandwichFeatures()
{
	//TR << "A SandwichFeatures constructor called" << endl;
}

SandwichFeatures::~SandwichFeatures()
{
	//TR << "A SandwichFeatures destructor called" << endl;
}

utility::vector1<std::string>
SandwichFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;

	// needed for SandwichFeatures!
	dependencies.push_back("ResidueFeatures"); // needed for ResidueSecondaryStructureFeatures
	dependencies.push_back("ResidueSecondaryStructureFeatures"); // needed for SecondaryStructureSegmentFeatures
	dependencies.push_back("SecondaryStructureSegmentFeatures");


	// not needed for SandwichFeatures, but for format_converter
	dependencies.push_back("ResidueFeatures"); // needed for ProteinResidueConformationFeatures & ResidueSecondaryStructureFeatures
	dependencies.push_back("PoseConformationFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");

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
	Column sheet_PK_id ("sheet_PK_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column sheet_id ("sheet_id", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*no autoincrement*/);
	// <history> changed into 'null-possible' because of sandwich

	Column sheet_antiparallel ("sheet_antiparallel", DbDataTypeOP( new DbText() ), true /* could be null at first, eventually it will not be null though*/, false /*no autoincrement*/);
	// A: antiparallel
	// P_or_mix: parallel or mix of antiparallel and parallel
	Column num_of_sheets_that_surround_this_sheet ("num_of_sheets_that_surround_this_sheet", DbDataTypeOP( new DbInteger() ), true /* could be null */, false /*no autoincrement*/);

	// unique key of original PDB file
	Column struct_id             ("struct_id",              DbDataTypeOP( new DbBigInt() ),    false /*not null*/, false /*don't autoincrement*/);

	// ForeignKey
	Column segment_id ("segment_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);

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
	sheet.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols_beta;
	fkey_reference_cols_beta.push_back("struct_id");
	fkey_reference_cols_beta.push_back("segment_id");

	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(segment_id);

	sheet.add_foreign_key(ForeignKey(fkey_cols, "secondary_structure_segments", fkey_reference_cols_beta, true /*defer*/));

	sheet.write(db_session);
	/****** <end> writing sheet ******/


	/****** <begin> writing sw_can_by_sh (sandwich candidate by sheets) ******/
	// Columns
	// id of sw_can_by_sh
	//unique
	Column sw_can_by_sh_PK_id ("sw_can_by_sh_PK_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	//Column tag ("tag", new DbText(), false /*not null*/, false /*no autoincrement*/);
	Column tag ("tag", DbDataTypeOP( new DbText() ), true /*could be null*/, false /*no autoincrement*/);

	// could be redundant
	Column sw_can_by_sh_id ("sw_can_by_sh_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);
	//Column sw_can_by_sh_id ("sw_can_by_sh_id", new DbInteger(), true /*could be null*/, false /*no autoincrement*/);

	// could be redundant
	Column strand_num ("strand_num", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

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
	sw_can_by_sh.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols_sh;
	fkey_reference_cols_sh.push_back("struct_id");
	fkey_reference_cols_sh.push_back("sheet_PK_id");
	//fkey_reference_cols_sh.push_back("sheet_id"); // sqlite3 OK but psql "ERROR: cppdb::posgresql: statement execution failed : ERROR:  there is no unique constraint matching given keys for referenced table "sheet""

	utility::vector1<Column> fkey_cols_sh;
	fkey_cols_sh.push_back(struct_id);
	fkey_cols_sh.push_back(sheet_id);

	sw_can_by_sh.add_foreign_key(ForeignKey(fkey_cols_sh, "sheet", fkey_reference_cols_sh, true /*defer*/));

	// add column which is neither PrimaryKey nor ForeignKey
	sw_can_by_sh.add_column(strand_num);

	sw_can_by_sh.write(db_session);
	/****** <end> writing sw_can_by_sh ******/


	/****** <begin> writing sandwich (sandwich candidate by components such as strands, loops, helices) ******/

	// Columns
	// id of sandwich

	//unique
	Column sandwich_PK_id ("sandwich_PK_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

	// could be redundant
	Column sandwich_bs_id ("sandwich_bs_id", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column long_strand_id ("long_strand_id", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column strand_edge ("strand_edge", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	// edge strand
	// core strand

	Column num_strands_in_each_sw ("num_strands_in_each_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column num_edge_strands_in_each_sw ("num_edge_strands_in_each_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column intra_sheet_con_id ("intra_sheet_con_id", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column inter_sheet_con_id ("inter_sheet_con_id", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column loop_kind ("loop_kind", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	// starting_loop
	// hairpin_loop (intra-sheet loop)
	// inter_sheet_loop
	// ending_loop

	Column LR ("LR", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);

	Column canonical_LR ("canonical_LR", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	// T, -> true, canonical chiral
	// F, -> false, non-canonical chiral
	// U, -> uncertain, this loop-size with this condition has no definite canonical chiral reference in the first place!

	Column turn_type ("turn_type", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);

	Column i_AA ("i_AA", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column i_p1_AA ("i_p1_AA", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column i_p2_AA ("i_p2_AA", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column i_p3_AA ("i_p3_AA", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);

	Column canonical_turn_AA ("canonical_turn_AA", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);

	Column PA_by_preceding_E ("PA_by_preceding_E", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column PA_by_following_E ("PA_by_following_E", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);

	Column cano_PA ("cano_PA", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	// T, -> true, canonical PA
	// F, -> false, non-canonical PA
	// U, -> uncertain, this loop-size with this condition has no definite canonical PA reference in the first place!

	Column heading_direction ("heading_direction", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column parallel_EE ("parallel_EE", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column cano_parallel_EE ("cano_parallel_EE", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column component_size ("component_size", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column A  ("A", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column C  ("C", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column D  ("D", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column E  ("E", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column F  ("F", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column G  ("G", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column H  ("H", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column I  ("I", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column K  ("K", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column L  ("L", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column M  ("M", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column N  ("N", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column P  ("P", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column Q  ("Q", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column R  ("R", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column S  ("S", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column T  ("T", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column V  ("V", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column W  ("W", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column Y  ("Y", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column A_core_heading  ("A_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column A_surface_heading  ("A_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column C_core_heading  ("C_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column C_surface_heading  ("C_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);


	Column D_core_heading  ("D_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column D_surface_heading  ("D_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column E_core_heading  ("E_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column E_surface_heading  ("E_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);


	Column F_core_heading  ("F_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column F_surface_heading  ("F_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column G_core_heading  ("G_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column G_surface_heading  ("G_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column H_core_heading  ("H_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column H_surface_heading  ("H_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column I_core_heading  ("I_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column I_surface_heading  ("I_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column K_core_heading  ("K_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column K_surface_heading  ("K_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column L_core_heading  ("L_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column L_surface_heading  ("L_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column M_core_heading  ("M_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column M_surface_heading  ("M_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column N_core_heading  ("N_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column N_surface_heading  ("N_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);


	Column P_core_heading  ("P_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column P_surface_heading  ("P_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column Q_core_heading  ("Q_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column Q_surface_heading  ("Q_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column R_core_heading  ("R_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column R_surface_heading  ("R_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);


	Column S_core_heading  ("S_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column S_surface_heading  ("S_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column T_core_heading  ("T_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column T_surface_heading  ("T_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);


	Column V_core_heading  ("V_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column V_surface_heading  ("V_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column W_core_heading  ("W_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column W_surface_heading  ("W_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column Y_core_heading  ("Y_core_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);
	Column Y_surface_heading  ("Y_surface_heading", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);


	Column H_percentage ("H_percentage", DbDataTypeOP( new DbReal() ), true /*could be null*/, false /*don't autoincrement*/);
	Column E_percentage ("E_percentage", DbDataTypeOP( new DbReal() ), true /*could be null*/, false /*don't autoincrement*/);
	Column L_percentage ("L_percentage", DbDataTypeOP( new DbReal() ), true /*could be null*/, false /*don't autoincrement*/);

	Column number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands ("number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands ("number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands", DbDataTypeOP( new DbInteger() ), true /*could be null*/, false /*don't autoincrement*/);

	Column residue_begin("residue_begin", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);
	Column residue_end  ("residue_end", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*don't autoincrement*/);


	Column number_of_hydrophobic_res ("number_of_hydrophobic_res", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_hydrophilic_res ("number_of_hydrophilic_res", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column number_of_CGP ("number_of_CGP", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column ratio_hydrophobic_philic_of_sw_in_percent ("ratio_hydrophobic_philic_of_sw_in_percent", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);

	Column number_of_RK_in_sw ("number_of_RK_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_DE_in_sw ("number_of_DE_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column net_charge_of_sw ("net_charge_of_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column number_of_core_heading_FWY_in_sw ("number_of_core_heading_FWY_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column ratio_of_core_heading_FWY_in_sw ("ratio_of_core_heading_FWY_in_sw", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_core_heading_W_in_sw ("number_of_core_heading_W_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_core_heading_L_in_core_strands_in_sw ("number_of_core_heading_L_in_core_strands_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_core_heading_W_in_core_strands_in_sw ("number_of_core_heading_W_in_core_strands_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column number_of_core_heading_Y_in_core_strands_in_sw ("number_of_core_heading_Y_in_core_strands_in_sw", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column avg_dihedral_angle_between_core_strands_across_facing_sheets ("avg_dihedral_angle_between_core_strands_across_facing_sheets", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);
	Column sw_res_size ("sw_res_size", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);

	Column multimer_is_suspected ("multimer_is_suspected", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	Column avg_b_factor_CB_at_each_component ("avg_b_factor_CB_at_each_component", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);

	Column topology_candidate ("topology_candidate", DbDataTypeOP( new DbText() ), true /* could be null*/, false /*no autoincrement*/);
	// "fnIII" or "not_fnIII" or "not_known_topology"

	Column min_dis_between_sheets_by_all_res ("min_dis_between_sheets_by_all_res", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);
	// all residues of not only edge strands, but all strands

	Column min_dis_between_sheets_by_cen_res ("min_dis_between_sheets_by_cen_res", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);
	// central residues of not only edge strands, but all strands

	Column avg_dis_between_sheets_by_cen_res ("avg_dis_between_sheets_by_cen_res", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);
	// central residues of not only edge strands, but all strands

	Column num_PRO_in_starting_loop_and_1st_3rd_inter_sheet_loop ("num_PRO_in_starting_loop_and_1st_3rd_inter_sheet_loop", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column num_PRO_in_starting_loop ("num_PRO_in_starting_loop", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column num_PRO_in_1st_inter_sheet_loop ("num_PRO_in_1st_inter_sheet_loop", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column num_PRO_in_3rd_inter_sheet_loop ("num_PRO_in_3rd_inter_sheet_loop", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column weighted_num_PRO_prevent ("weighted_num_PRO_prevent", DbDataTypeOP( new DbInteger() ), true /* could be null*/, false /*no autoincrement*/);
	Column shortest_dis_between_facing_aro_in_sw ("shortest_dis_between_facing_aro_in_sw", DbDataTypeOP( new DbReal() ), true /* could be null*/, false /*no autoincrement*/);

	// Schema
	// PrimaryKey
	utility::vector1<Column> primary_key_columns_sandwich;
	primary_key_columns_sandwich.push_back(struct_id);
	primary_key_columns_sandwich.push_back(sandwich_PK_id);

	// table name
	Schema sandwich("sandwich",  PrimaryKey(primary_key_columns_sandwich));

	// add column which is neither PrimaryKey nor ForeignKey
	sandwich.add_column(tag);
	sandwich.add_column(sw_can_by_sh_id);
	sandwich.add_column(sheet_id);

	sandwich.add_column(sheet_antiparallel);
	sandwich.add_column(sandwich_bs_id);
	sandwich.add_column(long_strand_id);
	sandwich.add_column(strand_edge);
	sandwich.add_column(num_strands_in_each_sw);
	sandwich.add_column(num_edge_strands_in_each_sw);

	sandwich.add_column(intra_sheet_con_id);
	sandwich.add_column(inter_sheet_con_id);
	sandwich.add_column(loop_kind); // better to be located right after intra_sheet_con_id/inter_sheet_con_id for better readability

	sandwich.add_column(LR);
	sandwich.add_column(canonical_LR);

	sandwich.add_column(turn_type);
	sandwich.add_column(i_AA);
	sandwich.add_column(i_p1_AA);
	sandwich.add_column(i_p2_AA);
	sandwich.add_column(i_p3_AA);

	sandwich.add_column(canonical_turn_AA);

	sandwich.add_column(PA_by_preceding_E);
	sandwich.add_column(PA_by_following_E);
	sandwich.add_column(cano_PA);
	sandwich.add_column(heading_direction); // nega,posi,meet ...
	sandwich.add_column(parallel_EE);
	sandwich.add_column(cano_parallel_EE);

	sandwich.add_column(A);
	sandwich.add_column(C);
	sandwich.add_column(D);
	sandwich.add_column(E);
	sandwich.add_column(F);
	sandwich.add_column(G);
	sandwich.add_column(H);
	sandwich.add_column(I);
	sandwich.add_column(K);
	sandwich.add_column(L);
	sandwich.add_column(M);
	sandwich.add_column(N);

	sandwich.add_column(P);

	sandwich.add_column(Q);

	sandwich.add_column(R);
	sandwich.add_column(S);
	sandwich.add_column(T);

	sandwich.add_column(V);

	sandwich.add_column(W);

	sandwich.add_column(Y);

	sandwich.add_column(A_core_heading);
	sandwich.add_column(A_surface_heading);
	sandwich.add_column(C_core_heading);
	sandwich.add_column(C_surface_heading);


	sandwich.add_column(D_core_heading);
	sandwich.add_column(D_surface_heading);
	sandwich.add_column(E_core_heading);
	sandwich.add_column(E_surface_heading);


	sandwich.add_column(F_core_heading);
	sandwich.add_column(F_surface_heading);

	sandwich.add_column(G_core_heading);
	sandwich.add_column(G_surface_heading);


	sandwich.add_column(H_core_heading);
	sandwich.add_column(H_surface_heading);
	sandwich.add_column(I_core_heading);
	sandwich.add_column(I_surface_heading);

	sandwich.add_column(K_core_heading);
	sandwich.add_column(K_surface_heading);

	sandwich.add_column(L_core_heading);
	sandwich.add_column(L_surface_heading);

	sandwich.add_column(M_core_heading);
	sandwich.add_column(M_surface_heading);

	sandwich.add_column(N_core_heading);
	sandwich.add_column(N_surface_heading);

	sandwich.add_column(P_core_heading);
	sandwich.add_column(P_surface_heading);

	sandwich.add_column(Q_core_heading);
	sandwich.add_column(Q_surface_heading);

	sandwich.add_column(R_core_heading);
	sandwich.add_column(R_surface_heading);

	sandwich.add_column(S_core_heading);
	sandwich.add_column(S_surface_heading);
	sandwich.add_column(T_core_heading);
	sandwich.add_column(T_surface_heading);

	sandwich.add_column(V_core_heading);
	sandwich.add_column(V_surface_heading);
	sandwich.add_column(W_core_heading);
	sandwich.add_column(W_surface_heading);

	sandwich.add_column(Y_core_heading);
	sandwich.add_column(Y_surface_heading);


	sandwich.add_column(H_percentage);
	sandwich.add_column(E_percentage);
	sandwich.add_column(L_percentage);

	sandwich.add_column(number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands);
	sandwich.add_column(number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands);

	sandwich.add_column(number_of_hydrophobic_res); // A,V,I,L,M,F,Y,W
	sandwich.add_column(number_of_hydrophilic_res); // R,H,K,D,E,S,T,N,Q
	sandwich.add_column(number_of_CGP); // C,G,P
	sandwich.add_column(ratio_hydrophobic_philic_of_sw_in_percent); // (no_hydrophobic/no_hydrophilic)*100

	sandwich.add_column(number_of_RK_in_sw); // R,K
	sandwich.add_column(number_of_DE_in_sw); // D,E
	sandwich.add_column(net_charge_of_sw);
	sandwich.add_column(number_of_core_heading_FWY_in_sw);
	sandwich.add_column(ratio_of_core_heading_FWY_in_sw);
	sandwich.add_column(number_of_core_heading_W_in_sw);
	sandwich.add_column(number_of_core_heading_L_in_core_strands_in_sw);
	sandwich.add_column(number_of_core_heading_W_in_core_strands_in_sw);
	sandwich.add_column(number_of_core_heading_Y_in_core_strands_in_sw);
	sandwich.add_column(avg_dihedral_angle_between_core_strands_across_facing_sheets);
	sandwich.add_column(sw_res_size);
	sandwich.add_column(multimer_is_suspected);
	sandwich.add_column(avg_b_factor_CB_at_each_component);
	sandwich.add_column(topology_candidate);
	sandwich.add_column(min_dis_between_sheets_by_all_res);
	sandwich.add_column(min_dis_between_sheets_by_cen_res);
	sandwich.add_column(avg_dis_between_sheets_by_cen_res);
	sandwich.add_column(shortest_dis_between_facing_aro_in_sw);

	sandwich.add_column(num_PRO_in_starting_loop_and_1st_3rd_inter_sheet_loop);
	sandwich.add_column(num_PRO_in_starting_loop);
	sandwich.add_column(num_PRO_in_1st_inter_sheet_loop);
	sandwich.add_column(num_PRO_in_3rd_inter_sheet_loop);
	sandwich.add_column(weighted_num_PRO_prevent);

	sandwich.add_column(component_size);


	// ForeignKey
	sandwich.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	utility::vector1<Column> residue_begin_fkey_cols;
	residue_begin_fkey_cols.push_back(struct_id);
	residue_begin_fkey_cols.push_back(residue_begin);

	sandwich.add_foreign_key(ForeignKey(residue_begin_fkey_cols, "residues", fkey_reference_cols, true /*defer*/));

	utility::vector1<Column> residue_end_fkey_cols;
	residue_end_fkey_cols.push_back(struct_id);
	residue_end_fkey_cols.push_back(residue_end);

	sandwich.add_foreign_key(ForeignKey(residue_end_fkey_cols, "residues", fkey_reference_cols, true /*defer*/));

	sandwich.write(db_session);
	/****** <end> writing sandwich ******/


	/****** <begin> writing rkde_in_strands ******/

	// Columns

	// unique_primary_key
	Column rkde_in_strands_PK_id ("rkde_in_strands_PK_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

	// may not be unique
	Column residue_number ("residue_number", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);
	Column residue_type ("residue_type", DbDataTypeOP( new DbText() ), false /*not null*/, false /*no autoincrement*/);

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
	rkde_in_strands.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	rkde_in_strands.write(db_session);
	/****** <end> rkde_in_strands ******/


	/****** <begin> writing rkde ******/

	// Columns

	// unique_primary_key
	Column rkde_PK_id ("rkde_PK_id", DbDataTypeOP( new DbInteger() ), false /*not null*/, false /*no autoincrement*/);

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
	rkde.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

	rkde.write(db_session);
	/****** <end> rkde ******/

}

core::scoring::ScoreFunctionOP
SandwichFeatures::generate_scorefxn( bool fullatom ) {
	core::scoring::ScoreFunctionOP scorefxn( nullptr );
	if ( fullatom ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	}
	return scorefxn;
}


void
SandwichFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
)
{
	if ( tag->hasOption( "min_num_strands_to_deal" ) ) {
		min_num_strands_to_deal_ = tag->getOption<Size>("min_num_strands_to_deal", 4); //needed for development or custom use
	} else {
		min_num_strands_to_deal_ = 4; // needed for unit_test or default use
	}
	// At least 4 strands should be in pdb file


	if ( tag->hasOption( "max_num_strands_to_deal" ) ) {
		max_num_strands_to_deal_ = tag->getOption<Size>("max_num_strands_to_deal", 140); //needed for development or custom use
	} else {
		max_num_strands_to_deal_ = 140; // needed for unit_test or default use
	}
	// example: (in all chains of 1FE8) There are 132 strands, it took ~ 7 cpu minutes to process



	if ( tag->hasOption( "min_res_in_strand" ) ) {
		min_res_in_strand_ = tag->getOption<Size>("min_res_in_strand", 2); //needed for development or custom use
	} else {
		min_res_in_strand_ = 2; // needed for unit_test or default use
	}
	// definition: minimum number of residues in a strand, for edge strand definition & analysis
	// example: 4=< is recommended (in 1A8M) min_res_in_strand = 2, (in 1PMY) min_res_in_strand = 3


	if ( tag->hasOption( "min_CA_CA_dis" ) ) {
		min_CA_CA_dis_ = tag->getOption<Real>("min_CA_CA_dis", 3.5); //needed for development or custom use
	} else {
		min_CA_CA_dis_ = 3.5; // needed for unit_test or default use
	}
	// definition: minimum CA_CA_distance between strands in same sheet
	// example: (in 1A8M) 'min_CA_CA_dis_= 3.5', (in 1KIT) 'min_CA_CA_dis_= 4.0'


	if ( tag->hasOption( "max_CA_CA_dis" ) ) {
		max_CA_CA_dis_ = tag->getOption<Real>("max_CA_CA_dis", 6.2); //needed for development or custom use
	} else {
		max_CA_CA_dis_ = 6.2; // needed for unit_test or default use
	}
	// example: (in 1A8M) 'max_CA_CA_dis_= 6.2', (in 1KIT) 'max_CA_CA_dis_= 5.7'



	min_N_O_dis_between_two_sheets_ = tag->getOption<Real>("min_N_O_dis_between_sheets", 3.3);
	// definition: min distance between bb N and bb O between two sheets
	// example: (in c.10.0, 4.7 Angstrom exists between O and N of edge strands of two sheets) this is a canonical sw
	// example: (in c.128.0, 3.1 Angstrom exists between O and N of edge strands of two sheets) this is not a canonical sw

	min_N_H_O_angle_between_two_sheets_ = tag->getOption<Real>("min_N_H_O_angle_between_two_sheets", 154.0);
	// definition: minimum N-H-O angle between bb N and bb O between two sheets
	// example: (in c.30.0, 147.2 degree angle exists between O and N of edge strands of two sheets) this is a canonical sw
	// example: (in c.26.0, 155.5 degree angle exists between O and N of edge strands of two sheets) this is not a canonical sw
	// example: (in c.128.0, 156.4 degree angle exists between O and N of edge strands of two sheets) this is not a canonical sw

	min_C_O_N_angle_ = tag->getOption<Real>("min_C_O_N_angle", 120.0);
	// example: (in 1L0Q chain A), 138 is the smallest C_O_N_angle (C and O from one sheet, N from other sheet)
	min_sheet_dis_ = tag->getOption<Real>("min_sheet_dis", 7.0);
	// definition: minimum CA_CA_distance between strands in different sheets to constitute a sandwich
	// 7 Angstrom seems OK
	max_sheet_dis_ = tag->getOption<Real>("max_sheet_dis", 15.0);
	// definition: maximum CA_CA_distance between strands in different sheets to constitute a sandwich
	// 15 Angstrom seems OK
	max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_ = tag->getOption<Real>("max_sheet_angle_with_cen_res", 130.0);
	// definition: Maximum angle between sheets (CA and CA) with two terminal residues and one central residue
	// example: I need to remove non-facing sheets

	min_sheet_angle_by_four_term_cen_res_ = tag->getOption<Real>("min_sheet_angle_by_four_term_cen_res", 25.0);
	// definition: Minimum angle between sheets (CA and CA)
	// defined by: 3 middle residues in each edge strand among 4 middle residues
	// usage: used in judge_facing "angle_1 > min_sheet_angle_by_four_term_cen_res_"
	// [example] during repopulated de novo design of 3_1L9N, 27.2 of min_sheet_angle was observed for ideal-looking decoy

	max_sheet_angle_by_four_term_cen_res_ = tag->getOption<Real>("max_sheet_angle_by_four_term_cen_res", 150.0);
	// definition: Maximum angle between sheets (CA and CA)
	// In [1TEN] even 155 degree comes from same sheet!

	min_sheet_torsion_cen_res_ = tag->getOption<Real>("min_sheet_torsion_cen_res", -150.0);
	// definition: "Minimum torsion between sheets (CA and CA) with respect to terminal central residues in each beta-sheet
	// explanation: with respect to central residues, one torsion angles of 1TEN is 84.9
	max_sheet_torsion_cen_res_ = tag->getOption<Real>("max_sheet_torsion_cen_res", 150.0);
	// definition: "maximum torsion between sheets (CA and CA) with respect to terminal central residues in each beta-sheet
	// usage: used in judge_facing "torsion_i_j < max_sheet_torsion_cen_res_"
	min_num_strands_in_sheet_ = tag->getOption<core::Size>("min_num_strands_in_sheet", 2);
	//  definition: a sheet with < 3 strands will be ignored
	// usage: if (num_strands_i < min_num_strands_in_sheet_)
	min_inter_sheet_dis_CA_CA_ = tag->getOption<Real>("min_inter_sheet_dis_CA_CA", 4.0);
	// example: (in 12E8) the distance between S52 and G64 is 4.2 A
	max_inter_sheet_dis_CA_CA_ = tag->getOption<Real>("max_inter_sheet_dis_CA_CA", 24.0);
	// example: (in 2WYF) the distance between val and gln is 5.8 A
	// usage: shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_
	// example: (in 1TEN) shortest_avg_dis_inter_sheet between sheet 1 and 2 = 11.6 A (and these sheets should be a sandwich)
	// example: (in 1A64_chain_A) shortest_avg_dis_inter_sheet between sheet 1 and 2 = 25 A (and these sheets should not be a sandwich)
	// example: (in 1A1M) the average distance between sheet 1 and 4 > 20 A (but these sheets should be a sandwich)
	// example: (in 1ASQ) the average distance between sheet 1 and 4 > 20 A (and these sheets should not be a sandwich)

	max_inter_strand_angle_to_not_be_same_direction_strands_ = tag->getOption<Real>("max_inter_strand_angle_to_not_be_same_direction_strands", 120.0);
	// usage:  if (angle_start_res_being_middle > max_inter_strand_angle_to_not_be_same_direction_strands_)
	// example: (in 1BQB chain A) 121 is possible (but this should be excluded as same direction strand)
	// example: (in 1A0Q chain L) 127 is possible for 3-6-9 angle (but this should be excluded as same direction strand)
	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_ = tag->getOption<Real>("max_abs_inter_strand_dihedral_to_not_be_same_direction_strands", 100.0);
	// usage: if ( std::abs(torsion_between_strands) > max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_)
	// example: (in 1U3J chain A) 105 is possible for 4-7-11-13 dihedral angle (but this should be excluded as same direction strand)
	// example: (in 1QAC chain A) 128.5 is possible for 4-7-10-13 dihedral angle (but this should be excluded as same direction strand)
	// example: (in 1A3R chain L) 130 is possible for 4-7-10-14 dihedral angle (but this should be excluded as same direction strand)
	max_num_sw_per_pdb_ = tag->getOption<core::Size>("max_num_sw_per_pdb", 100);
	// definition: maximum number of sandwiches to be extracted per a pdb file
	check_N_to_C_direction_by_ = tag->getOption<string>("check_N_to_C_direction_by", "PE");
	// definition: check N->C going direction by option
	// PE: preceding residue in beta-strand
	// FE: following residue in
	check_canonicalness_cutoff_ = tag->getOption<Real>("check_canonicalness_cutoff", 80.0);
	// definition: cutoff to determine canonicalness of L/R, P/A and directionality


	///////// An option that takes the longest time ///////
	count_AA_with_direction_ = tag->getOption<bool>("count_AA_with_direction", false);
	// definition: if true, count AA considering direction too!
	// <note> if true, it takes more time, but ~50 sandwiches this can be used within ~ minutes.


	///////// strictness options ///////
	do_not_write_resfile_of_sandwich_that_has_non_canonical_LR_ = tag->getOption<bool>("do_not_write_resfile_of_sandwich_that_has_non_canonical_LR", false);
	// definition: if true, exclude sandwich_that_has_non_canonical_LR

	exclude_desinated_pdbs_ = tag->getOption<bool>("exclude_desinated_pdbs", false);
	// definition: if true, exclude certain designated pdbs

	exclude_sandwich_that_has_near_backbone_atoms_between_sheets_ = tag->getOption<bool>("exclude_sandwich_that_has_near_backbone_atoms_between_sheets", false);
	// definition: if true, exclude sandwich_that_has_near_backbone_atoms_between_sheets

	exclude_sandwich_that_has_non_canonical_LR_ = tag->getOption<bool>("exclude_sandwich_that_has_non_canonical_LR", false);
	// definition: if true, exclude sandwich_that_has_non_canonical_LR

	exclude_sandwich_that_has_non_canonical_properties_ = tag->getOption<bool>("exclude_sandwich_that_has_non_canonical_properties", false);

	exclude_sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw_ = tag->getOption<bool>("exclude_sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw", false);
	// definition: if true, exclude sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw_

	exclude_sandwich_that_is_linked_w_same_direction_strand_ = tag->getOption<bool>("exclude_sandwich_that_is_linked_w_same_direction_strand", false);
	// definition: if true, exclude a sandwich that is linked with same_direction_strand
	// Rationale of default=true (1)
	// If true, it is useful to exclude [1QAC]_chain_A, [2v33]_chain_A which is a canonical sandwich but linked by same direction strands between sheets

	exclude_sandwich_that_is_suspected_to_have_not_facing_2_sheets_ = tag->getOption<bool>("exclude_sandwich_that_is_suspected_to_have_not_facing_2_sheets", true);
	// definition: if true, exclude that_is_suspected_to_have_not_facing_2_sheets

	exclude_sandwich_with_SS_bond_ = tag->getOption<bool>("exclude_sandwich_with_SS_bond", false);


	max_starting_loop_size_ = tag->getOption<core::Size>("max_starting_loop_size", 6);
	// definition: maximum starting loop size to extract
	max_ending_loop_size_ = tag->getOption<core::Size>("max_ending_loop_size", 6);
	// definition: maximum ending loop size to extract
	max_E_in_extracted_sw_loop_ = tag->getOption<core::Size>("max_E_in_extracted_sw_loop", 10);
	// definition: maximum allowable number of E residues in extracted sandwich loop
	// usefulness: If used, it is useful to exclude [1LOQ] which is a beta-propeller

	//max_H_in_extracted_sw_loop_ = tag->getOption<core::Size>("max_H_in_extracted_sw_loop", 9);
	max_H_in_extracted_sw_loop_ = tag->getOption<core::Size>("max_H_in_extracted_sw_loop", 100);
	// definition: maximum allowable number of helix residues in extracted sandwich loop
	// example: 0 would be ideal, but then only ~10% of sandwiches will be extracted among CATH classified sandwiches instead even when same_direction_strand linking sw is allowed!

	no_helix_in_pdb_ = tag->getOption<bool>("no_helix_in_pdb", false);
	// if true, ignore any pdb that has helix

	inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_ = tag->getOption<Real>("inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets", 13.0);
	// definition: within this distance, sheets are considered to be too near each other
	// example: (in 1LOQ) inter-sheet distances are 11.5~14.1
	// Rationale of default value=13 Angstron
	// it is useful to exclude [1LOQ] which is beta-propeller and [3BVT] which is a stacked sandwich
	// but it also excludes [2V33] which has two canonical sandwiches near each other and [1W8O] which is a canonical sandwich near a single beta-sheet

	allowed_deviation_for_turn_type_id_ = tag->getOption<Real>("allowed_deviation_for_turn_type_id", 40.0);


	primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping_ = tag->getOption<int>("primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping", 2);

	primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping_ = tag->getOption<int>("primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping", 0);

	primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping_ = tag->getOption<int>("primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping", 0);

	primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping_ = tag->getOption<int>("primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping", 2);


	/// electrostatic_interactions related
	distance_cutoff_for_electrostatic_interactions_ = tag->getOption<Real>("distance_cutoff_for_electrostatic_interactions", 7.0);
	//"N-O bridges follow only one, that is, they have at least a pair of side-chain functional-group nitrogen and oxygen atoms within 4A distance, but the side-chain functional-group centroids are > 4A apart."
	// source: 2002_Close-Range Electrostatic Interactions in Proteins

	CB_b_factor_cutoff_for_electrostatic_interactions_ = tag->getOption<Real>("CB_b_factor_cutoff_for_electrostatic_interactions", 10000);
	// CB_b_factor_cutoff_for_electrostatic_interactions_ = tag->getOption<Real>("CB_b_factor_cutoff_for_electrostatic_interactions", 60);
	//"Values of 60 or greater may imply disorder (for example, free movement of a side chain or alternative side-chain conformations). Values of 20 and 5 correspond to uncertainties of 0.5 and 0.25 angstroms, respectively."
	// source: http://spdbv.vital-it.ch/TheMolecularLevel/SPVTut/text/STut09aTN.html

	min_primary_seq_distance_diff_for_electrostatic_interactions_ = tag->getOption<core::Size>("min_primary_seq_distance_diff_for_electrostatic_interactions", 4);
	// rationale for default value: I hypothesize that electrostatic interaction between 38E and 41K of Tencon do not stabilize that much


	///////// development options ///////
	do_not_connect_sheets_by_loops_ = tag->getOption<bool>("do_not_connect_sheets_by_loops", false);
	// definition: if true, don't connect sheets by loops
	extract_sandwich_ = tag->getOption<bool>("extract_sandwich", true);


	///////// proline options ///////
	wt_for_pro_in_starting_loop_ = tag->getOption<Real>("wt_for_pro_in_starting_loop", 20);
	wt_for_pro_in_1st_inter_sheet_loop_ = tag->getOption<Real>("wt_for_pro_in_1st_inter_sheet_loop", 10);
	wt_for_pro_in_3rd_inter_sheet_loop_ = tag->getOption<Real>("wt_for_pro_in_3rd_inter_sheet_loop", 1);


	///////// writing options ///////
	write_all_info_files_ = tag->getOption<bool>("write_all_info_files", false);
	// definition: if true, write all below

	write_AA_kind_files_ = tag->getOption<bool>("write_AA_kind_files", false);
	// definition: if true, write files that have amino acid kinds


	///////// AA distribution ///////
	write_loop_AA_distribution_files_ = tag->getOption<bool>("write_loop_AA_distribution_files", false);
	// definition: if true, write files that have amino acid distributions without directions

	write_strand_AA_distribution_files_ = tag->getOption<bool>("write_strand_AA_distribution_files", false);
	// definition: if true, write files that have amino acid distributions with directions


	write_beta_sheet_capping_info_ = tag->getOption<bool>("write_beta_sheet_capping_info", false);
	// reference: 2008_beta-Sheet capping- Signals that initiate and terminate beta-sheet formation, Journal of Structural Biology FarzadFard et al.,

	write_chain_B_resnum_ = tag->getOption<bool>("write_chain_B_resnum", false);
	// if true, write chain_B_resnum file for InterfaceAnalyzer

	write_electrostatic_interactions_of_all_residues_ = tag->getOption<bool>("write_electrostatic_interactions_of_all_residues", false);
	write_electrostatic_interactions_of_all_residues_in_a_strand_ = tag->getOption<bool>("write_electrostatic_interactions_of_all_residues_in_a_strand", false);
	write_electrostatic_interactions_of_surface_residues_in_a_strand_ = tag->getOption<bool>("write_electrostatic_interactions_of_surface_residues_in_a_strand", false);

	write_heading_directions_of_all_AA_in_a_strand_ = tag->getOption<bool>("write_heading_directions_of_all_AA_in_a_strand", false);


	write_p_aa_pp_files_ = tag->getOption<bool>("write_p_aa_pp_files", false);
	// definition: if true, write p_aa_pp_files
	write_phi_psi_of_all_ = tag->getOption<bool>("write_phi_psi_of_all", false);
	// definition: if true, write phi_psi_file
	write_phi_psi_of_E_ = tag->getOption<bool>("write_phi_psi_of_E", false);
	// definition: if true, write phi_psi_file
	write_rama_at_AA_to_files_ = tag->getOption<bool>("write_rama_at_AA_to_files", false);
	// definition: if true, write write_rama_at_AA_to_files

	write_resfile_ = tag->getOption<bool>("write_resfile", false);
	// definition: if true, write resfile automatically

	write_resfile_NOT_FWY_on_surface_ = tag->getOption<bool>("write_resfile_NOT_FWY_on_surface", false);

	write_resfile_when_seq_rec_is_bad_ = tag->getOption<bool>("write_resfile_when_seq_rec_is_bad", false);
	//      definition: if true, write resfile as if score fn's seq_rec_is_bad

	write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_ = tag->getOption<bool>("write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands", false);
	// definition: if true, write resfile to_minimize_too_many_core_heading_FWY_on_core_strands

	write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_ = tag->getOption<bool>("write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands", false);
	// definition: if true, write resfile to_minimize_too_many_core_heading_FWY_on_edge_strands

	write_resfile_to_minimize_too_much_hydrophobic_surface_ = tag->getOption<bool>("write_resfile_to_minimize_too_much_hydrophobic_surface", false);
	// definition: if true, write resfile to_minimize_too_much_hydrophobic_surface


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////    SandwichFeatures   ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief collect all the feature data for the pose
core::Size
SandwichFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	try // added try-catch according to rosetta/main/source/stubs/Application.cc
{
		TR.Info << "======================= <begin> report_features of SandwichFeatures =========================" << endl;

		unit_test_pass_identifier = true; // assumed to pass unit_test now

		string tag = get_tag(struct_id, db_session);
		bool canonical_sw_extracted_from_this_pdb_file = false;

		// <begin> exclude_desinated_pdbs
		if ( exclude_desinated_pdbs_ ) {
			bool this_pdb_should_be_excluded = check_whether_this_pdb_should_be_excluded(tag);
			if ( this_pdb_should_be_excluded ) {
				TR.Info << "Exit early since this pdb should be excluded " << endl;
				return 0;
			}
		}
		// <end> exclude_desinated_pdbs

		// utility::vector1< pose::PoseOP > singlechain_poses;
		// singlechain_poses = pose.split_by_chain();
		//  TR << "singlechain_poses.size(): " << singlechain_poses.size() << endl;
		//
		// for(Size pose_i=1; pose_i<=singlechain_poses.size(); pose_i++)
		// {
		//  pose::Pose pose ( singlechain_poses[ pose_i ] );


		pose::Pose dssp_pose ( pose ); //copy of pose, since the original pose is called as 'const'
		core::scoring::dssp::Dssp dssp( dssp_pose );
		dssp.insert_ss_into_pose( dssp_pose );


		if ( no_helix_in_pdb_ ) {
			bool helix_existence = check_helix_existence(dssp_pose);
			if ( helix_existence ) {
				TR.Info << "Exit early since this pdb has helix " << endl;
				return 0;
			}
		}

		if ( write_all_info_files_ ) {
			write_AA_kind_files_ = true;
			write_beta_sheet_capping_info_ = true;

			write_strand_AA_distribution_files_ = true;
			write_loop_AA_distribution_files_ = true;

			write_chain_B_resnum_ = true;
			write_electrostatic_interactions_of_surface_residues_in_a_strand_ = true;
			write_electrostatic_interactions_of_all_residues_in_a_strand_ = true;
			write_electrostatic_interactions_of_all_residues_ = true;
			write_heading_directions_of_all_AA_in_a_strand_ = true;

			write_p_aa_pp_files_ = true;
			write_phi_psi_of_all_ = true;
			write_phi_psi_of_E_ = true;

			write_rama_at_AA_to_files_ = true;
			write_resfile_ = true;
			write_resfile_to_minimize_too_much_hydrophobic_surface_ = true;
			write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_ = true;
			write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_ = true;
			write_resfile_when_seq_rec_is_bad_ = true;
		}

		if ( exclude_sandwich_that_has_non_canonical_properties_ ) {
			//exclude_sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw_ = true;
			exclude_sandwich_that_has_non_canonical_LR_ = true;
		}

		pose::Pose pose_w_center_000 = pose;
		if ( count_AA_with_direction_ || write_resfile_ ) {
			pose_w_center_000.center();
		}


		Size sheet_PK_id_counter=1; //initial value
		Size sw_can_by_sh_PK_id_counter=1; //initial value

		Size rkde_in_strands_PK_id_counter=1; //initial value
		Size rkde_PK_id_counter=1; //initial value

		Size sw_can_by_sh_id_counter=1; //initial value
		Size sandwich_PK_id_counter=1; //initial value
		Size intra_sheet_con_id_counter=1; //initial value
		Size inter_sheet_con_id_counter=1; //initial value

		utility::vector1<SandwichFragment> all_strands = get_full_strands(struct_id, db_session);


		if ( (static_cast<core::Size>(all_strands.size()) < min_num_strands_to_deal_) || (static_cast<core::Size>(all_strands.size()) > max_num_strands_to_deal_) ) {

			TR.Info << "Exit early since all_strands.size(): " << all_strands.size() << endl;
			return 0;
		}

		// <begin> assignment of strands into the sheet & define very first sheet ("1")
		bool first_sheet_assigned = false;
		for ( Size i=1; i<all_strands.size() && !first_sheet_assigned; ++i ) { // I don't need the last strand since this double for loops are exhaustive search for all pairs of strands
			if ( all_strands[i].get_size() < min_res_in_strand_ ) {
				continue;
			}
			for ( Size j=i+1; j<=all_strands.size() && !first_sheet_assigned; ++j ) { // I need the last strand for this second for loop
				if ( all_strands[j].get_size() < min_res_in_strand_ ) {
					continue;
				}
				SandwichFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
				SandwichFragment temp_strand_j(all_strands[j].get_start(), all_strands[j].get_end());

				// <begin> anti-parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand
				Size return_of_find_sheet_antiparallel(0); // temporary 'false' designation
				Size return_of_find_sheet_parallel(0); // temporary 'false' designation
				return_of_find_sheet_antiparallel =
					find_sheet (
					pose, temp_strand_i, temp_strand_j, true,
					min_CA_CA_dis_,
					max_CA_CA_dis_,
					min_C_O_N_angle_,
					false); //care_smaller_sheet
				if ( return_of_find_sheet_antiparallel == 999 ) {
					break;
				}
				if ( !return_of_find_sheet_antiparallel ) {
					return_of_find_sheet_parallel =
						find_sheet (
						pose, temp_strand_i, temp_strand_j, false,
						min_CA_CA_dis_,
						max_CA_CA_dis_,
						min_C_O_N_angle_,
						true); //care_smaller_sheet for parallel search
				}
				if ( return_of_find_sheet_parallel == 999 ) { // since these two strands are too distant to each other, there is virtually no chance to be sheet!
					break;
				}
				if ( return_of_find_sheet_antiparallel || return_of_find_sheet_parallel ) {
					first_sheet_assigned = true;
					WriteToDB_sheet (
						struct_id,
						db_session,
						sheet_PK_id_counter,
						1, // sheet_id
						get_segment_id( // segment_id
						struct_id,
						db_session,
						i));

					sheet_PK_id_counter++;
					WriteToDB_sheet (
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
		for ( Size i=1; i<=all_strands.size(); ++i ) {
			if ( all_strands[i].get_size() >= min_res_in_strand_ ) { // the length of this beta strand > min_res_in_strand_
				bool strand_i_is_in_any_sheet = check_whether_strand_i_is_in_sheet(struct_id, db_session, get_segment_id(
					struct_id,
					db_session,
					i));
				if ( !strand_i_is_in_any_sheet ) { //  this strand is not in any sheet, so this strand needs to be assigned into any sheet
					utility::vector1<SandwichFragment> cur_strands = get_current_strands_in_sheet(struct_id, db_session);
					bool sheet_unassigned = true;
					for ( Size j=1; j<=cur_strands.size() && sheet_unassigned == true; ++j ) {
						SandwichFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
						SandwichFragment temp_strand_j(cur_strands[j].get_start(), cur_strands[j].get_end());

						Size return_of_find_sheet_antiparallel(0); // temporary 'false' designation
						Size return_of_find_sheet_parallel(0); // temporary 'false' designation

						return_of_find_sheet_antiparallel =
							find_sheet (
							pose, temp_strand_i, temp_strand_j, true,
							min_CA_CA_dis_,
							max_CA_CA_dis_,
							min_C_O_N_angle_,
							false); //care_smaller_sheet

						if ( return_of_find_sheet_antiparallel == 999 ) {
							break; // too distant strands
						}

						if ( !return_of_find_sheet_antiparallel ) {
							return_of_find_sheet_parallel =
								find_sheet (
								pose, temp_strand_i, temp_strand_j, false,
								min_CA_CA_dis_,
								max_CA_CA_dis_,
								min_C_O_N_angle_,
								true); //care_smaller_sheet for parallel finding
						}

						if ( return_of_find_sheet_parallel == 999 ) {
							break; // too distant strands
						}

						if ( return_of_find_sheet_antiparallel || return_of_find_sheet_parallel ) {
							sheet_unassigned = false;
							WriteToDB_sheet (struct_id, db_session, sheet_PK_id_counter, cur_strands[j].get_sheet_id(), get_segment_id(
								struct_id,
								db_session,
								i)); // struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
							sheet_PK_id_counter++;
						}
					} //for(Size j=1; j<=cur_strands.size(); ++j)

					if ( sheet_unassigned ) {
						Size max_sheet_id = get_max_sheet_id(struct_id, db_session);
						WriteToDB_sheet (struct_id, db_session, sheet_PK_id_counter, max_sheet_id+1, get_segment_id(
							struct_id,
							db_session,
							i)); // struct_id, db_session, sheet_PK_id_counter, sheet_id, segment_id
						sheet_PK_id_counter++;
					}
				}
			} else { //all_strands[i].get_size() >= min_res_in_strand_ // all_strands[i].get_size() < min_res_in_strand_  // "this strand is too small, assign it into '99999' sheet"
				WriteToDB_sheet (struct_id, db_session, sheet_PK_id_counter,
					99999, //sheet_id
					get_segment_id( //// segment_id
					struct_id,
					db_session,
					i));
				sheet_PK_id_counter++;
			}// all_strands[i].get_size() < min_res_in_strand_
		}
		// <end> assignment of strand into rest sheets (other than "1")


		// <begin> redefine sheet id
		bool sheet_id_changed = true; //temp bool
		while ( sheet_id_changed )
				{
			sheet_id_changed = change_sheet_id_if_possible(
				struct_id,
				db_session,
				pose,
				min_CA_CA_dis_,
				max_CA_CA_dis_,
				min_C_O_N_angle_);
		}
		// <end> redefine sheet id


		// <begin> see_whether_sheet_is_antiparallel
		utility::vector1<core::Size> all_distinct_sheet_ids = get_distinct_sheet_id_from_sheet_table(struct_id, db_session);
		for ( Size i=1; i<=all_distinct_sheet_ids.size(); i++ ) {
			if ( all_distinct_sheet_ids[i] == 99999 ) { //all_strands[i].get_size() < min_res_in_strand_
				continue;
			}
			Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

			TR << "i: " << i << endl;
			TR << "num_strands_i: " << num_strands_i << endl;

			// <begin> unit_test specific for 3b83
			if ( (i == 1) && (num_strands_i != 3) ) {
				unit_test_pass_identifier = false;
			}
			if ( (i == 2) && (num_strands_i != 5) ) {
				unit_test_pass_identifier = false;
			}
			// <end> unit_test specific for 3b83

			if ( num_strands_i < min_num_strands_in_sheet_ ) {
				continue;
			}
			string sheet_is_antiparallel = see_whether_sheet_is_antiparallel(
				struct_id,
				db_session,
				pose,
				all_distinct_sheet_ids[i], //sheet id
				min_CA_CA_dis_,
				max_CA_CA_dis_,
				min_C_O_N_angle_
			);
			WriteToDB_sheet_antiparallel(struct_id, db_session, all_distinct_sheet_ids[i], sheet_is_antiparallel);
		}
		// <end> see_whether_sheet_is_antiparallel


		if ( !extract_sandwich_ ) {
			TR.Info << "Exit early since the user doesn't request to extract beta-sandwich " << endl;
			return 0;
		}

		/////////////////// <begin> assignment of sheet into sandwich_candidate_by_sheet (sw_can_by_sh)
		TR.Info << "<begin> assignment of sheet into sandwich_candidate_by_sheet (sw_can_by_sh) " << endl;

		// <begin> assign num_of_sheets_that_surround_this_sheet only ONCE!
		for ( Size i=1; i <= all_distinct_sheet_ids.size(); i++ ) {
			Size num_of_sheets_that_surround_this_sheet =
				cal_num_of_sheets_that_surround_this_sheet(
				struct_id, db_session, dssp_pose, all_distinct_sheet_ids,
				all_distinct_sheet_ids[i],
				min_num_strands_in_sheet_,
				inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_);
			WriteToDB_number_of_sheets_that_surround_this_sheet (struct_id, db_session, all_distinct_sheet_ids[i], num_of_sheets_that_surround_this_sheet);
		}
		// <end> assign num_of_sheets_that_surround_this_sheet only ONCE!

		Size sheet_j_that_will_be_used_for_pairing_with_sheet_i = 0; // temp value
		for ( Size i=1; i <= all_distinct_sheet_ids.size()-1; i++ ) { // now I check all possible combinations
			if ( all_distinct_sheet_ids[i] == sheet_j_that_will_be_used_for_pairing_with_sheet_i ) { // useful to exclude sheet later
				continue;
			}
			if ( all_distinct_sheet_ids[i] == 99999 ) { // all_strands[i].get_size() < min_res_in_strand_
				continue;
			}
			Size num_strands_in_i_sheet = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

			if ( num_strands_in_i_sheet < min_num_strands_in_sheet_ ) {
				continue;
			}
			Size num_of_sheets_that_surround_this_sheet = get_num_of_sheets_that_surround_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);

			if ( num_of_sheets_that_surround_this_sheet > 1 ) {
				continue; // i
			}
			Real minimum_avg_dis_between_sheets = 9999; //temp value
			Real avg_dis_between_sheets;

			// <begin> identify sheet_j that_will_be_used_for_pairing_with_sheet_i to be a sandwich
			// I need to iterate all 'j' for loop to find the closest sheet from sheet_id (i)
			for ( Size j=i+1; j<=all_distinct_sheet_ids.size(); j++ ) {
				//TR << "all_distinct_sheet_ids[j]: " << all_distinct_sheet_ids[j] << endl;
				if ( all_distinct_sheet_ids[j] == sheet_j_that_will_be_used_for_pairing_with_sheet_i ) { // useful to exclude sheet later
					continue;
				}
				if ( all_distinct_sheet_ids[j] == 99999 ) { // all_strands[i].get_size() < min_res_in_strand_
					continue;
				}
				Size num_strands_in_j_sheet = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]); // struct_id, db_session, sheet_id

				if ( num_strands_in_j_sheet < min_num_strands_in_sheet_ ) {
					continue;
				}

				Size num_of_sheets_that_surround_this_sheet = get_num_of_sheets_that_surround_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]);

				if ( num_of_sheets_that_surround_this_sheet > 1 ) {
					continue;  // j
				}


				/////////////////// DO NOT ERASE THESE TRACERS ///////////
				//    TR << "Now a preliminary candidate of sheet pair to be a sw is identified after checking that these sheets are not surrounded by more than 1 other sheet" << endl;
				//    TR.Info << "sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
				//    TR.Info << "sheet_id (all_distinct_sheet_ids[j]): " << all_distinct_sheet_ids[j] << endl;
				/////////////////// DO NOT ERASE THESE TRACERS ///////////


				utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
				utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[j]);


				// <begin> check whether strands are too distant to each other
				bool these_2_sheets_are_too_distant = false; // temporary 'false' designation
				avg_dis_between_sheets = 0;
				for ( Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_distant; ++ii ) { // this '&& !these_2_sheets_are_too_distant' is needed for better performance
					vector<Real> array_avg_dis_sheet_i_j;
					for ( Size jj=1; jj<=strands_from_sheet_j.size(); ++jj ) {
						SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
						SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());

						Real avg_dis_strands = get_avg_dis_strands (pose, temp_strand_i, temp_strand_j);
						array_avg_dis_sheet_i_j.push_back(avg_dis_strands);
					}

					Real sum = 0;
					for ( double k : array_avg_dis_sheet_i_j ) {
						sum += k;
					}
					avg_dis_between_sheets = sum / array_avg_dis_sheet_i_j.size();

					Real shortest_avg_dis_inter_sheet = 9999;

					for ( Size kk=1; kk<=strands_from_sheet_j.size(); ++kk ) {
						if ( shortest_avg_dis_inter_sheet > array_avg_dis_sheet_i_j[kk-1] ) {
							shortest_avg_dis_inter_sheet = array_avg_dis_sheet_i_j[kk-1];
						}
					}
					//TR.Info << "shortest_avg_dis_inter_sheet: " << shortest_avg_dis_inter_sheet << endl;
					if ( shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_ ) {
						// TR.Info << "shortest_avg_dis_inter_sheet > max_inter_sheet_dis_CA_CA_" << endl;
						these_2_sheets_are_too_distant = true;
					}
				} // for(Size ii=1; ii<=strands_from_sheet_i.size(); ++ii)
				// <end> check whether strands are too distant to each other


				//TR.Info << "avg_dis_between_sheets: " << avg_dis_between_sheets << endl;
				if ( these_2_sheets_are_too_distant ) {
					continue; // continue j sheet loop
				}

				if ( minimum_avg_dis_between_sheets > avg_dis_between_sheets ) {
					minimum_avg_dis_between_sheets = avg_dis_between_sheets;

					// <begin> check_strand_too_closeness
					bool these_2_sheets_are_too_close = false; // temporary 'false' designation

					for ( Size ii=1; ii<=strands_from_sheet_i.size() && !these_2_sheets_are_too_close; ++ii ) { // && !these_2_sheets_are_too_close NEEDED for better performance
						for ( Size jj=1; jj<=strands_from_sheet_j.size() && !these_2_sheets_are_too_close; ++jj ) { // && !these_2_sheets_are_too_close NEEDED for better performance
							SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
							SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());

							bool are_strands_too_close =
								check_strand_too_closeness (
								pose,
								temp_strand_i,
								temp_strand_j,
								min_inter_sheet_dis_CA_CA_);
							if ( are_strands_too_close ) {
								//  TR.Info << "these two sheets are too close when I calculate its distance by their strands" << endl;
								these_2_sheets_are_too_close = true;
							}
						}
					}
					if ( these_2_sheets_are_too_close ) {
						continue; // continue in j sheet loop
					}
					// <end> check_strand_too_closeness

					sheet_j_that_will_be_used_for_pairing_with_sheet_i = all_distinct_sheet_ids[j]; // use all_distinct_sheet_ids[j] to pair with all_distinct_sheet_ids[i]
				}


			} // for(Size j=i+1; j<=all_distinct_sheet_ids.size(); ++j)

			// <end> identify sheet_j_that_will_be_used_for_pairing_with_sheet_i to be a sandwich


			if ( sheet_j_that_will_be_used_for_pairing_with_sheet_i == 0 ) {
				continue; // continue i sheet 'for' loop
			}

			/////////////////// DO NOT ERASE THESE TRACERS ///////////
			//   TR.Info << "Now a real pair of candidates between the closest sheets is identified " << endl;
			//   TR.Info << "candidate 1: sheet_id (all_distinct_sheet_ids[i]): " << all_distinct_sheet_ids[i] << endl;
			//   TR.Info << "candidate 2: sheet_id (sheet_j_that_will_be_used_for_pairing_with_sheet_i): " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << endl;
			/////////////////// DO NOT ERASE THESE TRACERS ///////////

			// <begin> check_sw_by_distance
			bool found_sandwich_w_these_2_sheets = false; // temporary 'false' designation
			bool chance_of_being_sandwich_w_these_2_sheets = true; // temporary 'true' designation
			utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, all_distinct_sheet_ids[i]);
			utility::vector1<SandwichFragment> strands_from_sheet_j = get_full_strands_from_sheet(struct_id, db_session, sheet_j_that_will_be_used_for_pairing_with_sheet_i);

			//   TR.Info << "<begin> check_sw_by_distance" << endl;
			while ( !found_sandwich_w_these_2_sheets && chance_of_being_sandwich_w_these_2_sheets )
					{
				for ( Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich_w_these_2_sheets; ++ii ) {
					for ( Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich_w_these_2_sheets; ++jj ) {
						SandwichFragment temp_strand_i(strands_from_sheet_i[ii].get_start(), strands_from_sheet_i[ii].get_end());
						SandwichFragment temp_strand_j(strands_from_sheet_j[jj].get_start(), strands_from_sheet_j[jj].get_end());

						Real return_of_check_sw_by_dis_anti = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, true, min_sheet_dis_, max_sheet_dis_);
						Real return_of_check_sw_by_dis_parallel = check_sw_by_dis (pose, temp_strand_i, temp_strand_j, false, min_sheet_dis_, max_sheet_dis_);

						if ( return_of_check_sw_by_dis_anti == -999 || return_of_check_sw_by_dis_parallel == -999 ) {
							TR.Info << "these sheets will not be sandwich ever because these are too close or distant to each other!" << endl;
							chance_of_being_sandwich_w_these_2_sheets = false;
						}

						if ( return_of_check_sw_by_dis_anti != -99 || return_of_check_sw_by_dis_parallel != -99 ) {
							TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " are in the ideal distance range" << endl;
							found_sandwich_w_these_2_sheets = true;
							chance_of_being_sandwich_w_these_2_sheets = false; // these are sandwich, but no more sheet search is needed! (this "false" false assignment is needed (confirmed! by experiment))
						}
					} //for(Size jj=1; jj<=strands_from_sheet_j.size() && chance_of_being_sandwich_w_these_2_sheets; ++jj)
				} //for(Size ii=1; ii<=strands_from_sheet_i.size() && chance_of_being_sandwich_w_these_2_sheets; ++ii)
				break; // no sandwich here
			} //while (!found_sandwich_w_these_2_sheets && chance_of_being_sandwich_w_these_2_sheets)
			//   TR.Info << "<end> check_sw_by_distance" << endl;
			// <end> check_sw_by_distance

			if ( !found_sandwich_w_these_2_sheets ) {
				continue;
			}


			if ( exclude_sandwich_that_is_suspected_to_have_not_facing_2_sheets_ ) {
				int facing =
					judge_facing(
					struct_id,
					db_session,
					pose,
					all_distinct_sheet_ids[i], // sheet_i
					sheet_j_that_will_be_used_for_pairing_with_sheet_i, // sheet_j
					min_CA_CA_dis_,
					max_CA_CA_dis_,
					min_sheet_angle_by_four_term_cen_res_,
					max_sheet_angle_by_four_term_cen_res_,
					min_sheet_torsion_cen_res_,
					max_sheet_torsion_cen_res_,
					max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_);
				// if false, these two strand_pairs are linear to each other or do not face each other properly

				if ( facing == 0 ) {
					TR.Info << "sheet " << all_distinct_sheet_ids[i] << " and sheet " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << " seem not facing to each other" << endl;
					continue; // skip this sheet
				} else if ( facing == -99 ) {
					TR.Info << "at least one sheet (either " << all_distinct_sheet_ids[i] << " or " << sheet_j_that_will_be_used_for_pairing_with_sheet_i << ")  may be a beta-barrel like sheet_id = 1 in 1N8O" << endl;
					continue; // skip this sheet
				}
			}

			TR.Info << "writing into 'sandwich candidate by sheet'" << endl;
			WriteToDB_sw_can_by_sh (struct_id, db_session, sw_can_by_sh_PK_id_counter, tag, sw_can_by_sh_id_counter, all_distinct_sheet_ids[i], strands_from_sheet_i.size());
			sw_can_by_sh_PK_id_counter++;

			WriteToDB_sw_can_by_sh (struct_id, db_session, sw_can_by_sh_PK_id_counter, tag, sw_can_by_sh_id_counter, sheet_j_that_will_be_used_for_pairing_with_sheet_i, strands_from_sheet_j.size());
			sw_can_by_sh_PK_id_counter++;

			sw_can_by_sh_id_counter++;

		} // for(Size i=1; i<all_distinct_sheet_ids.size(); ++i)
		/////////////////// <end> assignment of sheet into sw_can_by_sh


		/////////////////// <begin> fill a table 'sandwich' by secondary_structure_segments
		TR.Info << "<begin> fill a table 'sandwich' by secondary_structure_segments" << endl;

		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh = prepare_WriteToDB_sandwich(struct_id, db_session);
		// It retrieves all beta-strands of sandwich_candidate_by_sheets, it does not make sandwich_by_components

		if ( bs_of_sw_can_by_sh.size() == 0 ) {
			TR.Info << endl << "no beta segment in sandwich_by_sheet " << endl;
			TR.Info << "(maybe these two sheets do not face each other <OR> " << endl;
			TR.Info << "there are only < " << min_num_strands_in_sheet_ << " number of strands in one sheet <OR> " << endl;
			TR.Info << "these are too distant sheets <OR> " << endl;
			TR.Info << "this is a beta barrel <OR> " << endl;
			TR.Info << "\"non-canonical\" like 1MSP"") " << endl << endl;
			TR.Info << "<Exit-Done> for this pdb including extraction of sandwich" << endl;
			return 0;
		}

		sandwich_PK_id_counter =
			Run_WriteToDB_sandwich(
			tag,
			dssp_pose,
			bs_of_sw_can_by_sh,
			max_num_sw_per_pdb_,
			struct_id, // needed argument
			db_session,
			min_CA_CA_dis_,
			max_CA_CA_dis_,
			sandwich_PK_id_counter);

		write_phi_psi_of_each_residue_to_a_file(
			tag,
			dssp_pose,
			bs_of_sw_can_by_sh,
			write_phi_psi_of_E_,
			write_phi_psi_of_all_,
			max_num_sw_per_pdb_,
			struct_id, // needed argument
			db_session,
			min_CA_CA_dis_,
			max_CA_CA_dis_);

		if ( count_AA_with_direction_ ) {
			//// <begin> count AA with direction
			for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii ) { // per each beta-strand
				if ( bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_ ) {
					break;
				}
				WriteToDB_sandwich_by_AA_w_direction (struct_id, db_session, pose, pose_w_center_000, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			}
			//// <end> count AA with direction
		} //count_AA_with_direction_

		/////////////////// <end> fill a table 'sandwich' by secondary_structure_segments
		TR.Info << "<end> fill a table 'sandwich' by secondary_structure_segments" << endl;


		if ( do_not_connect_sheets_by_loops_ ) {
			TR.Info << "<Exit> do_not_connect_sheets_by_loops_:" << do_not_connect_sheets_by_loops_ << endl;
			return 0;
		}

		/////////////////// <begin> update beta-hairpin or inter_sheet_connecting_loops (2nd judgement whether each sandwich_by_sheet_id can become a sandwich_by_components)

		// get_distinct(sw_can_by_sh_id)
		utility::vector1<core::Size> vec_sw_can_by_sh_id =  get_vec_of_sw_can_by_sh_id(struct_id, db_session);

		for ( Size ii=1; ii<=vec_sw_can_by_sh_id.size(); ii++ ) {
			// I think that mostly vec_sw_can_by_sh_id.size() = just 1
			TR << "See whether 'sw_candidate_by_sheets_id = " << vec_sw_can_by_sh_id[ii] << " ' be a canonical beta-sandwich or not" << endl;

			chance_of_being_canonical_sw = true; // not yet decided fate whether this could be canonical sandwich or not, but assumed to be true for now
			Size size_sandwich_PK_id =
				get_size_sandwich_PK_id(
				struct_id,
				db_session,
				vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
			);

			bool bool_proper_num_helix_in_loop = true;
			bool bool_proper_num_E_in_loop = true;
			Size former_start_res = 0; //temporary

			// this 'jj' is used for 'for' iteration purpose only and this 'for' loop iterates only for connecting sheets/strands
			for ( Size jj=1; jj<=size_sandwich_PK_id-1; ++jj ) {
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

				if ( start_res == 0 ) { // no proper retrieval of start_res
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
				for ( Size kk=start_res+1; kk<=next_start_res-1; ++kk ) {
					char res_ss( dssp_pose.secstruct( kk ) ) ;
					if ( res_ss == 'H' ) {
						helix_num += 1;
					}
				}
				if ( helix_num > max_H_in_extracted_sw_loop_ ) {
					if ( TR.Debug.visible() ) {
						TR.Debug << "helix_num: " << helix_num << ", max_H_in_extracted_sw_loop_: " << max_H_in_extracted_sw_loop_ << std::endl;
					}
					bool_proper_num_helix_in_loop = false;
				}
				//////////// <end> check whether there is a helix as a loop in this extracted sandwich candidate


				//////////// <begin> check whether there is a strand as a loop in this extracted sandwich candidate
				Size E_num = 0;
				for ( Size kk=start_res+1; kk<=next_start_res-1; ++kk ) {
					char res_ss( dssp_pose.secstruct( kk ) ) ;
					if ( res_ss == 'E' ) {
						E_num += 1;
					}
				}
				if ( E_num > max_E_in_extracted_sw_loop_ ) {
					TR << "E_num > max_E_in_extracted_sw_loop_ " << endl;
					bool_proper_num_E_in_loop = false;
				}
				//////////// <end> check whether there is a strand as a loop in this extracted sandwich candidate


				if ( !bool_proper_num_helix_in_loop || !bool_proper_num_E_in_loop ) {
					delete_this_sw_can_by_sh_id_from_sw_by_comp(
						struct_id,
						db_session,
						vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);
					chance_of_being_canonical_sw = false;
					break; // break jj 'for' loop
				}
				//////////// <begin> check numbers of helix and strand residues in loops


				//string LR = check_LR(pose, start_res+1, next_start_res-1);
				string LR = check_LR(dssp_pose, start_res+1, next_start_res-1);

				std::pair<string, string> PA = check_PA(dssp_pose, start_res+1, next_start_res-1);
				string PA_by_preceding_E = PA.first;
				string PA_by_following_E = PA.second;

				// check whethere it is positive, negative, away or meet
				string heading_direction = check_heading_direction(dssp_pose, start_res+1, next_start_res-1, check_N_to_C_direction_by_);

				if ( heading_direction == "except" ) {
					TR.Info << "Exit-Exception:: check_N_to_C_direction_by should be either PF or FE !" << endl;
					return 0;
				}

				string parallel_EE;
				if ( heading_direction == "posi" || heading_direction == "nega" ) {
					parallel_EE = "P_EE";
				} else {
					parallel_EE = "A_EE"; // to keep chacracter size be same as "parallel"
				}

				Size loop_size = (next_start_res-1) - (start_res+1) + 1;
				if ( sheet_id_of_start_res == sheet_id_of_next_start_res ) {
					// this loop is a beta-hairpin loop (that connects sheets as intra-sheet way)
					//     TR << "start_res: " << start_res << endl;
					//     TR << "next_start_res: " << next_start_res << endl;
					bool loop_is_surrounded_by_same_direction_strands = check_whether_same_direction_strands_connect_two_sheets_or_a_loop(struct_id, db_session, pose, start_res, next_start_res);
					//TR.Info << "loop_is_surrounded_by_same_direction_strands: " << loop_is_surrounded_by_same_direction_strands << endl;

					string canonical_LR = "-";
					string cano_PA =  "-";
					string cano_parallel_EE =  "-";
					string loop_kind =  "loop_connecting_same_direction_strands_within_same_sheet";

					bool hairpin_connects_short_strand = check_whether_hairpin_connects_short_strand(struct_id, db_session, start_res, next_start_res);

					if ( !loop_is_surrounded_by_same_direction_strands ) {
						canonical_LR = check_canonicalness_of_LR(loop_size, true, LR); // loop_size, intra_sheet bool, LR

						if ( do_not_write_resfile_of_sandwich_that_has_non_canonical_LR_ && canonical_LR == "F_LR" ) {
							TR.Info << "This sandwich has non-canonical hairpin chirality " << endl;
							TR.Info << "So don't write resfile for this sandwich for subsequent design" << endl;
							write_resfile_ = false;
							TR.Info << "this pdb file fails to pass ChiralitySandwichFilter, this procedure is obviously risky when extracting multiple sandwiches from raw pdb files, so use this ChiralitySandwichFilter only when each pdb file has 1 sandwich" << endl;
							//return 1; // which means this sandwich fail to pass ChiralitySandwichFilter
							//       TR.Info << "this sandwich has non-canonical hairpin chirality, so delete " << vec_sw_can_by_sh_id[ii] << " sw_can_by_sh_id" << endl;
							//      delete_this_sw_can_by_sh_id_from_sw_by_comp(
							//       struct_id,
							//       db_session,
							//       vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
							//       );
							//      chance_of_being_canonical_sw = false;
						}

						if ( exclude_sandwich_that_has_non_canonical_LR_ && canonical_LR == "F_LR" ) {
							TR.Info << "This sandwich has non-canonical hairpin chirality " << endl;
							delete_this_struct_id(
								struct_id,
								db_session
							);
							return 1;
							//      chance_of_being_canonical_sw = false;
						}

						cano_PA = check_canonicalness_of_PA(loop_size, true, PA_by_preceding_E, PA_by_following_E, check_canonicalness_cutoff_); // loop_size, intra_sheet bool, 2 PAs
						cano_parallel_EE = check_canonicalness_of_parallel_EE(loop_size, true, parallel_EE); // loop_size, intra_sheet bool, parallel_EE

						if ( hairpin_connects_short_strand ) { // so it is not recommended for LR/PA analysis
							loop_kind =  "hairpin_connecting_short_strand";
						} else { // this hairpin-loop is ideal for LR/PA
							loop_kind =  "hairpin";
						}

					}

					WriteToDB_sheet_connectivity(
						struct_id,
						db_session,
						pose,
						sandwich_PK_id_counter,
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

					if ( !loop_is_surrounded_by_same_direction_strands && !hairpin_connects_short_strand ) { // so it is recommended for LR/PA analysis
						if ( ((next_start_res-1) - (start_res+1)) == 1 ) { // this loop is a 2-residues long hairpin
							string turn_type = WriteToDB_turn_type
								(pose,
								vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
								start_res+1, // new start_res for intra_sheet_con
								next_start_res-1, // new end_res for intra_sheet_con
								struct_id,
								db_session,
								allowed_deviation_for_turn_type_id_
							);

							WriteToDB_turn_AA
								(pose,
								vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
								start_res, // i
								struct_id,
								db_session,
								turn_type
							);
						}
					}

					sandwich_PK_id_counter++;
					intra_sheet_con_id_counter++;
				} else { // this loop connects sheets as inter-sheet way
					if ( exclude_sandwich_that_is_linked_w_same_direction_strand_ ) {
						bool sheets_are_connected_by_same_direction_strand = check_whether_same_direction_strands_connect_two_sheets_or_a_loop(struct_id, db_session, pose, start_res, next_start_res);
						if ( sheets_are_connected_by_same_direction_strand ) {
							TR.Info << "sheets_are_connected_by_same_direction_strand, so delete " << vec_sw_can_by_sh_id[ii] << " sw_can_by_sh_id" << endl;
							delete_this_sw_can_by_sh_id_from_sw_by_comp(
								struct_id,
								db_session,
								vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
							);
							chance_of_being_canonical_sw = false;
							break; // break jj 'for' loop
						}
					}


					string canonical_LR = check_canonicalness_of_LR(loop_size, false, LR); // loop_size, intra_sheet bool, LR
					string cano_PA = check_canonicalness_of_PA(loop_size, false, PA_by_preceding_E, PA_by_following_E, check_canonicalness_cutoff_);
					// loop_size, intra_sheet bool, PA_ref_1, PA_ref_2, cutoff
					string cano_parallel_EE = check_canonicalness_of_parallel_EE(loop_size, false, parallel_EE);
					// loop_size, intra_sheet bool, parallel_EE
					WriteToDB_sheet_connectivity(
						struct_id,
						db_session,
						pose,
						sandwich_PK_id_counter,
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

					if ( ((next_start_res-1) - (start_res+1)) == 1 ) { // this loop is a 2-residues long inter-sheet_loop
						string turn_type = WriteToDB_turn_type
							(pose,
							vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
							start_res+1, // new start_res for intra_sheet_con
							next_start_res-1, // new end_res for intra_sheet_con
							struct_id,
							db_session,
							allowed_deviation_for_turn_type_id_
						);

						WriteToDB_turn_AA
							(pose,
							vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
							start_res, // i
							struct_id,
							db_session,
							turn_type
						);
					}

					sandwich_PK_id_counter++;
					inter_sheet_con_id_counter++;
				} // this loop connects sheets as inter-sheet way

			} // for(Size jj=1; (jj<=size_sandwich_PK_id-1) && (bool_no_helix_in_loop) && (bool_no_more_strand_in_loop); ++jj)
			/////////////////// <end> update beta-hairpin or inter_sheet_connecting_loops (2nd judgement whether each sandwich_by_sheet_id becomes sandwich_by_components)


			if ( exclude_sandwich_that_has_near_backbone_atoms_between_sheets_ ) {
				bool sheets_are_connected_with_near_bb_atoms =
					check_whether_sheets_are_connected_with_near_bb_atoms(
					struct_id,
					db_session,
					dssp_pose,
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id,
					min_N_O_dis_between_two_sheets_,
					min_N_H_O_angle_between_two_sheets_
				);
				if ( sheets_are_connected_with_near_bb_atoms ) {
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
			if ( (!unit_test_pass_identifier) || (!chance_of_being_canonical_sw) ) {
				unit_test_pass_identifier = false;
			}

			if ( chance_of_being_canonical_sw ) {
				canonical_sw_extracted_from_this_pdb_file = true;

				WriteToDB_starting_loop(
					struct_id,
					db_session,
					dssp_pose,
					sandwich_PK_id_counter, //sandwich_PK_id
					vec_sw_can_by_sh_id[ii], //sw_can_by_sh_id
					tag,
					max_starting_loop_size_);
				sandwich_PK_id_counter++;

				WriteToDB_ending_loop(
					struct_id,
					db_session,
					dssp_pose,
					sandwich_PK_id_counter, //sandwich_PK_id
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
					tag,
					max_starting_loop_size_);

				sandwich_PK_id_counter++;


				// <begin> mark beta-sandwiches that is not connected by continuous atoms like 1A78
				string sw_is_not_connected_with_continuous_atoms = check_whether_sw_is_not_connected_with_continuous_atoms(
					struct_id,
					db_session,
					dssp_pose,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);
				//TR.Info << "sw_is_not_connected_with_continuous_atoms: " << sw_is_not_connected_with_continuous_atoms << endl;

				WriteToDB_whether_sw_is_not_connected_with_continuous_atoms(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
					sw_is_not_connected_with_continuous_atoms);
				// <end> mark beta-sandwiches that is not connected by continuous atoms like 1A78


				WriteToDB_dssp_ratio_in_sw(struct_id, db_session, dssp_pose,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				if ( count_AA_with_direction_ ) {
					WriteToDB_number_of_core_heading_FWY_in_sw(struct_id, db_session,
						vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);

					WriteToDB_number_of_core_heading_W_in_sw(struct_id, db_session,
						vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);


					WriteToDB_number_of_core_heading_LWY_in_core_strands_in_sw(struct_id, db_session,
						vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
					);

					WriteToDB_ratio_of_core_heading_FWY_in_sw(struct_id, db_session,
						vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
						pose
					);

				}

				Size sw_res_size = WriteToDB_sw_res_size(struct_id, db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				WriteToDB_avg_b_factor_CB_at_each_component(struct_id, db_session, pose,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				WriteToDB_dihedral_angle_between_core_strands_across_facing_sheets(struct_id, db_session, pose,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				WriteToDB_hydrophobic_ratio_net_charge(struct_id, db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				WriteToDB_long_strand_id_in_each_sw(struct_id, db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				WriteToDB_min_avg_dis_between_sheets_by_cen_res (struct_id, db_session,
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
					dssp_pose,
					all_distinct_sheet_ids,
					min_num_strands_in_sheet_
				);

				WriteToDB_min_dis_between_sheets_by_all_res (struct_id, db_session,
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
					dssp_pose,
					all_distinct_sheet_ids
				);

				WriteToDB_number_of_edge_strands_in_each_sw(struct_id, db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				WriteToDB_number_strands_in_each_sw(struct_id, db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);

				/*

				// 'WriteToDB_prolines_that_seem_to_prevent_aggregation' is defined in WriteToDBFromSandwichFeatures.hh

				WriteToDB_prolines_that_seem_to_prevent_aggregation(
				struct_id,
				db_session,
				vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
				wt_for_pro_in_starting_loop_,
				wt_for_pro_in_1st_inter_sheet_loop_,
				wt_for_pro_in_3rd_inter_sheet_loop_
				);
				*/

				WriteToDB_topology_candidate(struct_id, db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);


				Real shortest_dis_between_facing_aro_in_sw =
					WriteToDB_shortest_dis_between_facing_aro_in_sw (
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
					pose,
					all_distinct_sheet_ids,
					min_num_strands_in_sheet_
				);

				if ( exclude_sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw_ && (shortest_dis_between_facing_aro_in_sw < 4.0) ) {
					TR.Info << "This sandwich has an non_canonical_shortest_dis_between_facing_aro_in_sw " << endl;
					delete_this_struct_id(
						struct_id,
						db_session
					);
					return 1;
				}

				if ( exclude_sandwich_with_SS_bond_ ) {
					bool this_sw_has_SS_bond = see_whether_this_sw_has_SS_bond(
						struct_id,
						db_session);
					if ( this_sw_has_SS_bond ) {
						TR.Info << "This sandwich has S-S bond(s) " << endl;
						// those sandwiches having S-S bond in core tend to have big holes (bad packing score), so exclude them
						delete_this_struct_id(
							struct_id,
							db_session
						);
						return 1;
					}
				}

				if ( write_AA_kind_files_ ) {
					write_AA_kind_to_a_file(
						tag,
						struct_id,
						db_session,
						vec_sw_can_by_sh_id[ii], // sw_can_by_sh_id
						sw_res_size
					);
				}

			} else { // chance_of_being_canonical_sw = false
				delete_this_sw_can_by_sh_id_from_sw_by_comp(
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);
			}

			// <begin> write_chain_B_resNum_to_a_file
			if ( write_chain_B_resnum_ && chance_of_being_canonical_sw ) {
				write_chain_B_resNum_to_a_file(
					tag,
					struct_id,
					db_session,
					vec_sw_can_by_sh_id[ii] // sw_can_by_sh_id
				);
			}
			// <end> write_chain_B_resNum_to_a_file

		} // per each sandwich_candidate_by_sheet_id

		// write_heading_direction_of_all_AA_in_a_strand_to_a_file
		if ( write_heading_directions_of_all_AA_in_a_strand_ && canonical_sw_extracted_from_this_pdb_file ) {
			write_heading_direction_of_all_AA_in_a_strand_to_a_file(
				tag,
				struct_id, // needed argument
				db_session,
				pose,
				bs_of_sw_can_by_sh);
		} //write_heading_direction_of_all_AA_in_a_strand_to_a_file

		// prepare_and_write_number_of_electrostatic_interactions_of_residues_to_files
		if ( canonical_sw_extracted_from_this_pdb_file ) {
			prepare_and_write_number_of_electrostatic_interactions_of_residues_to_files(
				tag,
				struct_id,
				db_session,
				pose,
				bs_of_sw_can_by_sh,
				distance_cutoff_for_electrostatic_interactions_,
				CB_b_factor_cutoff_for_electrostatic_interactions_,
				min_primary_seq_distance_diff_for_electrostatic_interactions_,
				write_electrostatic_interactions_of_surface_residues_in_a_strand_,
				write_electrostatic_interactions_of_all_residues_in_a_strand_,
				write_electrostatic_interactions_of_all_residues_,
				rkde_PK_id_counter,
				rkde_in_strands_PK_id_counter);
		} //prepare_and_write_number_of_electrostatic_interactions_of_residues_to_files


		// <begin> write_beta_sheet_capping_info_
		if ( write_beta_sheet_capping_info_ && canonical_sw_extracted_from_this_pdb_file ) {
			write_beta_sheet_capping_info_to_a_file(
				tag,
				pose,
				bs_of_sw_can_by_sh,
				primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping_,
				primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping_,
				primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping_,
				primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping_);
		}
		// <end> write_beta_sheet_capping_info_


		// <begin> write_AA_distribution_with_direction_to_a_file
		if ( write_strand_AA_distribution_files_ && canonical_sw_extracted_from_this_pdb_file ) {
			if ( !count_AA_with_direction_ ) {
				TR << "You did not turn on count_AA_with_direction for write_strand_AA_distribution_files so SandwicheFeature will not write strand_AA_distribution_files" << endl;
			} else {
				write_AA_distribution_with_direction_to_a_file(
					tag,
					struct_id,
					db_session);
			}
		}
		// <end> write_AA_distribution_with_direction_to_a_file


		// not needed now (was needed for just during development
		///*
		{
			core::scoring::ScoreFunctionOP centroid_scorefxn( generate_scorefxn( false //fullatom
				) );
			core::scoring::ScoreFunctionOP fullatom_scorefxn( generate_scorefxn( true //fullatom
				) );

			process_decoy( dssp_pose, pose.is_fullatom() ? *fullatom_scorefxn : *centroid_scorefxn);

			//TR.Info << "Current energy terms that we deal with:" << endl;
			//dssp_pose.energies().show_total_headers( std::cout );
			//std::cout << std::endl;
		}
		//*/

		// <begin> write p_aa_pp to a file
		if ( write_p_aa_pp_files_ && canonical_sw_extracted_from_this_pdb_file ) {
			write_p_aa_pp_of_AAs_to_a_file(
				tag,
				dssp_pose);
		}
		// <end> write p_aa_pp to a file


		// <begin> write rama of residue to a file
		if ( write_rama_at_AA_to_files_ && canonical_sw_extracted_from_this_pdb_file ) {
			write_rama_of_AAs_to_a_file(
				tag,
				dssp_pose);
		}
		// <end> write rama of residue to a file

		// <begin> write AA_distribution_without_direction to a file
		if ( write_loop_AA_distribution_files_ && canonical_sw_extracted_from_this_pdb_file ) {
			write_AA_distribution_without_direction_to_a_file(
				tag,
				struct_id,
				db_session);
		}
		// <end> write AA_distribution_without_direction to a file


		//// <begin> WriteToDB number_of_core_heading_charged_AAs/aro_AAs_in_a_pair_of_edge_strands
		if ( count_AA_with_direction_ && canonical_sw_extracted_from_this_pdb_file ) {
			WriteToDB_number_of_AAs_in_a_pair_of_edge_strands (
				struct_id,
				db_session,
				pose,
				bs_of_sw_can_by_sh,
				max_num_sw_per_pdb_,
				min_CA_CA_dis_,
				max_CA_CA_dis_);
		}
		//// <end> WriteToDB number_of_core_heading_charged_AAs/aro_AAs_in_a_pair_of_edge_strands


		// <begin> write resfile automatically
		if ( write_resfile_ && canonical_sw_extracted_from_this_pdb_file ) {
			// (07/10/2014) This automatic resfile generation seems useful only for design with ramping repulsions, for PackRotamersMover & OffRotamerPack, this kind of resfile seems not needed.
			if ( write_resfile_when_seq_rec_is_bad_ ) {
				write_resfile_to_a_file_when_seq_rec_is_bad(
					tag,
					struct_id, // needed argument
					db_session, // needed argument
					pose,
					bs_of_sw_can_by_sh,
					write_resfile_to_minimize_too_much_hydrophobic_surface_,
					write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_,
					write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_);
			}
			write_resfile_to_a_file(
				tag,
				struct_id,      // needed argument
				db_session,     // needed argument
				pose,
				bs_of_sw_can_by_sh,
				write_resfile_NOT_FWY_on_surface_);
		}
		// <end> write resfile automatically


		// development
		// <begin> assign SS to blueprint file
		ObjexxFCL::FArray1D_char dsspSS( pose.size() );
		dssp.dssp_reduced(dsspSS);

		//TR << "input PDB dssp assignment: (based on start structure)" << std::endl;
		for ( Size i = 1; i <= pose.size(); i++ ) {
			///TR << dsspSS(i);
		}
		//TR << std::endl;
		// <end> assign SS to blueprint file
		// development


		TR.Info << "<Exit-Done> for this pdb including extraction of sandwich" << endl;

	} // try

catch ( utility::excn::Exception const & e )
{
	std::cout << "caught exception during SandwichFeatures" << e.msg() << std::endl;
} // catch


	// for unit_test, somehow hardcoded right now
	if ( unit_test_pass_identifier ) {
		return 1; // passed well
	} else {
		return 0; // something's wrong!
	}

} //SandwichFeatures::report_features

std::string SandwichFeatures::type_name() const {
	return class_name();
}

std::string SandwichFeatures::class_name() {
	return "SandwichFeatures";
}

void SandwichFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("min_num_strands_to_deal", xsct_non_negative_integer, "Minimum number of strands in beta sandwich", "4")
		+ XMLSchemaAttribute::attribute_w_default("max_num_strands_to_deal", xsct_non_negative_integer, "Maximum number of strands to consider", "140")
		+ XMLSchemaAttribute::attribute_w_default("min_res_in_strand", xsct_non_negative_integer, "Minimum number of residues per beta strand to count it", "2")
		+ XMLSchemaAttribute::attribute_w_default("min_CA_CA_dis", xsct_real, "Minimum Calpha distance between strands (in Angstroms)", "3.5")
		+ XMLSchemaAttribute::attribute_w_default("max_CA_CA_dis", xsct_real, "Maximum Calpha distance between strands (in Angstroms)", "6.2")
		+ XMLSchemaAttribute::attribute_w_default("min_N_O_dis_between_sheets", xsct_real, "Minimum N-O distance between beta sheets", "3.3")
		+ XMLSchemaAttribute::attribute_w_default("min_N_H_O_angle_between_two_sheets", xsct_real, "Minimum N-H-O distance between two beta sheets", "154.0")
		+ XMLSchemaAttribute::attribute_w_default("min_C_O_N_angle", xsct_real, "Minimum C-O-N angle for a hydrogen bond between strands", "120.0")
		+ XMLSchemaAttribute::attribute_w_default("min_sheet_dis", xsct_real, "Minimum distance between beta sheets in a sandwich (in Angstroms)", "7.0")
		+ XMLSchemaAttribute::attribute_w_default("max_sheet_dis", xsct_real, "Maximum distance between beta sheets in a sandwich (in Angstroms)", "15.0")
		+ XMLSchemaAttribute::attribute_w_default("max_sheet_angle_with_cen_res", xsct_real, "Maximum angle of the plane of the beta sheet with its center residue (defines the twist of the sheet)", "130.0")
		+ XMLSchemaAttribute::attribute_w_default("min_sheet_angle_by_four_term_cen_res", xsct_real, "Minimum sheet angle defined by the four terminal residues and central residue", "25.0")
		+ XMLSchemaAttribute::attribute_w_default("max_sheet_angle_by_four_term_cen_res", xsct_real, "Maximum sheet angle defined by the four terminal residues and central residue", "150.0")
		+ XMLSchemaAttribute::attribute_w_default("min_sheet_torsion_cen_res", xsct_real, "Minimum sheet torsion about the center residue", "-150.0")
		+ XMLSchemaAttribute::attribute_w_default("max_sheet_torsion_cen_res", xsct_real, "Maximum sheet torsion about the central residue", "150.0")
		+ XMLSchemaAttribute::attribute_w_default("min_num_strands_in_sheet", xsct_non_negative_integer, "Minimum number of strands per beta sheet", "2")
		+ XMLSchemaAttribute::attribute_w_default("min_inter_sheet_dis_CA_CA", xsct_real, "Minimum Calpha distance between sheets", "4.0")
		+ XMLSchemaAttribute::attribute_w_default("max_inter_sheet_dis_CA_CA", xsct_real, "Maximum Calpha distance between sheets", "24.0")
		+ XMLSchemaAttribute::attribute_w_default("max_inter_strand_angle_to_not_be_same_direction_strands", xsct_real, "Maximum angle between strands at which the strands will still be considered parallel", "120.0")
		+ XMLSchemaAttribute::attribute_w_default("max_abs_inter_strand_dihedral_to_not_be_same_direction_strands", xsct_real, "Maximum absolute value of the interstrand dihedral angle for them to be considered antiparallel", "100.0")
		+ XMLSchemaAttribute::attribute_w_default("max_num_sw_per_pdb", xsct_non_negative_integer, "Maximum number of sandwiches per PDB", "100")
		+ XMLSchemaAttribute::attribute_w_default("check_N_to_C_direction_by", xs_string, "How to determine the N to C direction?", "PE")
		+ XMLSchemaAttribute::attribute_w_default("check_canonicalness_cutoff", xsct_real, "Cutoff for considering a sandwich canonical", "80.0")
		+ XMLSchemaAttribute::attribute_w_default("count_AA_with_direction", xsct_rosetta_bool, "Count amino acids in each direction?", "false")
		+ XMLSchemaAttribute::attribute_w_default("do_not_write_resfile_of_sandwich_that_has_non_canonical_LR", xsct_rosetta_bool, "Do not write resfile of beta sandwiches with non-canonical loop regions", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_desinated_pdbs", xsct_rosetta_bool, "Exclude designated PDBs", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_that_has_near_backbone_atoms_between_sheets", xsct_rosetta_bool, "Exclude beta sandwiches with nearby backbone atoms between the two sheets", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_that_has_non_canonical_LR", xsct_rosetta_bool, "Exclude sandwiches with noncanonical loop regions", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_that_has_non_canonical_properties", xsct_rosetta_bool, "Exclude noncanonical beta sandwiches", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw", xsct_rosetta_bool, "Exclude sandwiches with noncanonical shortest distances between aromatics in the sandwich", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_that_is_linked_w_same_direction_strand", xsct_rosetta_bool, "Exclude sanwiches linked with a parallel strand?", "false")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_that_is_suspected_to_have_not_facing_2_sheets", xsct_rosetta_bool, "Exclude sandwiches that do not have two facing sheets?", "true")
		+ XMLSchemaAttribute::attribute_w_default("exclude_sandwich_with_SS_bond", xsct_rosetta_bool, "Exclude disulfide bonded sandwiches?", "false")
		+ XMLSchemaAttribute::attribute_w_default("max_starting_loop_size", xsct_non_negative_integer, "Maximum starting loop size", "6")
		+ XMLSchemaAttribute::attribute_w_default("max_ending_loop_size", xsct_non_negative_integer, "Maximum ending loop size", "6")
		+ XMLSchemaAttribute::attribute_w_default("max_E_in_extracted_sw_loop", xsct_non_negative_integer, "Maximum strand size in extracted sandwich loop", "10")
		+ XMLSchemaAttribute::attribute_w_default("max_H_in_extracted_sw_loop", xsct_non_negative_integer, "Maximum strand size in the extracted sandwich loop", "100")
		+ XMLSchemaAttribute::attribute_w_default("no_helix_in_pdb", xsct_rosetta_bool, "Do not include sandwiches with helices in the PDB", "false")
		+ XMLSchemaAttribute::attribute_w_default("inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets", xsct_real, "Minimum cutoff for inter-sheet distance for a sheet to be associated with surrounding sheets. WIthin this distance, they are too close to each other.", "13.0")
		+ XMLSchemaAttribute::attribute_w_default("allowed_deviation_for_turn_type_id", xsct_real, "Maximum deviation for turn type id", "40.0")
		+ XMLSchemaAttribute::attribute_w_default("primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping", xs_integer, "maximum distance from beta sheet capping and N terminus", "2")
		+ XMLSchemaAttribute::attribute_w_default("primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping", xs_integer, "How many residues after the N terminal capping residues must beta sheet capping residues occur", "0")
		+ XMLSchemaAttribute::attribute_w_default("primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping", xs_integer, "How many residues before the C terminal capping residues must beta sheet capping residues occur", "0")
		+ XMLSchemaAttribute::attribute_w_default("primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping", xs_integer, "How many residues after the C terminal capping residues must beta sheet capping residues occur", "2")
		+ XMLSchemaAttribute::attribute_w_default("distance_cutoff_for_electrostatic_interactions", xsct_real, "Only score electrostatic interactions within this distance", "7.0")
		+ XMLSchemaAttribute::attribute_w_default("CB_b_factor_cutoff_for_electrostatic_interactions", xsct_real, "Cbeta B factor above which electrostatics will not be scored", "10000")
		+ XMLSchemaAttribute::attribute_w_default("min_primary_seq_distance_diff_electrostatic_interactions", xsct_non_negative_integer, "Minimum primary sequence distance between residues for electrostatics to be scored", "4")
		+ XMLSchemaAttribute::attribute_w_default("do_not_connect_sheets_by_loops", xsct_rosetta_bool, "Do not connect sheets with loops", "false")
		+ XMLSchemaAttribute::attribute_w_default("extract_sandwich", xsct_rosetta_bool, "Extract beta sandwiches", "true")
		+ XMLSchemaAttribute::attribute_w_default("wt_for_pro_in_starting_loop", xsct_real, "Keep wild type residue for prolines in initial loop", "20")
		+ XMLSchemaAttribute::attribute_w_default("wt_for_pro_in_1st_inter_sheet_loop", xsct_real, "Keep wild type residues for prolines in the first intersheet loop", "10")
		+ XMLSchemaAttribute::attribute_w_default("wt_for_pro_in_3rd_inter_sheet_loop", xsct_real, "Keep wild type prolines in the third intersheet loop", "1")
		+ XMLSchemaAttribute::attribute_w_default("write_all_info_files", xsct_rosetta_bool, "Write all information to files", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_AA_kind_files", xsct_rosetta_bool, "Write files specifying amino acid distribution", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_loop_AA_distribution_files", xsct_rosetta_bool, "Write files specifying loop amino acid distribution", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_strand_AA_distribution_files", xsct_rosetta_bool, "Write files specifying beta strand amino acid distribution", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_beta_sheet_capping_info", xsct_rosetta_bool, "Output information on beta sheet capping", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_chain_B_resnum", xsct_rosetta_bool, "Write chain_B_resnum file for InterfaceAnalyzer", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_electrostatic_interactions_of_all_residues", xsct_rosetta_bool, "Output info on all electrostatic interactions", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_electrostatic_interactions_of_all_residues_in_a_strand", xsct_rosetta_bool, "Output info for strand electrostatic interactions", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_electrostatic_interactions_of_surface_residues_in_a_strand", xsct_rosetta_bool, "Output surface electrostatic interactions from strands", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_heading_directions_of_all_AA_in_a_strand", xsct_rosetta_bool, "Output strand direction for all strand residues", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_p_aa_pp_files", xsct_rosetta_bool, "Output p_aa_pp score for residues", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_phi_psi_of_all", xsct_rosetta_bool, "Output all phi and psi angles", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_phi_psi_of_E", xsct_rosetta_bool, "Output strand backbone torsion angles", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_rama_at_AA_to_files", xsct_rosetta_bool, "Output all rama scores to file", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_resfile", xsct_rosetta_bool, "Output a resfile for this strand", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_resfile_NOT_FWY_on_surface", xsct_rosetta_bool, "Do not allow aromatics on the sandwich surface", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_resfile_when_seq_rec_is_bad", xsct_rosetta_bool, "Output resfile even if sequence recovery is bad", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands", xsct_rosetta_bool, "Write a resfile that will help prevent having too many internal-facing aromatics on core strands", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands", xsct_rosetta_bool, "Write a resfile that will help prevent having too many internal-facing aromatics on core strands", "false")
		+ XMLSchemaAttribute::attribute_w_default("write_resfile_to_minimize_too_much_hydrophobic_surface", xsct_rosetta_bool, "Write a resfile that will minimize surface hydrophobic residues", "false");

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Extracts and analyzes beta sandwiches.", attlist );

}

std::string SandwichFeaturesCreator::type_name() const {
	return SandwichFeatures::class_name();
}

protocols::features::FeaturesReporterOP
SandwichFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new SandwichFeatures );
}

void SandwichFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SandwichFeatures::provide_xml_schema( xsd );
}



} //namespace strand_assembly
} //namespace features
} //namespace protocols
