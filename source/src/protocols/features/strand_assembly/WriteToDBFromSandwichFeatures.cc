// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.cc
/// @brief Write to a DB after SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)

//Devel
#include <protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.hh>

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

static thread_local basic::Tracer TR( "protocols.features.strand_assembly.WriteToDBFromSandwichFeatures" );

//change_sheet_id_if_possible
bool
change_sheet_id_if_possible( // combine_sheets_if_possible
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Real min_C_O_N_angle_) {
	bool sheet_id_changed = false; // don't repeat change_sheet_id_if_possible
	Size max_sheet_id = get_max_sheet_id(struct_id, db_session);
	for ( Size i=1; i<=max_sheet_id-1; ++i ) {
		for ( Size j=i+1; j<=max_sheet_id; ++j ) {
			bool sheets_can_be_combined = see_whether_sheets_can_be_combined(
				struct_id,
				db_session,
				pose,
				i, //i_sheet
				j, //j_sheet
				min_CA_CA_dis_,
				max_CA_CA_dis_,
				min_C_O_N_angle_);

			if ( sheets_can_be_combined ) {
				if ( i<j ) {
					WriteToDB_sheet_id(
						struct_id,
						db_session,
						i, //new_sheet_id
						j); //old_sheet_id
				} else {
					WriteToDB_sheet_id(
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
} //change_sheet_id_if_possible


//delete_this_sw_can_by_sh_id_from_sw_by_comp
Size
delete_this_sw_can_by_sh_id_from_sw_by_comp(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	TR << "delete_this_sw_can_by_sh_id_from_sw_by_comp with sheet_id: " << sw_can_by_sh_id << endl;
	string select_string =
		"DELETE\t\n"
		"FROM\n"
		"\tsandwich\t\n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND (sw_can_by_sh_id = ?);";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	basic::database::safely_write_to_database(select_statement);
	return 0;
} //delete_this_sw_can_by_sh_id_from_sw_by_comp


//delete_this_struct_id
Size
delete_this_struct_id(
	StructureID struct_id,
	sessionOP db_session)
{
	TR.Info << "So exclude this sandwich " << endl;
	TR << "delete_this_struct_id with struct_id: " << struct_id << endl;
	string select_string =
		"DELETE\t\n"
		"FROM\n"
		"\tsandwich\t\n"
		"WHERE\n"
		"\t(struct_id = ?) ;";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	basic::database::safely_write_to_database(select_statement);
	return 0;
} //delete_this_struct_id


//prepare_WriteToDB_sandwich
// prepare_to_fill_sandwich
// <role> It retrieves all beta-strands of sandwich_candidate_by_sheets, it does not make sandwich_by_components
utility::vector1<SandwichFragment>
prepare_WriteToDB_sandwich(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
		"SELECT \n"
		"\tsw_sh.sw_can_by_sh_id AS sw_can_by_sh_id, \n"
		"\tsh.sheet_id AS sheet_id, \n"
		"\tsss.segment_id AS sandwich_bs_id, \n"
		"\tsss.residue_begin AS\tresidue_begin, \n"
		"\tsss.residue_end AS\tresidue_end \n"
		"FROM  \n"
		"\tsecondary_structure_segments AS sss, \n"
		"\tsheet AS sh, \n"
		"\tsw_can_by_sh AS sw_sh\n"
		"WHERE \n"
		"\t(sss.struct_id = ?) \n"
		"\tAND (sss.struct_id = sh.struct_id) \n"
		"\tAND (sss.struct_id = sw_sh.struct_id) \n"
		"\tAND (sss.dssp = 'E') \n"
		"\tAND (sw_sh.sheet_id = sh.sheet_id) \n"
		"\tAND (sh.segment_id = sss.segment_id);";


	/* (for sqlite3 try)
	SELECT
	sss.struct_id AS struct_id,
	sw_sh.sw_can_by_sh_id AS sw_can_by_sh_id,
	sh.sheet_id AS sheet_id,
	sss.segment_id AS sandwich_bs_id,
	sss.residue_begin AS residue_begin,
	sss.residue_end AS residue_end
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
	while ( res.next() )
			{
		Size sw_can_by_sh_id, sheet_id, sandwich_bs_id, residue_begin, residue_end;
		res >> sw_can_by_sh_id >> sheet_id >> sandwich_bs_id >> residue_begin >> residue_end;
		all_strands.push_back(SandwichFragment(sw_can_by_sh_id, sheet_id, sandwich_bs_id, residue_begin, residue_end));
	}
	return all_strands;
} //prepare_WriteToDB_sandwich


// WriteToDB_AA_to_terminal_loops
void
WriteToDB_AA_to_terminal_loops (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sandwich_PK_id_counter,
	Size sw_can_by_sh_id,
	string tag,
	bool starting_loop,
	Size residue_begin,
	Size residue_end)
{
	string loop_kind;

	if ( starting_loop ) {
		loop_kind = "starting_loop";
	} else { // ending_loop
		loop_kind = "ending_loop";
	}

	string insert = "INSERT INTO sandwich (struct_id, sandwich_PK_id, tag, sw_can_by_sh_id, loop_kind, component_size,\tresidue_begin, residue_end, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,?,?,?,\t?,?,?,\t?,?,\t?,?,?,?,\t?,?,?,\t?,?,?,?,?,?,?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	insert_stmt.bind(1, struct_id);
	insert_stmt.bind(2, sandwich_PK_id_counter);
	insert_stmt.bind(3, tag);
	insert_stmt.bind(4, sw_can_by_sh_id);
	insert_stmt.bind(5, loop_kind);
	Size loop_size = residue_end - residue_begin + 1;
	insert_stmt.bind(6, loop_size);
	insert_stmt.bind(7, residue_begin);
	insert_stmt.bind(8, residue_end);

	vector<Size> AA_vector = count_AA_wo_direction(dssp_pose, residue_begin, residue_end);
	insert_stmt.bind(9, AA_vector[0]); //R_num
	insert_stmt.bind(10, AA_vector[1]); //H_num
	insert_stmt.bind(11, AA_vector[2]); //K_num
	insert_stmt.bind(12, AA_vector[3]); //D
	insert_stmt.bind(13, AA_vector[4]); //E

	insert_stmt.bind(14, AA_vector[5]); //S
	insert_stmt.bind(15, AA_vector[6]); //T
	insert_stmt.bind(16, AA_vector[7]); //N
	insert_stmt.bind(17, AA_vector[8]); //Q
	insert_stmt.bind(18, AA_vector[9]); //C
	insert_stmt.bind(19, AA_vector[10]); //G
	insert_stmt.bind(20, AA_vector[11]); //P

	insert_stmt.bind(21, AA_vector[12]); //A
	insert_stmt.bind(22, AA_vector[13]); //V
	insert_stmt.bind(23, AA_vector[14]); //I
	insert_stmt.bind(24, AA_vector[15]); //L
	insert_stmt.bind(25, AA_vector[16]); //M
	insert_stmt.bind(26, AA_vector[17]); //F
	insert_stmt.bind(27, AA_vector[18]); //Y
	insert_stmt.bind(28, AA_vector[19]); //W

	basic::database::safely_write_to_database(insert_stmt);

} // WriteToDB_AA_to_terminal_loops

//WriteToDB_ending_loop
Size
WriteToDB_ending_loop (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sandwich_PK_id_counter,
	Size sw_can_by_sh_id,
	string tag,
	Size max_starting_loop_size_)
{
	string select_string =
		"SELECT\n"
		"\tmax(residue_end) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size ending_res_of_any_strand;
	while ( res.next() )
			{
		res >> ending_res_of_any_strand;
	}

	Size starting_res_of_ending_loop = 0 ;  // initial value=0 just to avoid build warning at rosetta trunk
	Size ending_res_of_ending_loop = 0;  // initial value=0 just to avoid build warning at rosetta trunk

	bool there_is_an_ending_loop = false;

	for ( Size ii = static_cast<Size>(ending_res_of_any_strand+1) ; ii <= dssp_pose.total_residue() && ii <= static_cast<Size>((static_cast<Size>(ending_res_of_any_strand) + static_cast<Size>(max_starting_loop_size_))); ii++ ) {
		//for( Size ii = (ending_res_of_any_strand+1) ; ii <= dssp_pose.total_residue() && ii <= (ending_res_of_any_strand + max_starting_loop_size_); ii++ )
		char res_ss( dssp_pose.secstruct( ii ) ) ;

		if ( res_ss == 'L' ) {
			Real dis_former_latter_AA = dssp_pose.residue(ii-1).atom("CA").xyz().distance(dssp_pose.residue(ii).atom("CA").xyz());
			//TR << "dis_former_latter_AA: " << dis_former_latter_AA << endl;
			if ( dis_former_latter_AA < 5.0 ) {
				if ( ii == ending_res_of_any_strand+1 ) {
					there_is_an_ending_loop = true;
					starting_res_of_ending_loop = ii;
				}
				ending_res_of_ending_loop = ii;
			} else {
				break;
			}
		} else {
			break;
		}
	}

	if ( !there_is_an_ending_loop ) {
		return 0;
	}

	WriteToDB_AA_to_terminal_loops (struct_id, db_session, dssp_pose, sandwich_PK_id_counter, sw_can_by_sh_id, tag, false, starting_res_of_ending_loop, ending_res_of_ending_loop);

	return 0;
} // WriteToDB_ending_loop


// WriteToDB_long_strand_id_in_each_sw
Size
WriteToDB_long_strand_id_in_each_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tresidue_begin \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(strand_edge=\'edge\' \n"
		"\tOR\tstrand_edge=\'core\') \n"
		"\tAND sw_can_by_sh_id = ? \n"
		"\tAND struct_id = ? ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, sw_can_by_sh_id);
	select_statement.bind(2, struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> vec_residue_begin_of_long_strand;
	while ( res.next() )
			{
		Size residue_begin_of_long_strand;
		res >> residue_begin_of_long_strand;
		vec_residue_begin_of_long_strand.push_back(residue_begin_of_long_strand);
	}

	for ( Size i=1; i<=vec_residue_begin_of_long_strand.size(); i++ ) {
		string update =
			"UPDATE sandwich set long_strand_id = ?\t"
			"WHERE\n"
			"\tsw_can_by_sh_id = ? \n"
			"\tAND\tresidue_begin = ? \n"
			"\tAND struct_id = ?;";

		statement update_statement(basic::database::safely_prepare_statement(update, db_session));

		update_statement.bind(1, i);
		update_statement.bind(2, sw_can_by_sh_id);
		update_statement.bind(3, vec_residue_begin_of_long_strand[i]);
		update_statement.bind(4, struct_id);

		basic::database::safely_write_to_database(update_statement);
	}

	return 0;
} // WriteToDB_long_strand_id_in_each_sw


// WriteToDB_avg_b_factor_CB_at_each_component
Size
WriteToDB_avg_b_factor_CB_at_each_component (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id)
{
	//// <begin> retrieve residue_begin, residue_end at each component
	string select_string =
		"SELECT\n"
		"\tresidue_begin, residue_end\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Size> vector_of_residue_begin;
	utility::vector1<Size> vector_of_residue_end;

	while ( res.next() )
			{
		Size residue_begin, residue_end;
		res >> residue_begin >> residue_end;
		vector_of_residue_begin.push_back(residue_begin);
		vector_of_residue_end.push_back(residue_end);
	}
	//// <end> retrieve residue_begin, residue_end at each component


	pose::PDBInfoCOP info = pose.pdb_info();
	if ( info ) {
		for ( Size i=1; i<=vector_of_residue_begin.size(); i++ ) {
			Real sum_of_b_factor_CB_at_each_component = 0;
			Size count_atoms = 0;
			for ( Size resid=vector_of_residue_begin[i]; resid<=vector_of_residue_end[i]; resid++ ) {
				Real B_factor_of_CB = info->temperature( resid, 5 ); // '5' atom will be 'H' for Gly
				sum_of_b_factor_CB_at_each_component = sum_of_b_factor_CB_at_each_component + B_factor_of_CB;
				count_atoms++;
				//for ( Size ii = 1; ii <= info->natoms( resid ); ++ii )
				//{
				// std::cout << "Temperature on " << resid << " " << ii << " " << info->temperature( resid, ii ) << std::endl;
				//}
			}
			Real avg_b_factor_CB_at_each_component = sum_of_b_factor_CB_at_each_component/count_atoms;

			// <begin> UPDATE sandwich table
			string insert =
				"UPDATE sandwich set \n"
				"avg_b_factor_CB_at_each_component = ? \n"
				"WHERE\n"
				"\t(struct_id = ?) \n"
				"\tAND\t(sw_can_by_sh_id = ?) \n"
				"\tAND\t(residue_begin = ?) ;";

			statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

			Real rounded_avg_b_factor = utility::round(avg_b_factor_CB_at_each_component);
			insert_stmt.bind(1, rounded_avg_b_factor);
			insert_stmt.bind(2, struct_id);
			insert_stmt.bind(3, sw_can_by_sh_id);
			insert_stmt.bind(4, vector_of_residue_begin[i]);

			basic::database::safely_write_to_database(insert_stmt);
			// <end> UPDATE sandwich table

		}
	}
	return 0;
} // WriteToDB_avg_b_factor_CB_at_each_component


//WriteToDB_dihedral_angle_between_core_strands_across_facing_sheets
Size
WriteToDB_dihedral_angle_between_core_strands_across_facing_sheets (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id)
{
	//// <begin> report average dihedral angles between core strands between facing sheets
	utility::vector1<Size> vec_sheet_id =  get_vec_distinct_sheet_id(struct_id, db_session, sw_can_by_sh_id);

	utility::vector1<SandwichFragment> all_strands_in_sheet_i = get_all_strands_in_sheet_i(struct_id, db_session, vec_sheet_id[1]);
	utility::vector1<SandwichFragment> all_strands_in_sheet_j = get_all_strands_in_sheet_i(struct_id, db_session, vec_sheet_id[2]);
	Real total_dihedral_angle_between_core_strands_across_facing_sheets = 0;
	Size count_dihedral_angle_between_core_strands_across_facing_sheets = 0;

	for ( Size i=1; i<=all_strands_in_sheet_i.size(); i++ ) {
		string edge = is_this_strand_at_edge_by_looking_db (struct_id, db_session, all_strands_in_sheet_i[i].get_start());
		if ( edge == "core" && all_strands_in_sheet_i[i].get_size() > 3 ) {
			for ( Size j=1; j<=all_strands_in_sheet_j.size(); j++ ) {
				edge = is_this_strand_at_edge_by_looking_db (struct_id, db_session, all_strands_in_sheet_j[j].get_start());
				if ( edge == "core" && all_strands_in_sheet_j[j].get_size() > 3 ) {
					Real dihedral = calculate_dihedral_w_4_resnums(pose, all_strands_in_sheet_i[i].get_start(), all_strands_in_sheet_i[i].get_end(),  all_strands_in_sheet_j[j].get_start(),  all_strands_in_sheet_j[j].get_end());
					//      TR << "all_strands_in_sheet_i[i].get_start(): " << all_strands_in_sheet_i[i].get_start() << endl;
					//      TR << "all_strands_in_sheet_i[i].get_end(): " << all_strands_in_sheet_i[i].get_end() << endl;
					//      TR << "all_strands_in_sheet_j[j].get_start(): " << all_strands_in_sheet_j[j].get_start() << endl;
					//      TR << "all_strands_in_sheet_j[j].get_end(): " << all_strands_in_sheet_j[j].get_end() << endl;
					//      TR << "dihedral: " << dihedral << endl;
					total_dihedral_angle_between_core_strands_across_facing_sheets = total_dihedral_angle_between_core_strands_across_facing_sheets + dihedral;
					count_dihedral_angle_between_core_strands_across_facing_sheets++;
				}
			}
		}
	}
	if ( count_dihedral_angle_between_core_strands_across_facing_sheets == 0 ) {
		return 0; // this sw is not representative to show avg dihedral angle
	}

	Real avg_dihedral_angle_between_core_strands_across_facing_sheets = total_dihedral_angle_between_core_strands_across_facing_sheets / count_dihedral_angle_between_core_strands_across_facing_sheets;
	Real rounded_dihedral = round_to_Real(avg_dihedral_angle_between_core_strands_across_facing_sheets);
	//// <end> report average dihedral angles between core strands between facing sheets


	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"avg_dihedral_angle_between_core_strands_across_facing_sheets = ? \n"
		"WHERE\n"
		"\t(sw_can_by_sh_id = ?) \n"
		"\tAND\t(struct_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	insert_stmt.bind(1, rounded_dihedral);
	insert_stmt.bind(2, sw_can_by_sh_id);
	insert_stmt.bind(3, struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;
} // WriteToDB_dihedral_angle_between_core_strands_across_facing_sheets


Size
WriteToDB_dssp_ratio_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sw_can_by_sh_id)
{
	Size H_num = 0;
	Size E_num = 0;
	Size L_num = 0;

	for ( Size ii=1; ii<=dssp_pose.total_residue(); ii++ ) {
		char res_ss( dssp_pose.secstruct( ii ) ) ;
		if ( res_ss == 'H' )  {   H_num += 1;  }
		if ( res_ss == 'E' )  {   E_num += 1;  }
		if ( res_ss == 'L' )  {   L_num += 1;  }
	}

	string update =
		"UPDATE sandwich set \n "
		"\tH_percentage = ?\t,"
		"\tE_percentage = ?\t,"
		"\tL_percentage = ?\t"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	Real rounded = round_to_Real(H_num*100/dssp_pose.total_residue());
	update_statement.bind(1, rounded);

	rounded = round_to_Real(E_num*100/dssp_pose.total_residue());
	update_statement.bind(2, rounded);

	rounded = round_to_Real(L_num*100/dssp_pose.total_residue());
	update_statement.bind(3, rounded);

	update_statement.bind(4, struct_id);
	update_statement.bind(5, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // WriteToDB_dssp_ratio_in_sw


Size
WriteToDB_hydrophobic_ratio_net_charge (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id)
{

	// <begin> sum number_of_AA
	string select_string =
		"SELECT\n"
		"\tsum(A+V+I+L+M+F+Y+W), \n"
		"\tsum(R+H+K+D+E+S+T+N+Q), \n"
		"\tsum(C+G+P), \n"
		"\tsum(R+K), \n"
		"\tsum(D+E) \n"
		"FROM\n"
		"\tsandwich\n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND (sw_can_by_sh_id = ?) ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	int number_of_hydrophobic_res, number_of_hydrophilic_res, number_of_CGP, number_of_RK_in_sw, number_of_DE_in_sw;
	while ( res.next() )
			{
		res >> number_of_hydrophobic_res >> number_of_hydrophilic_res >> number_of_CGP >> number_of_RK_in_sw >> number_of_DE_in_sw;
	}
	// <end> sum number_of_AA


	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"\tnumber_of_hydrophobic_res = ?\t,\t\n"
		"\tnumber_of_hydrophilic_res = ?\t,\t\n"
		"\tnumber_of_CGP = ?\t,\t\n"
		"\tratio_hydrophobic_philic_of_sw_in_percent = ?\t,\t\n"
		"\tnumber_of_RK_in_sw = ?\t,\t\n"
		"\tnumber_of_DE_in_sw = ?\t,\t\n"
		"\tnet_charge_of_sw = ?\t\n"
		"WHERE\n"
		"\t(sw_can_by_sh_id = ?) \n"
		"\tAND\t(struct_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	insert_stmt.bind(1, number_of_hydrophobic_res);
	insert_stmt.bind(2, number_of_hydrophilic_res);
	insert_stmt.bind(3, number_of_CGP);
	Real ratio_hydrophobic_philic_of_sw_in_percent = (number_of_hydrophobic_res*100)/(number_of_hydrophobic_res+number_of_hydrophilic_res);
	insert_stmt.bind(4, ratio_hydrophobic_philic_of_sw_in_percent);
	insert_stmt.bind(5, number_of_RK_in_sw);
	insert_stmt.bind(6, number_of_DE_in_sw);
	int net_charge_int = number_of_RK_in_sw - number_of_DE_in_sw; // Size net_charge may return like '18446744073709551612', so I don't use Size here
	// TR << "net_charge_int: " << net_charge_int << endl;
	insert_stmt.bind(7, net_charge_int); // Net charge of His at pH 7.4 is just '+0.11' according to http://www.bmolchem.wisc.edu/courses/spring503/503-sec1/503-1a.htm
	insert_stmt.bind(8, sw_can_by_sh_id);
	insert_stmt.bind(9, struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;

} // WriteToDB_hydrophobic_ratio_net_charge


//WriteToDB_min_avg_dis_between_sheets_by_cen_res
// as of 07/19/14, refactoring is not possible
Size
WriteToDB_min_avg_dis_between_sheets_by_cen_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size min_num_strands_in_sheet_)
{

	//// <begin> calculate minimum distance between sheets
	std::pair<Real, Real> minimum_average_distance_between_sheets_by_cen_res =
		cal_min_avg_dis_between_sheets_by_cen_res(
		struct_id,
		db_session,
		dssp_pose,
		all_distinct_sheet_ids,
		min_num_strands_in_sheet_);
	Real minimum_distance_between_sheets_by_cen_res = minimum_average_distance_between_sheets_by_cen_res.first;
	Real average_distance_between_sheets_by_cen_res = minimum_average_distance_between_sheets_by_cen_res.second;
	// [caution] these values are only meaningful when pdb files has only one beta-sandwich
	//// <end> calculate minimum distance between sheets

	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"\tmin_dis_between_sheets_by_cen_res = ? ,"
		"\tavg_dis_between_sheets_by_cen_res = ? "
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sw_can_by_sh_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	insert_stmt.bind(1, round_to_Real(minimum_distance_between_sheets_by_cen_res));
	insert_stmt.bind(2, round_to_Real(average_distance_between_sheets_by_cen_res));
	insert_stmt.bind(3, struct_id);
	insert_stmt.bind(4, sw_can_by_sh_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;
} // WriteToDB_min_avg_dis_between_sheets_by_cen_res


// WriteToDB_min_dis_between_sheets_by_all_res
Size
WriteToDB_min_dis_between_sheets_by_all_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids)
{

	//// <begin> calculate minimum distance between sheets
	float minimum_distance_between_sheets_by_all_res = cal_min_dis_between_sheets_by_all_res(struct_id, db_session, dssp_pose, all_distinct_sheet_ids);
	// [caution] these values are only meaningful when pdb files has only one beta-sandwich
	//// <end> calculate minimum distance between sheets

	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"\tmin_dis_between_sheets_by_all_res = ? \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sw_can_by_sh_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	insert_stmt.bind(1, round_to_Real(minimum_distance_between_sheets_by_all_res));
	//insert_stmt.bind(1, round_to_float(minimum_distance_between_sheets_by_all_res));
	insert_stmt.bind(2, struct_id);
	insert_stmt.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;
} // WriteToDB_min_dis_between_sheets_by_all_res


// WriteToDB_number_edge_strands_in_each_sw
Size
WriteToDB_number_of_edge_strands_in_each_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tcount(*) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstrand_edge=\'edge\' \n"
		"\tAND sw_can_by_sh_id = ? \n"
		"\tAND struct_id = ? ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, sw_can_by_sh_id);
	select_statement.bind(2, struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_edge_strands_in_each_sw;
	while ( res.next() )
			{
		res >> num_edge_strands_in_each_sw;
	}

	string update =
		"UPDATE sandwich set num_edge_strands_in_each_sw = ?\t"
		"WHERE\n"
		"\tsw_can_by_sh_id = ? \n"
		"\tAND struct_id = ?;";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, num_edge_strands_in_each_sw);
	update_statement.bind(2, sw_can_by_sh_id);
	update_statement.bind(3, struct_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // WriteToDB_number_edge_strands_in_each_sw


// WriteToDB_number_of_AAs_in_a_pair_of_edge_strands
Size
WriteToDB_number_of_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	sessionOP db_session,
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	Size max_num_sw_per_pdb_,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_)
{
	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii ) {
		// check_whether_sw_by_sh_id_still_alive is refactored
		bool sw_by_sh_id_still_alive = check_whether_sw_by_sh_id_still_alive (struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id());
		if ( !sw_by_sh_id_still_alive ) {
			continue;
		}

		if ( bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_ ) {
			break;
		}
		// get_sheet_antiparallel_info is refactored
		//string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());
		string strand_is_at_edge = is_this_strand_at_edge (
			pose,
			struct_id,
			db_session,
			bs_of_sw_can_by_sh[ii].get_sheet_id(),
			bs_of_sw_can_by_sh[ii].get_start(),
			bs_of_sw_can_by_sh[ii].get_end(),
			min_CA_CA_dis_,
			max_CA_CA_dis_);
		if ( strand_is_at_edge == "edge" ) {
			//get_current_bs_id_and_closest_edge_bs_id_in_different_sheet is refactored
			std::pair<Size, Size> current_bs_id_and_closest_edge_bs_id_in_different_sheet = get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (struct_id, db_session, pose, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), bs_of_sw_can_by_sh[ii].get_sheet_id(), bs_of_sw_can_by_sh[ii].get_start(), bs_of_sw_can_by_sh[ii].get_end());
			Size current_bs_id = current_bs_id_and_closest_edge_bs_id_in_different_sheet.first;
			Size closest_bs_id = current_bs_id_and_closest_edge_bs_id_in_different_sheet.second;

			WriteToDB_number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands (struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), current_bs_id, closest_bs_id);
			WriteToDB_number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands (struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(), current_bs_id, closest_bs_id);

		}
	}
	return 0;
} // WriteToDB_number_of_AAs_in_a_pair_of_edge_strands

Size
WriteToDB_number_of_core_heading_LWY_in_core_strands_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tsum(L_core_heading),\tsum(W_core_heading),\tsum(Y_core_heading) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (strand_edge = 'core') \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_core_heading_L_in_core_strands_in_sw;
	Size number_of_core_heading_W_in_core_strands_in_sw;
	Size number_of_core_heading_Y_in_core_strands_in_sw;
	while ( res.next() )
			{
		res >> number_of_core_heading_L_in_core_strands_in_sw >> number_of_core_heading_W_in_core_strands_in_sw >> number_of_core_heading_Y_in_core_strands_in_sw;
	}

	string update =
		"UPDATE sandwich set\n"
		"\tnumber_of_core_heading_L_in_core_strands_in_sw = ?\t,\n"
		"\tnumber_of_core_heading_W_in_core_strands_in_sw = ?\t,\n"
		" \tnumber_of_core_heading_Y_in_core_strands_in_sw = ?\t\n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, number_of_core_heading_L_in_core_strands_in_sw);
	update_statement.bind(2, number_of_core_heading_W_in_core_strands_in_sw);
	update_statement.bind(3, number_of_core_heading_Y_in_core_strands_in_sw);
	update_statement.bind(4, struct_id);
	update_statement.bind(5, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);
	return 0;
} // WriteToDB_number_of_core_heading_LWY_in_core_strands_in_sw


void
WriteToDB_number_of_core_heading_FWY_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tsum(F_core_heading+W_core_heading+Y_core_heading) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_core_heading_FWY_in_sw;
	while ( res.next() )
			{
		res >> number_of_core_heading_FWY_in_sw;
	}

	string update =
		"UPDATE sandwich set number_of_core_heading_FWY_in_sw = ?\t"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, number_of_core_heading_FWY_in_sw);
	update_statement.bind(2, struct_id);
	update_statement.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

} // WriteToDB_number_of_core_heading_FWY_in_sw


void
WriteToDB_number_of_core_heading_W_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tsum(W_core_heading) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_core_heading_W_in_sw;
	while ( res.next() )
			{
		res >> number_of_core_heading_W_in_sw;
	}

	string update =
		"UPDATE sandwich set number_of_core_heading_W_in_sw = ?\t"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, number_of_core_heading_W_in_sw);
	update_statement.bind(2, struct_id);
	update_statement.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

} // WriteToDB_number_of_core_heading_W_in_sw


Size
WriteToDB_number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size current_bs_id,
	Size closest_bs_id)
{
	// <begin> sum numbers of inward-pointing-AAs in current_bs_id and closest_bs_id
	string select_string =
		"SELECT\n"
		"\tsum(R_core_heading + H_core_heading + K_core_heading + D_core_heading + E_core_heading ) \n"
		"FROM\n"
		"\tsandwich\n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND (sw_can_by_sh_id = ?) \n"
		"\tAND ((sandwich_bs_id = ?) or (sandwich_bs_id = ?));";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	select_statement.bind(3, current_bs_id);
	select_statement.bind(4, closest_bs_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands;
	while ( res.next() )
			{
		res >> number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands ;
	}
	// <end> sum numbers of inward-pointing-AAs in current_bs_id and closest_bs_id

	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands = ? \n"
		"WHERE\n"
		"\t(sandwich_bs_id = ?)\t\n"
		"\tAND (sw_can_by_sh_id = ?) \n"
		"\tAND\t(struct_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	insert_stmt.bind(1, number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands);
	insert_stmt.bind(2, current_bs_id);
	insert_stmt.bind(3, sw_can_by_sh_id);
	insert_stmt.bind(4, struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;

} // WriteToDB_number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands


// WriteToDB_number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands
Size
WriteToDB_number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size current_bs_id,
	Size closest_bs_id)
{
	// <begin> sum numbers of inward-pointing-aro_AAs in current_bs_id and closest_bs_id
	string select_string =
		"SELECT\n"
		"\tsum(F_core_heading + Y_core_heading + W_core_heading) \n"
		"FROM\n"
		"\tsandwich\n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND (sw_can_by_sh_id = ?) \n"
		"\tAND ((sandwich_bs_id = ?) or (sandwich_bs_id = ?));";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	select_statement.bind(3,current_bs_id);
	select_statement.bind(4,closest_bs_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands;
	while ( res.next() )
			{
		res >> number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands ;
	}
	// <end> sum numbers of inward-pointing-AAs in current_bs_id and closest_bs_id

	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands = ? \n"
		"WHERE\n"
		"\t(sandwich_bs_id = ?)\t\n"
		"\tAND (sw_can_by_sh_id = ?) \n"
		"\tAND\t(struct_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	insert_stmt.bind(1, number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands);
	insert_stmt.bind(2, current_bs_id);
	insert_stmt.bind(3, sw_can_by_sh_id);
	insert_stmt.bind(4, struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;

} // WriteToDB_number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands


//WriteToDB_number_of_sheets_that_surround_this_sheet
Size
WriteToDB_number_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id,
	Size num_of_sheets_that_surround_this_sheet)
{
	string select_string =
		"UPDATE sheet set num_of_sheets_that_surround_this_sheet = ?\t"
		"WHERE\n"
		"\t(sheet_id = ?) \n"
		"\tAND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));

	select_statement.bind(1, num_of_sheets_that_surround_this_sheet);
	select_statement.bind(2, sheet_id);
	select_statement.bind(3, struct_id);

	basic::database::safely_write_to_database(select_statement);
	return 0;
} //WriteToDB_number_of_sheets_that_surround_this_sheet


// WriteToDB_number_strands_in_each_sw
Size
WriteToDB_number_strands_in_each_sw // it includes even 'short_edge_strands'
(StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tcount(*) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstrand_edge is not null \n"
		"\tAND sw_can_by_sh_id = ? \n"
		"\tAND struct_id = ? ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, sw_can_by_sh_id);
	select_statement.bind(2, struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size num_strands_in_each_sw;
	while ( res.next() )
			{
		res >> num_strands_in_each_sw;
	}

	string update =
		"UPDATE sandwich set num_strands_in_each_sw = ?\t"
		"WHERE\n"
		"\tsw_can_by_sh_id = ? \n"
		"\tAND struct_id = ?;";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, num_strands_in_each_sw);
	update_statement.bind(2, sw_can_by_sh_id);
	update_statement.bind(3, struct_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // WriteToDB_number_strands_in_each_sw

// WriteToDB_prolines_that_seem_to_prevent_aggregation
// reference: 2002, Sequence conservation in Ig-like domains- the role of highly conserved proline residues in the fibronectin type III superfamily
// Jane Clarke group
Size
WriteToDB_prolines_that_seem_to_prevent_aggregation (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Real wt_for_pro_in_starting_loop_,
	Real wt_for_pro_in_1st_inter_sheet_loop_,
	Real wt_for_pro_in_3rd_inter_sheet_loop_)
{
	string select_pro_in_starting_loop =
		"SELECT\n"
		"\tP \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?)\t\n"
		"\tAND (loop_kind = ?) ;";

	statement select_statement(basic::database::safely_prepare_statement(select_pro_in_starting_loop, db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	select_statement.bind(3, "starting_loop");
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_pro_in_starting_loop;
	while ( res.next() )
			{
		res >> number_of_pro_in_starting_loop;
	}


	string select_pro_in_1st_inter_sheet_loop =
		"SELECT\n"
		"\tP \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?)\t\n"
		"\tAND (loop_kind = ?) \n"
		" \tlimit 1\t;";

	statement select_statement_pro2(basic::database::safely_prepare_statement(select_pro_in_1st_inter_sheet_loop, db_session));
	select_statement_pro2.bind(1, struct_id);
	select_statement_pro2.bind(2, sw_can_by_sh_id);
	select_statement_pro2.bind(3, "loop_connecting_two_sheets");
	result res_pro2(basic::database::safely_read_from_database(select_statement_pro2));

	Size number_of_pro_in_1st_inter_sheet_loop;
	while ( res_pro2.next() )
			{
		res_pro2 >> number_of_pro_in_1st_inter_sheet_loop;
	}


	string select_pro_in_3rd_inter_sheet_loop =
		"SELECT\n"
		"\tP \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?)\t\n"
		"\tAND (loop_kind = ?) \n"
		" \tlimit 3\t;";

	statement select_statement_pro3 (basic::database::safely_prepare_statement(select_pro_in_3rd_inter_sheet_loop, db_session));
	select_statement_pro3.bind(1, struct_id);
	select_statement_pro3.bind(2, sw_can_by_sh_id);
	select_statement_pro3.bind(3, "loop_connecting_two_sheets");
	result res_pro3(basic::database::safely_read_from_database(select_statement_pro3));

	utility::vector1<Size> vec_can_number_of_pro_in_3rd_inter_sheet_loop;
	while ( res_pro3.next() )
			{
		Size num_PRO_in_inter_sheet_loop;
		res_pro3 >> num_PRO_in_inter_sheet_loop;
		vec_can_number_of_pro_in_3rd_inter_sheet_loop.push_back(num_PRO_in_inter_sheet_loop);
	}

	Size number_of_pro_in_3rd_inter_sheet_loop = vec_can_number_of_pro_in_3rd_inter_sheet_loop[3];


	string update =
		"UPDATE sandwich set \n"
		"num_PRO_in_starting_loop_and_1st_3rd_inter_sheet_loop = ? , \n"
		"num_PRO_in_starting_loop\t=\t? , \n"
		"num_PRO_in_1st_inter_sheet_loop = ? , \n"
		"num_PRO_in_3rd_inter_sheet_loop\t=\t?, \n"
		"weighted_num_PRO_prevent\t=\t? \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, (number_of_pro_in_starting_loop + number_of_pro_in_1st_inter_sheet_loop + number_of_pro_in_3rd_inter_sheet_loop));
	update_statement.bind(2, number_of_pro_in_starting_loop);
	update_statement.bind(3, number_of_pro_in_1st_inter_sheet_loop);
	update_statement.bind(4, number_of_pro_in_3rd_inter_sheet_loop);
	update_statement.bind(5, (number_of_pro_in_starting_loop*wt_for_pro_in_starting_loop_ + number_of_pro_in_1st_inter_sheet_loop*wt_for_pro_in_1st_inter_sheet_loop_ + number_of_pro_in_3rd_inter_sheet_loop*wt_for_pro_in_3rd_inter_sheet_loop_));
	update_statement.bind(6, struct_id);
	update_statement.bind(7, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // WriteToDB_prolines_that_seem_to_prevent_aggregation


// WriteToDB_ratio_of_core_heading_FWY_in_sw
Size
WriteToDB_ratio_of_core_heading_FWY_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Pose const & pose)
{
	string select_string =
		"SELECT\n"
		"\tsum(F_core_heading+W_core_heading+Y_core_heading) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size number_of_core_heading_FWY_in_sw;
	while ( res.next() )
			{
		res >> number_of_core_heading_FWY_in_sw;
	}

	Real ratio_of_core_heading_FWY_in_sw = round_to_Real((static_cast<Real>(number_of_core_heading_FWY_in_sw)*100)/static_cast<Real>(pose.total_residue()));

	string update =
		"UPDATE sandwich set ratio_of_core_heading_FWY_in_sw = ?\t"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, ratio_of_core_heading_FWY_in_sw);
	update_statement.bind(2, struct_id);
	update_statement.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // WriteToDB_ratio_of_core_heading_FWY_in_sw


//WriteToDB_rkde
Size
WriteToDB_rkde(
	StructureID struct_id,
	sessionOP db_session,
	Size rkde_PK_id_counter,
	string tag,
	Size residue_number,
	string residue_type)
{
	string insert = "INSERT INTO rkde (struct_id, rkde_PK_id, tag, residue_number, residue_type)  VALUES (?,?,?,\t?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	insert_stmt.bind(1, struct_id);
	insert_stmt.bind(2, rkde_PK_id_counter);
	insert_stmt.bind(3, tag);
	insert_stmt.bind(4, residue_number);
	insert_stmt.bind(5, residue_type);
	basic::database::safely_write_to_database(insert_stmt);

	return 0;

} //WriteToDB_rkde


// WriteToDB_rkde_in_strands
Size
WriteToDB_rkde_in_strands(
	StructureID struct_id,
	sessionOP db_session,
	Size rkde_in_strands_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size residue_number,
	string residue_type,
	string heading_direction)
{
	string insert = "INSERT INTO rkde_in_strands (struct_id, rkde_in_strands_PK_id, tag, sw_can_by_sh_id, residue_number, residue_type,\theading_direction)  VALUES (?,?,?,?,\t?,?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	insert_stmt.bind(1, struct_id);
	insert_stmt.bind(2, rkde_in_strands_PK_id_counter);
	insert_stmt.bind(3, tag);
	insert_stmt.bind(4, sw_can_by_sh_id);
	insert_stmt.bind(5, residue_number);
	insert_stmt.bind(6, residue_type);
	insert_stmt.bind(7, heading_direction); // surface or core
	basic::database::safely_write_to_database(insert_stmt);

	return 0;
} //update_rkde_in_strands

//WriteToDB_sheet
Size
WriteToDB_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_PK_id_counter,
	Size sheet_id,
	Size segment_id)
{
	string sheet_insert_i =
		"INSERT INTO sheet (sheet_PK_id, sheet_id, struct_id, segment_id)  VALUES (?,?,?,?);";
	statement sheet_insert_i_stmt(basic::database::safely_prepare_statement(sheet_insert_i,db_session));

	sheet_insert_i_stmt.bind(1, sheet_PK_id_counter);
	sheet_insert_i_stmt.bind(2, sheet_id);
	sheet_insert_i_stmt.bind(3, struct_id);
	sheet_insert_i_stmt.bind(4, segment_id);
	basic::database::safely_write_to_database(sheet_insert_i_stmt);
	return 0;
} //WriteToDB_sheet

//WriteToDB_sheet_antiparallel
Size
WriteToDB_sheet_antiparallel(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id,
	string antiparallel)
{
	string select_string =
		"UPDATE sheet set sheet_antiparallel = ?\t"
		"WHERE\n"
		"\t(sheet_id = ?) \n"
		"\tAND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));

	select_statement.bind(1,antiparallel);
	select_statement.bind(2,sheet_id);
	select_statement.bind(3,struct_id);

	basic::database::safely_write_to_database(select_statement);
	return 0;
} //WriteToDB_sheet_antiparallel

//WriteToDB_sheet_id
Size
WriteToDB_sheet_id(
	StructureID struct_id,
	sessionOP db_session,
	Size new_sheet_id,
	Size old_sheet_id)
{
	string select_string =
		"UPDATE sheet set sheet_id = ?\t"
		"WHERE\n"
		"\t(sheet_id = ?) \n"
		"\tAND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,new_sheet_id);
	select_statement.bind(2,old_sheet_id);
	select_statement.bind(3,struct_id);
	basic::database::safely_write_to_database(select_statement);
	return 0;
} // WriteToDB_sheet_id


//WriteToDB_sheet_connectivity
//update intra/inter_sheet_connection part
Size
WriteToDB_sheet_connectivity(
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	Size sandwich_PK_id_counter,
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

	if ( loop_kind == "inter_sheet" ) { // this loop connects by a inter_sheet way
		insert =
			"INSERT INTO sandwich (struct_id, sandwich_PK_id, tag, sw_can_by_sh_id, sheet_id,\tLR, canonical_LR, PA_by_preceding_E, PA_by_following_E,\tcano_PA,\theading_direction, parallel_EE, cano_parallel_EE,\tcomponent_size,\tresidue_begin, residue_end, loop_kind, inter_sheet_con_id, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,\t?,?,?,\t?,?,\t?,?,?,?,\t?,?,?,\t?,?,?,?,?,?,?,?);";
		con_id = inter_sheet_con_id;
	} else { // // this loop connects by a intra_sheet way

		sheet_id = identify_sheet_id_by_residue_end(struct_id,
			db_session,
			start_res-1); // residue_end of preceding strand

		insert =
			"INSERT INTO sandwich (struct_id, sandwich_PK_id, tag, sw_can_by_sh_id, sheet_id,\tLR, canonical_LR, PA_by_preceding_E, PA_by_following_E,\tcano_PA,\theading_direction, parallel_EE, cano_parallel_EE,\tcomponent_size,\tresidue_begin, residue_end, loop_kind, intra_sheet_con_id, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,\t?,?,?,\t?,?,\t?,?,?,?,\t?,?,?,\t?,?,?,?,?,?,?,?);";
		con_id = intra_sheet_con_id;
	}


	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	insert_stmt.bind(1, struct_id);
	insert_stmt.bind(2, sandwich_PK_id_counter);
	insert_stmt.bind(3, tag);
	insert_stmt.bind(4, sw_can_by_sh_id);
	insert_stmt.bind(5, sheet_id);
	insert_stmt.bind(6, LR);
	insert_stmt.bind(7, canonical_LR);
	insert_stmt.bind(8, PA_by_preceding_E);
	insert_stmt.bind(9, PA_by_following_E);
	insert_stmt.bind(10, cano_PA);
	insert_stmt.bind(11, heading_direction);
	insert_stmt.bind(12, parallel_EE);
	insert_stmt.bind(13, cano_parallel_EE);
	insert_stmt.bind(14, loop_size);
	insert_stmt.bind(15, start_res);
	insert_stmt.bind(16, end_res);
	insert_stmt.bind(17, loop_kind);
	insert_stmt.bind(18, con_id);

	vector<Size> AA_vector = count_AA_wo_direction(pose, start_res, end_res);
	insert_stmt.bind(19, AA_vector[0]); //R_num
	insert_stmt.bind(20, AA_vector[1]); //H_num
	insert_stmt.bind(21, AA_vector[2]); //K_num

	insert_stmt.bind(22, AA_vector[3]); //D
	insert_stmt.bind(23, AA_vector[4]); //E

	insert_stmt.bind(24, AA_vector[5]); //S
	insert_stmt.bind(25, AA_vector[6]); //T
	insert_stmt.bind(26, AA_vector[7]); //N
	insert_stmt.bind(27, AA_vector[8]); //Q

	insert_stmt.bind(28, AA_vector[9]); //C
	insert_stmt.bind(29, AA_vector[10]); //G
	insert_stmt.bind(30, AA_vector[11]); //P

	insert_stmt.bind(31, AA_vector[12]); //A
	insert_stmt.bind(32, AA_vector[13]); //V
	insert_stmt.bind(33, AA_vector[14]); //I
	insert_stmt.bind(34, AA_vector[15]); //L
	insert_stmt.bind(35, AA_vector[16]); //M
	insert_stmt.bind(36, AA_vector[17]); //F
	insert_stmt.bind(37, AA_vector[18]); //Y
	insert_stmt.bind(38, AA_vector[19]); //W

	basic::database::safely_write_to_database(insert_stmt);

	return 0;
} //WriteToDB_sheet_connectivity


//WriteToDB_shortest_dis_between_facing_aro_in_sw
Real
WriteToDB_shortest_dis_between_facing_aro_in_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Pose const & pose,
	// Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size min_num_strands_in_sheet_)
{

	//// <begin> calculate minimum distance between sheets
	float shortest_dis_between_facing_aro_in_sw =
		cal_shortest_dis_between_facing_aro_in_sw(
		struct_id,
		db_session,
		pose,
		all_distinct_sheet_ids,
		min_num_strands_in_sheet_);
	//// <end> calculate minimum distance between sheets

	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"\tshortest_dis_between_facing_aro_in_sw = ? "
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sw_can_by_sh_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));
	if ( (shortest_dis_between_facing_aro_in_sw == 999.0) || (shortest_dis_between_facing_aro_in_sw == 9999.0) ) {
		TR << "this sandwich either lacks aromatic residue in at least 1 sheet or has > 2 sheets, so don't report shortest_dis_between_facing_aro_in_sw" << endl;
		return 0;
	}

	insert_stmt.bind(1, round_to_Real(shortest_dis_between_facing_aro_in_sw));
	insert_stmt.bind(2, struct_id);
	insert_stmt.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return round_to_Real(shortest_dis_between_facing_aro_in_sw);
} // WriteToDB_shortest_dis_between_facing_aro_in_sw


//WriteToDB_starting_loop
Size
WriteToDB_starting_loop (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sandwich_PK_id_counter,
	Size sw_can_by_sh_id,
	string tag,
	Size max_starting_loop_size_)
{
	string select_string =
		"SELECT\n"
		"\tmin(residue_begin) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	int starting_res_of_any_strand;
	while ( res.next() )
			{
		res >> starting_res_of_any_strand;
	}

	Size starting_res_of_starting_loop = 0; // initial value=0 just to avoid build warning at rosetta trunk
	Size ending_res_of_starting_loop = 0 ;  // initial value=0 just to avoid build warning at rosetta trunk

	bool there_is_a_starting_loop = false;

	// I used to use Size type for i, starting_res_of_any_strand and static_cast<Size> for max_starting_loop_size_, but as of 05/02/2013, simple math like 'static_cast<Size>(starting_res_of_any_strand) - static_cast<Size>(max_starting_loop_size_)' doesn't work, so I change to use int
	for ( int ii = starting_res_of_any_strand-1;
			(ii >= 1) && (ii >= (starting_res_of_any_strand) - (static_cast<int>(max_starting_loop_size_) ));
			ii-- ) {
		char res_ss( dssp_pose.secstruct( ii ) ) ;

		if ( res_ss == 'L' ) {
			Real dis_former_latter_AA = dssp_pose.residue(ii+1).atom("CA").xyz().distance(dssp_pose.residue(ii).atom("CA").xyz());
			if ( dis_former_latter_AA < 5.0 ) {
				if ( ii == starting_res_of_any_strand-1 ) {
					there_is_a_starting_loop = true;
					ending_res_of_starting_loop = ii;
				}
				starting_res_of_starting_loop = ii;
			} else {
				break;
			}
		} else {
			break;
		}
	}

	if ( !there_is_a_starting_loop ) {
		return 0;
	}

	WriteToDB_AA_to_terminal_loops (struct_id, db_session, dssp_pose, sandwich_PK_id_counter, sw_can_by_sh_id, tag, true, starting_res_of_starting_loop, ending_res_of_starting_loop);

	return 0;
} // WriteToDB_starting_loop

// WriteToDB_sw_res_size
Size
WriteToDB_sw_res_size (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
		"SELECT\n"
		"\tmin(residue_begin), max(residue_end) \n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size starting_res_of_sw;
	Size ending_res_of_sw;
	while ( res.next() )
			{
		res >> starting_res_of_sw >> ending_res_of_sw;
	}
	Size sw_res_size = ending_res_of_sw - starting_res_of_sw + 1;

	string update =
		"UPDATE sandwich set sw_res_size = ?\t"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1, sw_res_size);
	update_statement.bind(2, struct_id);
	update_statement.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

	return sw_res_size;
} // WriteToDB_sw_res_size



// (begin) Run_WriteToDB_sandwich
Size
Run_WriteToDB_sandwich(
	string tag,
	Pose & dssp_pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	Size max_num_sw_per_pdb_,
	StructureID struct_id, // needed argument
	sessionOP db_session,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Size sandwich_PK_id_counter
)
{
	for ( Size ii=1; ii<=bs_of_sw_can_by_sh.size(); ++ii ) {
		if ( bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id() > max_num_sw_per_pdb_ ) {
			break;
		}
		string sheet_antiparallel = get_sheet_antiparallel_info(struct_id, db_session, bs_of_sw_can_by_sh[ii].get_sheet_id());

		string strand_is_at_edge = is_this_strand_at_edge (
			dssp_pose,
			struct_id,
			db_session,
			bs_of_sw_can_by_sh[ii].get_sheet_id(),
			bs_of_sw_can_by_sh[ii].get_start(),
			bs_of_sw_can_by_sh[ii].get_end(),
			min_CA_CA_dis_,
			max_CA_CA_dis_);

		Size component_size = bs_of_sw_can_by_sh[ii].get_size();
		WriteToDB_sandwich (struct_id,
			db_session,
			dssp_pose,
			sandwich_PK_id_counter,
			tag,
			bs_of_sw_can_by_sh[ii].get_sw_can_by_sh_id(),
			bs_of_sw_can_by_sh[ii].get_sheet_id(),
			sheet_antiparallel,
			bs_of_sw_can_by_sh[ii].get_strand_id(),
			strand_is_at_edge,
			component_size,
			bs_of_sw_can_by_sh[ii].get_start(),
			bs_of_sw_can_by_sh[ii].get_end());
		sandwich_PK_id_counter++;
	}
	return sandwich_PK_id_counter;
} // Run_WriteToDB_sandwich


//WriteToDB_sandwich
Size
WriteToDB_sandwich (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sandwich_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size sheet_id,
	string sheet_antiparallel,
	Size sandwich_bs_id,
	string strand_is_at_edge,
	Size component_size,
	Size residue_begin,
	Size residue_end)
{
	string insert =
		"INSERT INTO sandwich (struct_id, sandwich_PK_id, tag, sw_can_by_sh_id, sheet_id, sheet_antiparallel, sandwich_bs_id, strand_edge, component_size, residue_begin, residue_end, R,H,K, D,E, S,T,N,Q, C,G,P, A,V,I,L,M,F,Y,W)  VALUES (?,?,?,?,?,\t?,?,?,?,?,\t?,\t?,?,?,\t?,?,\t?,?,?,?,\t?,?,?,\t?,?,?,?,?,?,?,?);";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	insert_stmt.bind(1, struct_id);
	insert_stmt.bind(2, sandwich_PK_id_counter);
	insert_stmt.bind(3, tag);
	insert_stmt.bind(4, sw_can_by_sh_id);
	insert_stmt.bind(5, sheet_id);
	insert_stmt.bind(6, sheet_antiparallel);
	insert_stmt.bind(7, sandwich_bs_id); //bs_id
	insert_stmt.bind(8, strand_is_at_edge);
	insert_stmt.bind(9, component_size);
	insert_stmt.bind(10, residue_begin);
	insert_stmt.bind(11, residue_end);
	vector<Size> AA_vector = count_AA_wo_direction(pose,  residue_begin, residue_end);
	insert_stmt.bind(12, AA_vector[0]); //R_num
	insert_stmt.bind(13, AA_vector[1]); //H_num
	insert_stmt.bind(14, AA_vector[2]); //K_num
	insert_stmt.bind(15, AA_vector[3]); //D
	insert_stmt.bind(16, AA_vector[4]); //E

	insert_stmt.bind(17, AA_vector[5]); //S
	insert_stmt.bind(18, AA_vector[6]); //T
	insert_stmt.bind(19, AA_vector[7]); //N
	insert_stmt.bind(20, AA_vector[8]); //Q
	insert_stmt.bind(21, AA_vector[9]); //C
	insert_stmt.bind(22, AA_vector[10]); //G
	insert_stmt.bind(23, AA_vector[11]); //P

	insert_stmt.bind(24, AA_vector[12]); //A
	insert_stmt.bind(25, AA_vector[13]); //V
	insert_stmt.bind(26, AA_vector[14]); //I
	insert_stmt.bind(27, AA_vector[15]); //L
	insert_stmt.bind(28, AA_vector[16]); //M
	insert_stmt.bind(29, AA_vector[17]); //F
	insert_stmt.bind(30, AA_vector[18]); //Y
	insert_stmt.bind(31, AA_vector[19]); //W

	basic::database::safely_write_to_database(insert_stmt);
	return 0;
} //WriteToDB_sandwich


// WriteToDB_sandwich_by_AA_w_direction
Size
WriteToDB_sandwich_by_AA_w_direction (
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
		"UPDATE sandwich set \n"
		"R_core_heading = ? , \n"
		"R_surface_heading\t=\t? , \n"
		"H_core_heading = ? , \n"
		"H_surface_heading\t=\t? , \n"
		"K_core_heading = ? , \n"
		"K_surface_heading\t=\t? , \n"
		"D_core_heading = ? , \n"
		"D_surface_heading\t=\t? , \n"
		"E_core_heading = ? , \n"
		"E_surface_heading\t=\t? , \n"
		"S_core_heading = ? , \n"
		"S_surface_heading\t=\t? , \n"
		"T_core_heading = ? , \n"
		"T_surface_heading\t=\t? , \n"
		"N_core_heading = ? , \n"
		"N_surface_heading\t=\t? , \n"
		"Q_core_heading = ? , \n"
		"Q_surface_heading\t=\t? , \n"
		"C_core_heading = ? , \n"
		"C_surface_heading\t=\t? , \n"
		"G_core_heading = ? , \n"
		"G_surface_heading\t=\t? , \n"
		"P_core_heading = ? , \n"
		"P_surface_heading\t=\t? , \n"
		"A_core_heading = ? , \n"
		"A_surface_heading\t=\t? , \n"
		"V_core_heading\t= ? , \n"
		"V_surface_heading\t=\t? ,\t\n"
		"I_core_heading\t= ? , \n"
		"I_surface_heading\t=\t? , \n"
		"L_core_heading\t= ? , \n"
		"L_surface_heading\t=\t? , \n"
		"M_core_heading = ? , \n"
		"M_surface_heading\t=\t? , \n"
		"F_core_heading\t= ? , \n"
		"F_surface_heading\t=\t? ,\t\n"
		"Y_core_heading\t= ? , \n"
		"Y_surface_heading\t=\t? , \n"
		"W_core_heading\t= ? , \n"
		"W_surface_heading\t=\t? \n"
		"WHERE\n"
		"\tresidue_begin = ?\t\n"
		"\tAND\tstruct_id = ?;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	vector<Size> AA_vector = count_AA_w_direction(struct_id, db_session, pose, pose_w_center_000, sw_can_by_sh_id, sheet_id, residue_begin, residue_end);

	insert_stmt.bind(1, AA_vector[0]); // R_core_heading
	insert_stmt.bind(2, AA_vector[1]); // R_surface_heading
	insert_stmt.bind(3, AA_vector[2]); // H_core_heading
	insert_stmt.bind(4, AA_vector[3]); //
	insert_stmt.bind(5, AA_vector[4]); // K_core_heading
	insert_stmt.bind(6, AA_vector[5]); //
	insert_stmt.bind(7, AA_vector[6]); // D_core_heading
	insert_stmt.bind(8, AA_vector[7]); //
	insert_stmt.bind(9, AA_vector[8]); // E_core_heading
	insert_stmt.bind(10, AA_vector[9]); //

	insert_stmt.bind(11, AA_vector[10]); // S_core_heading
	insert_stmt.bind(12, AA_vector[11]); // S_surface_heading
	insert_stmt.bind(13, AA_vector[12]); // T_core_heading
	insert_stmt.bind(14, AA_vector[13]); // T_surface_heading
	insert_stmt.bind(15, AA_vector[14]); // N_core_heading
	insert_stmt.bind(16, AA_vector[15]); //
	insert_stmt.bind(17, AA_vector[16]); // Q_core_heading
	insert_stmt.bind(18, AA_vector[17]); //

	insert_stmt.bind(19, AA_vector[18]); // C_core_heading
	insert_stmt.bind(20, AA_vector[19]); // C_surface_heading
	insert_stmt.bind(21, AA_vector[20]); // G_core_heading
	insert_stmt.bind(22, AA_vector[21]); //
	insert_stmt.bind(23, AA_vector[22]); // P_core_heading
	insert_stmt.bind(24, AA_vector[23]); //

	insert_stmt.bind(25, AA_vector[24]); // A_core_heading
	insert_stmt.bind(26, AA_vector[25]); // A_surface_heading
	insert_stmt.bind(27, AA_vector[26]); // V_core_heading
	insert_stmt.bind(28, AA_vector[27]); // V_surface_heading
	insert_stmt.bind(29, AA_vector[28]); // I_core_heading
	insert_stmt.bind(30, AA_vector[29]); //
	insert_stmt.bind(31, AA_vector[30]); // L_core_heading
	insert_stmt.bind(32, AA_vector[31]); //

	insert_stmt.bind(33, AA_vector[32]); // M_core_heading
	insert_stmt.bind(34, AA_vector[33]); // M_surface_heading
	insert_stmt.bind(35, AA_vector[34]); // F_core_heading
	insert_stmt.bind(36, AA_vector[35]); // F_surface_heading
	insert_stmt.bind(37, AA_vector[36]); // Y_core_heading
	insert_stmt.bind(38, AA_vector[37]); //
	insert_stmt.bind(39, AA_vector[38]); // W_core_heading
	insert_stmt.bind(40, AA_vector[39]); //

	insert_stmt.bind(41, residue_begin);
	insert_stmt.bind(42, struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	return 0;
} // WriteToDB_sandwich_by_AA_w_direction

//WriteToDB_sw_can_by_sh
Size
WriteToDB_sw_can_by_sh (
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

	statement sw_can_by_sh_insert_stmt(basic::database::safely_prepare_statement(sw_can_by_sh_insert, db_session));

	sw_can_by_sh_insert_stmt.bind(1, struct_id);
	sw_can_by_sh_insert_stmt.bind(2, sw_can_by_sh_PK_id_counter);
	sw_can_by_sh_insert_stmt.bind(3, tag);
	sw_can_by_sh_insert_stmt.bind(4, sw_can_by_sh_id_counter);
	sw_can_by_sh_insert_stmt.bind(5, sheet_id);
	sw_can_by_sh_insert_stmt.bind(6, num_strands_from_sheet); // number of strands of sheet
	basic::database::safely_write_to_database(sw_can_by_sh_insert_stmt);
	return 0;
} // WriteToDB_sw_can_by_sh


Size
WriteToDB_topology_candidate (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id)
{
	// Warning: this is NOT a strict fnIII topology determinator, this function is useful only to identify fnIII topology beta-sandwich from so many beta-sandwiches.
	// So final human inspection is required to confirm fnIII eventually after this function. So column name is 'topology_candidate' instead of 'topology'
	// Also, be sure to run with individual sandwich in each pdb file to better identify fnIII because sheet_id is not necessarily final one when each pdb file has multiple sandwiches (starting 06/20/2014)
	// (starting 06/20/2014), all 103 fnIII identified out of 2,805 WT beta-sandwiches are confirmed to be correct even by manual inspection
	// Changed from "fn3" to "fnIII" because "fnIII" seems more correct according to "Manipulating the stability of fibronectin type III domains by protein engineering, Sean P Ng, K S Billings, L G Randles and Jane Clarke, Nanotechnology 19 (2008) 384023" (on 06/13/14)


	//// <begin> retrieve long_strand_id
	string select_string_all =
		"SELECT\n"
		"\tlong_strand_id\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(strand_edge=\'edge\' \n"
		"\tOR\tstrand_edge=\'core\') \n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_all(basic::database::safely_prepare_statement(select_string_all, db_session));
	select_statement_all.bind(1, struct_id);
	select_statement_all.bind(2, sw_can_by_sh_id);
	result res_all(basic::database::safely_read_from_database(select_statement_all));

	utility::vector1<int> vector_of_long_strand_id;

	while ( res_all.next() )
			{
		int long_strand_id;
		res_all >> long_strand_id;
		vector_of_long_strand_id.push_back(long_strand_id);
	}
	//// <end> retrieve long_strand_id


	//// <begin> retrieve long_strand_id from sheet_id = 1
	string select_string =
		"SELECT\n"
		"\tlong_strand_id\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sheet_id = 1)\t\n"
		"\tAND\t(strand_edge=\'edge\' \n"
		"\tOR\tstrand_edge=\'core\') \n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string, db_session));
	select_statement.bind(1, struct_id);
	select_statement.bind(2, sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<int> vector_of_long_strand_id_from_sheet_1;

	while ( res.next() )
			{
		int long_strand_id_from_sheet_1;
		res >> long_strand_id_from_sheet_1;
		vector_of_long_strand_id_from_sheet_1.push_back(long_strand_id_from_sheet_1);
	}
	//// <end> retrieve long_strand_id from sheet_id = 1


	string topology_candidate;

	//// <begin> check whether there is 'P_or_mix' sheet like 1F56
	string select_string_sheet =
		"SELECT\n"
		"\tcount(*)\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sheet_antiparallel = 'P_or_mix')\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_sheet(basic::database::safely_prepare_statement(select_string_sheet, db_session));
	select_statement_sheet.bind(1, struct_id);
	select_statement_sheet.bind(2, sw_can_by_sh_id);
	result res_sheet(basic::database::safely_read_from_database(select_statement_sheet));

	while ( res_sheet.next() )
			{
		Size num_P_or_mix;
		res_sheet >> num_P_or_mix;
		if ( num_P_or_mix != 0 ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check whether there is 'P_or_mix' sheet like 1F56


	//// <begin> check long_strand_id = 1
	string select_string_1 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 1)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_1(basic::database::safely_prepare_statement(select_string_1, db_session));
	select_statement_1.bind(1, struct_id);
	select_statement_1.bind(2, sw_can_by_sh_id);
	result res_1(basic::database::safely_read_from_database(select_statement_1));

	while ( res_1.next() )
			{
		string strand_edge;
		res_1 >> strand_edge;
		if ( strand_edge != "edge" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 1


	//// <begin> check long_strand_id = 2
	string select_string_2 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 2)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_2(basic::database::safely_prepare_statement(select_string_2, db_session));
	select_statement_2.bind(1, struct_id);
	select_statement_2.bind(2, sw_can_by_sh_id);
	result res_2(basic::database::safely_read_from_database(select_statement_2));

	while ( res_2.next() )
			{
		string strand_edge;
		res_2 >> strand_edge;
		if ( strand_edge != "core" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 2


	//// <begin> check long_strand_id = 3
	string select_string_3 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 3)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_3(basic::database::safely_prepare_statement(select_string_3, db_session));
	select_statement_3.bind(1, struct_id);
	select_statement_3.bind(2, sw_can_by_sh_id);
	result res_3(basic::database::safely_read_from_database(select_statement_3));

	while ( res_3.next() )
			{
		string strand_edge;
		res_3 >> strand_edge;
		if ( strand_edge != "core" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 3

	//// <begin> check long_strand_id = 4
	string select_string_4 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 4)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_4(basic::database::safely_prepare_statement(select_string_4, db_session));
	select_statement_4.bind(1, struct_id);
	select_statement_4.bind(2, sw_can_by_sh_id);
	result res_4(basic::database::safely_read_from_database(select_statement_4));

	while ( res_4.next() )
			{
		string strand_edge;
		res_4 >> strand_edge;
		if ( strand_edge != "edge" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 4


	//// <begin> check long_strand_id = 5
	string select_string_5 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 5)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_5(basic::database::safely_prepare_statement(select_string_5, db_session));
	select_statement_5.bind(1, struct_id);
	select_statement_5.bind(2, sw_can_by_sh_id);
	result res_5(basic::database::safely_read_from_database(select_statement_5));

	while ( res_5.next() )
			{
		string strand_edge;
		res_5 >> strand_edge;
		if ( strand_edge != "edge" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 5


	//// <begin> check long_strand_id = 6
	string select_string_6 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 6)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_6(basic::database::safely_prepare_statement(select_string_6, db_session));
	select_statement_6.bind(1, struct_id);
	select_statement_6.bind(2, sw_can_by_sh_id);
	result res_6(basic::database::safely_read_from_database(select_statement_6));

	while ( res_6.next() )
			{
		string strand_edge;
		res_6 >> strand_edge;
		if ( strand_edge != "core" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 6


	//// <begin> check long_strand_id = 7
	string select_string_7 =
		"SELECT\n"
		"\tstrand_edge\t\n"
		"FROM\n"
		"\tsandwich \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(long_strand_id = 7)\t\n"
		"\tAND\t(sw_can_by_sh_id = ?);";

	statement select_statement_7(basic::database::safely_prepare_statement(select_string_7, db_session));
	select_statement_7.bind(1, struct_id);
	select_statement_7.bind(2, sw_can_by_sh_id);
	result res_7(basic::database::safely_read_from_database(select_statement_7));

	while ( res_7.next() )
			{
		string strand_edge;
		res_7 >> strand_edge;
		if ( strand_edge != "edge" ) {
			topology_candidate = "not_fnIII";
		}
	}
	//// <end> check long_strand_id = 7


	//TR << "vector_of_long_strand_id_from_sheet_1.size(): " << vector_of_long_strand_id_from_sheet_1.size() << std::endl;
	std::cout << vector_of_long_strand_id_from_sheet_1.size() << std::endl;

	if ( (vector_of_long_strand_id.size() != 7) || (vector_of_long_strand_id_from_sheet_1.size() != 3) ) {
		topology_candidate = "not_known_topology";
	} else {
		if ( topology_candidate == "not_fnIII" ) {} // this sandwich didn't pass fnIII check based on 'edge'/'core' check
		else if ( vector_of_long_strand_id_from_sheet_1[1] == 1 && vector_of_long_strand_id_from_sheet_1[2] == 2 && vector_of_long_strand_id_from_sheet_1[3] == 5 ) {
			topology_candidate = "fnIII";
		} else {
			topology_candidate = "not_fnIII";
		}
	}

	// <begin> UPDATE sandwich table
	string insert =
		"UPDATE sandwich set \n"
		"topology_candidate = ? \n"
		"WHERE\n"
		"\t(struct_id = ?) \n"
		"\tAND\t(sw_can_by_sh_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert, db_session));

	insert_stmt.bind(1, topology_candidate);
	insert_stmt.bind(2, struct_id);
	insert_stmt.bind(3, sw_can_by_sh_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sandwich table

	return 0;
} // WriteToDB_topology_candidate


void
WriteToDB_turn_AA(
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size i,
	StructureID struct_id,
	sessionOP db_session,
	string turn_type)
{
	string canonical_turn_AA = "F_canonical_turn_AA";
	if ( turn_type == "I" ) {
		if ( pose.residue_type(i).name3() == "LEU" || pose.residue_type(i).name3() == "ALA" || pose.residue_type(i).name3() == "GLY" || pose.residue_type(i).name3() == "PRO" ||
				pose.residue_type(i).name3() == "THR" || pose.residue_type(i).name3() == "SER" || pose.residue_type(i).name3() == "GLU" || pose.residue_type(i).name3() == "ASN" ||
				pose.residue_type(i).name3() == "ASP" ) {
			if ( pose.residue_type(i+1).name3() == "LEU" || pose.residue_type(i+1).name3() == "ALA" || pose.residue_type(i+1).name3() == "PRO" ||
					pose.residue_type(i+1).name3() == "THR" || pose.residue_type(i+1).name3() == "SER" || pose.residue_type(i+1).name3() == "GLU" ||
					pose.residue_type(i+1).name3() == "ASP" || pose.residue_type(i+1).name3() == "LYS" ) {
				if ( pose.residue_type(i+2).name3() == "ALA" || pose.residue_type(i+2).name3() == "GLY" || pose.residue_type(i+2).name3() == "THR" ||
						pose.residue_type(i+2).name3() == "SER" || pose.residue_type(i+2).name3() == "GLU" || pose.residue_type(i+2).name3() == "ASN" ||
						pose.residue_type(i+2).name3() == "ASP" || pose.residue_type(i+2).name3() == "LYS" ) {
					if ( pose.residue_type(i+3).name3() == "VAL" || pose.residue_type(i+3).name3() == "LEU" || pose.residue_type(i+3).name3() == "ALA" ||
							pose.residue_type(i+3).name3() == "GLY" || pose.residue_type(i+3).name3() == "THR" || pose.residue_type(i+3).name3() == "SER" ||
							pose.residue_type(i+3).name3() == "GLU" || pose.residue_type(i+3).name3() == "ASP" || pose.residue_type(i+3).name3() == "LYS" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else if ( turn_type == "II" ) {
		if ( pose.residue_type(i).name3() == "VAL" || pose.residue_type(i).name3() == "LEU" || pose.residue_type(i).name3() == "ALA" || pose.residue_type(i).name3() == "GLY" ||
				pose.residue_type(i).name3() == "TYR" || pose.residue_type(i).name3() == "PRO" || pose.residue_type(i).name3() == "GLU" || pose.residue_type(i).name3() == "LYS" ) {
			if ( pose.residue_type(i+1).name3() == "ALA" || pose.residue_type(i+1).name3() == "PRO" || pose.residue_type(i+1).name3() == "SER" ||
					pose.residue_type(i+1).name3() == "GLU" || pose.residue_type(i+1).name3() == "LYS" ) {
				if ( pose.residue_type(i+2).name3() == "GLY" ) {
					if ( pose.residue_type(i+3).name3() == "VAL" || pose.residue_type(i+3).name3() == "ALA" || pose.residue_type(i+3).name3() == "SER" ||
							pose.residue_type(i+3).name3() == "GLU" || pose.residue_type(i+3).name3() == "LYS" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else if ( turn_type == "VIII" ) {
		if ( pose.residue_type(i).name3() == "GLY" || pose.residue_type(i).name3() == "PRO" ) {
			if ( pose.residue_type(i+1).name3() == "PRO" || pose.residue_type(i+1).name3() == "ASP" ) {
				if ( pose.residue_type(i+2).name3() == "VAL" || pose.residue_type(i+2).name3() == "LEU" || pose.residue_type(i+2).name3() == "ASN" || pose.residue_type(i+2).name3() == "ASP" ) {
					if ( pose.residue_type(i+3).name3() == "PRO" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else if ( turn_type == "I_prime" ) {
		if ( pose.residue_type(i).name3() == "ILE" || pose.residue_type(i).name3() == "VAL" || pose.residue_type(i).name3() == "LEU" ||
				pose.residue_type(i).name3() == "ALA" || pose.residue_type(i).name3() == "TYR" || pose.residue_type(i).name3() == "THR" ||
				pose.residue_type(i).name3() == "SER" || pose.residue_type(i).name3() == "ASP" || pose.residue_type(i).name3() == "LYS" ) {
			if ( pose.residue_type(i+1).name3() == "PRO" || pose.residue_type(i+1).name3() == "GLY"  || pose.residue_type(i+1).name3() == "HIS" ||
					pose.residue_type(i+1).name3() == "ASN" || pose.residue_type(i+1).name3() == "ASP" ) {
				if ( pose.residue_type(i+2).name3() == "GLY" ) {
					if ( pose.residue_type(i+3).name3() == "VAL" || pose.residue_type(i+3).name3() == "GLU" || pose.residue_type(i+3).name3() == "ASN" ||
							pose.residue_type(i+3).name3() == "LYS" || pose.residue_type(i+3).name3() == "ARG" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else if ( turn_type == "II_prime" ) {
		// TR << "pose.residue_type(i).name3(): " << pose.residue_type(i).name3() << endl;
		if ( pose.residue_type(i).name3() == "PHE" || pose.residue_type(i).name3() == "VAL" || pose.residue_type(i).name3() == "LEU" ||
				pose.residue_type(i).name3() == "ALA" || pose.residue_type(i).name3() == "GLY" ||
				pose.residue_type(i).name3() == "TYR" || pose.residue_type(i).name3() == "THR" ||
				pose.residue_type(i).name3() == "SER" || pose.residue_type(i).name3() == "HIS" || pose.residue_type(i).name3() == "GLU" ||
				pose.residue_type(i).name3() == "ASN" ||
				pose.residue_type(i).name3() == "GLN" || pose.residue_type(i).name3() == "ASP" || pose.residue_type(i).name3() == "ARG" ) {
			// TR << "pose.residue_type(i+1).name3(): " << pose.residue_type(i+1).name3() << endl;
			if ( pose.residue_type(i+1).name3() == "GLY" ) {
				// TR << "pose.residue_type(i+2).name3(): " << pose.residue_type(i+2).name3() << endl;
				if ( pose.residue_type(i+2).name3() == "LEU" || pose.residue_type(i+2).name3() == "ALA" || pose.residue_type(i+2).name3() == "GLY" ||
						pose.residue_type(i+2).name3() == "PRO" ||
						pose.residue_type(i+2).name3() == "SER" || pose.residue_type(i+2).name3() == "GLU" ||
						pose.residue_type(i+2).name3() == "ASN" || pose.residue_type(i+2).name3() == "ASP" || pose.residue_type(i+2).name3() == "LYS" ) {
					// TR << "pose.residue_type(i+3).name3(): " << pose.residue_type(i+3).name3() << endl;
					if ( pose.residue_type(i+3).name3() == "PHE" || pose.residue_type(i+3).name3() == "VAL" || pose.residue_type(i+3).name3() == "LEU" ||
							pose.residue_type(i+3).name3() == "ALA" ||
							pose.residue_type(i+3).name3() == "GLY" || pose.residue_type(i+3).name3() == "TYR" || pose.residue_type(i+3).name3() == "THR" ||
							pose.residue_type(i+3).name3() == "SER" || pose.residue_type(i+3).name3() == "GLU" ||
							pose.residue_type(i+3).name3() == "ASN" || pose.residue_type(i+3).name3() == "GLN" || pose.residue_type(i+3).name3() == "LYS" || pose.residue_type(i+3).name3() == "ARG" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else if ( turn_type == "VIa1" || turn_type == "VIa2" ) {
		if ( pose.residue_type(i).name3() == "PHE" || pose.residue_type(i).name3() == "VAL" || pose.residue_type(i).name3() == "THR" || pose.residue_type(i).name3() == "HIS" ||
				pose.residue_type(i).name3() == "ASN" ) {
			if ( pose.residue_type(i+1).name3() == "ILE" || pose.residue_type(i+1).name3() == "SER" || pose.residue_type(i+1).name3() == "ASN" ) {
				if ( pose.residue_type(i+2).name3() == "PRO" ) {
					if ( pose.residue_type(i+3).name3() == "GLY" || pose.residue_type(i+3).name3() == "THR" || pose.residue_type(i+3).name3() == "HIS" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else if ( turn_type == "VIb" ) {
		if ( pose.residue_type(i).name3() == "PHE" || pose.residue_type(i).name3() == "GLY" || pose.residue_type(i).name3() == "THR" || pose.residue_type(i).name3() == "SER" ) {
			if ( pose.residue_type(i+1).name3() == "LEU" || pose.residue_type(i+1).name3() == "TYR" || pose.residue_type(i+1).name3() == "THR"  || pose.residue_type(i+1).name3() == "GLU" ) {
				if ( pose.residue_type(i+2).name3() == "PRO" ) {
					if ( pose.residue_type(i+3).name3() == "PHE" || pose.residue_type(i+3).name3() == "ALA" || pose.residue_type(i+3).name3() == "TYR" ||
							pose.residue_type(i+3).name3() == "THR" || pose.residue_type(i+3).name3() == "LYS" ) {
						canonical_turn_AA = "T_canonical_turn_AA";
					}
				}

			}
		}
	} else { //(turn_type == 'IV')
		canonical_turn_AA = "uncertain_canonical_turn_AA_since_turn_type_eq_IV";
	}


	string select_string =
		"UPDATE sandwich set \n"
		"i_AA = ? , \n"
		"i_p1_AA\t=\t? , \n"
		"i_p2_AA\t=\t? , \n"
		"i_p3_AA\t=\t? ,\t\n"
		"canonical_turn_AA\t=\t? \n"
		"WHERE\n"
		"\t(sw_can_by_sh_id = ?) \n"
		"\tAND\t(residue_begin = ?) \n"
		"\tAND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, pose.residue_type(i).name3());
	select_statement.bind(2, pose.residue_type(i+1).name3());
	select_statement.bind(3, pose.residue_type(i+2).name3());
	select_statement.bind(4, pose.residue_type(i+3).name3());
	select_statement.bind(5, canonical_turn_AA);
	select_statement.bind(6, sw_can_by_sh_id);
	select_statement.bind(7, (i+1));
	select_statement.bind(8, struct_id);
	basic::database::safely_write_to_database(select_statement);

} //WriteToDB_turn_AA


string
WriteToDB_turn_type(
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size first_res,
	Size second_res,
	StructureID struct_id,
	sessionOP db_session,
	Real allowed_deviation_for_turn_type_id_)
{
	Real first_phi = pose.phi(first_res);
	Real first_psi = pose.psi(first_res);

	Real second_phi = pose.phi(second_res);
	Real second_psi = pose.psi(second_res);

	// I use mean dihedral values in Protein Science (1994), 3:2207-2216 "A revised set of potentials for beta-turn formation in proteins" by Hutchinson and Thornton
	// worth to be referred http://en.wikipedia.org/wiki/Turn_(biochemistry)#Hairpins
	// I didn't use Brian's BetaTurnDetectionFeatures since I don't understand it fully
	string turn_type = "uncertain";

	if ( (first_phi > (-64-allowed_deviation_for_turn_type_id_) && first_phi < (-64+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (-27-allowed_deviation_for_turn_type_id_) && first_psi < (-27+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (-90-allowed_deviation_for_turn_type_id_) && second_phi < (-90+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (-7-allowed_deviation_for_turn_type_id_) && second_psi < (-7+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "I";
	} else if ( (first_phi > (-60-allowed_deviation_for_turn_type_id_) && first_phi < (-60+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (131-allowed_deviation_for_turn_type_id_) && first_psi < (131+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (84-allowed_deviation_for_turn_type_id_) && second_phi < (84+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (1-allowed_deviation_for_turn_type_id_) && second_psi < (1+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "II";
	} else if ( (first_phi > (-72-allowed_deviation_for_turn_type_id_) && first_phi < (-72+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (-33-allowed_deviation_for_turn_type_id_) && first_psi < (-33+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (-123-allowed_deviation_for_turn_type_id_) && second_phi < (-123+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (121-allowed_deviation_for_turn_type_id_) && second_psi < (121+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "VIII";
	} else if ( (first_phi > (55-allowed_deviation_for_turn_type_id_) && first_phi < (55+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (38-allowed_deviation_for_turn_type_id_) && first_psi < (38+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (78-allowed_deviation_for_turn_type_id_) && second_phi < (78+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (6-allowed_deviation_for_turn_type_id_) && second_psi < (6+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "I_prime";
	} else if ( (first_phi > (60-allowed_deviation_for_turn_type_id_) && first_phi < (60+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (-126-allowed_deviation_for_turn_type_id_) && first_psi < (-126+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (-91-allowed_deviation_for_turn_type_id_) && second_phi < (-91+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (1-allowed_deviation_for_turn_type_id_) && second_psi < (1+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "II_prime";
	} else if ( (first_phi > (-64-allowed_deviation_for_turn_type_id_) && first_phi < (-64+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (142-allowed_deviation_for_turn_type_id_) && first_psi < (142+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (-93-allowed_deviation_for_turn_type_id_) && second_phi < (-93+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (5-allowed_deviation_for_turn_type_id_) && second_psi < (5+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "VIa1";
	} else if ( (first_phi > (-132-allowed_deviation_for_turn_type_id_) && first_phi < (-132+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (139-allowed_deviation_for_turn_type_id_) && first_psi < (139+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (-80-allowed_deviation_for_turn_type_id_) && second_phi < (-80+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (-10-allowed_deviation_for_turn_type_id_) && second_psi < (-10+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "VIa2";
	} else if ( (first_phi > (-135-allowed_deviation_for_turn_type_id_) && first_phi < (-135+allowed_deviation_for_turn_type_id_))
			&& (first_psi > (131-allowed_deviation_for_turn_type_id_) && first_psi < (131+allowed_deviation_for_turn_type_id_))
			&& (second_phi > (-76-allowed_deviation_for_turn_type_id_) && second_phi < (-76+allowed_deviation_for_turn_type_id_))
			&& (second_psi > (157-allowed_deviation_for_turn_type_id_) && second_psi < (157+allowed_deviation_for_turn_type_id_)) ) {
		turn_type = "VIa2";
	} else {
		turn_type = "IV";
	}


	string select_string =
		"UPDATE sandwich set turn_type = ?\t"
		"WHERE\n"
		"\t(sw_can_by_sh_id = ?) \n"
		"\tAND\t(residue_begin = ?) \n"
		"\tAND (struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1, turn_type);
	select_statement.bind(2, sw_can_by_sh_id);
	select_statement.bind(3, first_res);
	select_statement.bind(4, struct_id);
	basic::database::safely_write_to_database(select_statement);

	return turn_type;
} // WriteToDB_turn_type


Size
WriteToDB_whether_sw_is_not_connected_with_continuous_atoms (
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	string sw_is_not_connected_with_continuous_atoms)
{
	string update =
		"UPDATE sandwich set multimer_is_suspected = ?\t"
		"WHERE\n"
		"\tstruct_id = ? \n"
		"\tAND (sw_can_by_sh_id = ?);";

	statement update_statement(basic::database::safely_prepare_statement(update, db_session));

	update_statement.bind(1,sw_is_not_connected_with_continuous_atoms);
	update_statement.bind(2,struct_id);
	update_statement.bind(3,sw_can_by_sh_id);

	basic::database::safely_write_to_database(update_statement);

	return 0;
} // WriteToDB_whether_sw_is_not_connected_with_continuous_atoms


} //namespace strand_assembly
} //namespace features
} //namespace protocols
