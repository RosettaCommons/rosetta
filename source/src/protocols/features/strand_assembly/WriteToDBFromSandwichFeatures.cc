// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.cc
/// @brief Write to a file after SandwichFeatures
/// @author Doo Nam Kim
/// @overview

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

static basic::Tracer TR("protocols.features.strand_assembly.WriteToDBFromSandwichFeatures");




Size
WriteToDB_avg_b_factor_CB_at_each_component	(
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


	pose::PDB_InfoCOP info = pose.pdb_info();
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

			Real	rounded_avg_b_factor = utility::round(avg_b_factor_CB_at_each_component);
			insert_stmt.bind(1,	rounded_avg_b_factor);
			insert_stmt.bind(2,	struct_id);
			insert_stmt.bind(3,	sw_can_by_sh_id);
			insert_stmt.bind(4,	vector_of_residue_begin[i]);

			basic::database::safely_write_to_database(insert_stmt);
			// <end> UPDATE sw_by_components table

		}
	}
	return 0;
} //	WriteToDB_avg_b_factor_CB_at_each_component




Size
WriteToDB_hydrophobic_ratio_net_charge	(
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
	"	sum(R+K), \n"
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

	int number_of_hydrophobic_res,	number_of_hydrophilic_res,	number_of_CGP,	number_of_RK_in_sw,	number_of_DE_in_sw;
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
	int net_charge_int = number_of_RK_in_sw - number_of_DE_in_sw; // Size net_charge may return like '18446744073709551612', so I don't use Size here
	//	TR << "net_charge_int: " << net_charge_int << endl;
	insert_stmt.bind(7,	net_charge_int); // Net charge of His at pH 7.4 is just '+0.11' according to http://www.bmolchem.wisc.edu/courses/spring503/503-sec1/503-1a.htm
	insert_stmt.bind(8,	sw_can_by_sh_id);
	insert_stmt.bind(9,	struct_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sw_by_components table

	return 0;

} //	WriteToDB_hydrophobic_ratio_net_charge





Size
WriteToDB_min_dis_between_sheets_by_all_res	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Pose & dssp_pose,
	utility::vector1<Size>	all_distinct_sheet_ids)
{

	//// <begin> calculate minimum distance between sheets
	float	minimum_distance_between_sheets_by_all_res	=	cal_min_dis_between_sheets_by_all_res(struct_id,	db_session,	dssp_pose,	all_distinct_sheet_ids);
		// [caution] these values are only meaningful when pdb files has only	one beta-sandwich
	//// <end> calculate minimum distance between sheets

	// <begin> UPDATE sw_by_components table
	string insert =
	"UPDATE sw_by_components set \n"
	"	min_dis_between_sheets_by_all_res = ? \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND	(sw_can_by_sh_id = ?) ;";

	statement insert_stmt(basic::database::safely_prepare_statement(insert,	db_session));
	insert_stmt.bind(1,	round_to_Real(minimum_distance_between_sheets_by_all_res));
	//insert_stmt.bind(1,	round_to_float(minimum_distance_between_sheets_by_all_res));
	insert_stmt.bind(2,	struct_id);
	insert_stmt.bind(3,	sw_can_by_sh_id);

	basic::database::safely_write_to_database(insert_stmt);
	// <end> UPDATE sw_by_components table

	return 0;
} //	WriteToDB_min_dis_between_sheets_by_all_res






void
WriteToDB_turn_AA(
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
		//	TR << "pose.residue_type(i).name3(): " << pose.residue_type(i).name3() << endl;
		if (pose.residue_type(i).name3() == "PHE"	||	pose.residue_type(i).name3() == "VAL"	||	pose.residue_type(i).name3() == "LEU"	||
			pose.residue_type(i).name3() == "ALA"	||	pose.residue_type(i).name3() == "GLY"	||
			pose.residue_type(i).name3() == "TYR"	||	pose.residue_type(i).name3() == "THR"	||
			pose.residue_type(i).name3() == "SER"	||	pose.residue_type(i).name3() == "HIS"	||	pose.residue_type(i).name3() == "GLU"	||
			pose.residue_type(i).name3() == "ASN"	||
			pose.residue_type(i).name3() == "GLN"	||	pose.residue_type(i).name3() == "ASP"	||	pose.residue_type(i).name3() == "ARG")
		{
			//	TR << "pose.residue_type(i+1).name3(): " << pose.residue_type(i+1).name3() << endl;
			if (pose.residue_type(i+1).name3() == "GLY")
			{
				//	TR << "pose.residue_type(i+2).name3(): " << pose.residue_type(i+2).name3() << endl;
				if (pose.residue_type(i+2).name3() == "LEU"	||	pose.residue_type(i+2).name3() == "ALA"	||	pose.residue_type(i+2).name3() == "GLY"	||
					pose.residue_type(i+2).name3() == "PRO"	||
					pose.residue_type(i+2).name3() == "SER"	||	pose.residue_type(i+2).name3() == "GLU"	||
					pose.residue_type(i+2).name3() == "ASN"	||	pose.residue_type(i+2).name3() == "ASP"	|| pose.residue_type(i+2).name3() == "LYS")
				{
					//	TR << "pose.residue_type(i+3).name3(): " << pose.residue_type(i+3).name3() << endl;
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




} //namespace strand_assembly
} //namespace features
} //namespace protocols
