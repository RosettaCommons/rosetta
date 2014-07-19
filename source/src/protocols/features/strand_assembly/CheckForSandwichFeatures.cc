// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/CheckForSandwichFeatures.cc
/// @brief Check various properties for SandwichFeatures
/// @author Doo Nam Kim
/// @overview

//Devel
#include <protocols/features/strand_assembly/CheckForSandwichFeatures.hh>
//#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace core;
using namespace std;
using core::pose::Pose;
using cppdb::statement;
using cppdb::result;
using utility::sql_database::sessionOP;
using utility::vector1;

static basic::Tracer TR("protocols.features.strand_assembly.CheckForSandwichFeatures");


//cal_min_dis_between_sheets_by_all_res
float
cal_min_dis_between_sheets_by_all_res (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	utility::vector1<Size>	all_distinct_sheet_ids)
{
	float	min_dis_between_two_sheets_by_all_res	=	9999.0;
	int	appropriate_sheet_num	=	0;

	for(Size i=1; (i <= all_distinct_sheet_ids.size())	&&	(appropriate_sheet_num	!=	2); i++)
	{ // now I check all possible combinations
		if (all_distinct_sheet_ids[i] == 99999) { //all_strands[i].get_size() < min_res_in_strand_
			continue;
		}

		Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

		//if (num_strands_i < min_num_strands_in_sheet_)
		if (num_strands_i < 3)
		{
			continue;
		}
		appropriate_sheet_num++;
		for(Size j=i+1; (j <= all_distinct_sheet_ids.size())	&&	(appropriate_sheet_num	!=	2); j++)
		{
			if (all_distinct_sheet_ids[j] == 99999)
			{ //all_strands[j].get_size() < min_res_in_strand_
				continue;
			}
			Size num_strands_j = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]); // struct_id, db_session, sheet_id

			//if (num_strands_j < min_num_strands_in_sheet_)
			if (num_strands_j < 3)
			{
				continue;
			}
			appropriate_sheet_num++;

			if	(appropriate_sheet_num	==	2)
			{
				min_dis_between_two_sheets_by_all_res	=	cal_min_dis_between_two_sheets_by_all_res	(struct_id,	db_session,	dssp_pose,	all_distinct_sheet_ids[i],	all_distinct_sheet_ids[j]);
			}
		}
	}

	if	(appropriate_sheet_num	!=	2)
		// if there are not two sheets,	both avg_dis_between_sheets_by_cen_res and	min_dis_between_sheets_by_cen_res are not meaningful, so don't calculate these
	{
		return 9999.0; // default of min_dis_between_two_sheets_by_all_res
	}
	return min_dis_between_two_sheets_by_all_res;
} //cal_min_dis_between_sheets_by_all_res



//cal_min_dis_between_two_sheets_by_all_res
//std::pair<float, float>
float
cal_min_dis_between_two_sheets_by_all_res (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sheet_id_1,
	Size sheet_id_2)
{
	vector<Size>	vector_of_all_residues_in_sheet_1;
	vector_of_all_residues_in_sheet_1.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_all_residues_in_sheet_1	=	get_all_residues_in_this_sheet(struct_id, db_session,	sheet_id_1);

	vector<Size>	vector_of_all_residues_in_sheet_2;
	vector_of_all_residues_in_sheet_2.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_all_residues_in_sheet_2	=	get_all_residues_in_this_sheet(struct_id, db_session,	sheet_id_2);

	float min_dis_between_two_sheets_by_all_res = 9999;

	for (Size ii=0;	ii<vector_of_all_residues_in_sheet_1.size();	ii++)
	{
		for (Size	jj=0;	jj<vector_of_all_residues_in_sheet_2.size();	jj++)
		{
			Real current_distance = dssp_pose.residue(vector_of_all_residues_in_sheet_1[ii]).atom("CA").xyz().distance(dssp_pose.residue(vector_of_all_residues_in_sheet_2[jj]).atom("CA").xyz());
			if (current_distance < min_dis_between_two_sheets_by_all_res)
			{
				min_dis_between_two_sheets_by_all_res = current_distance;
			}
		}
	}
	//	TR << "min_dis between sheet: " << sheet_id_1 << " and sheet: " << sheet_id_2 << " is " << min_dis << endl;

	return min_dis_between_two_sheets_by_all_res;
} //cal_min_dis_between_two_sheets_by_all_res


//cal_min_avg_dis_between_two_sheets_by_cen_res
std::pair<float, float>
cal_min_avg_dis_between_two_sheets_by_cen_res (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sheet_id_1,
	Size sheet_id_2)
{
	vector<Size>	vector_of_cen_residues_in_sheet_1;
	vector_of_cen_residues_in_sheet_1.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_cen_residues_in_sheet_1	=	get_cen_residues_in_this_sheet(struct_id, db_session,	sheet_id_1);

	vector<Size>	vector_of_cen_residues_in_sheet_2;
	vector_of_cen_residues_in_sheet_2.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_cen_residues_in_sheet_2	=	get_cen_residues_in_this_sheet(struct_id, db_session,	sheet_id_2);

	float min_dis = 9999;
	float second_min_dis = 9999;
	float sum_dis = 0.0;
	float cal_num	=	0.0;

	for (Size ii=0;	ii<vector_of_cen_residues_in_sheet_1.size();	ii++){
		for (Size	jj=0;	jj<vector_of_cen_residues_in_sheet_2.size();	jj++){
			Real distance = dssp_pose.residue(vector_of_cen_residues_in_sheet_1[ii]).atom("CA").xyz().distance(dssp_pose.residue(vector_of_cen_residues_in_sheet_2[jj]).atom("CA").xyz());
			sum_dis	=	sum_dis	+	distance;
			cal_num	=	cal_num	+	1.0;
			if (distance < min_dis)
			{
				min_dis = distance;
			}
			if ((distance < second_min_dis)	&&	(distance	>	min_dis))
			{
				second_min_dis = distance;
			}
		}
	}
	//	TR << "min_dis between sheet: " << sheet_id_1 << " and sheet: " << sheet_id_2 << " is " << min_dis << endl;

	return std::make_pair(min_dis, ((min_dis	+	second_min_dis)/2));
		// (min_dis	+	second_min_dis)/2) is a better description of an average of min_dis since I deal with distance between central residues only
		// If I calculate average of min_dis with just 'sum_dis/cal_num', it will give unwanted average values, what I want is just average distance between "possible" shortest min_dis only
} //cal_min_avg_dis_between_two_sheets_by_cen_res


// check_whether_sw_by_sh_id_still_alive
bool
check_whether_sw_by_sh_id_still_alive(
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




//get_all_residues_in_this_sheet
vector<Size>
get_all_residues_in_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
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

	vector<Size> vector_of_all_residues;
	for(Size i=0; i<vector_of_residue_begin.size(); i++ )
	{
		for(Size all_resnum_in_this_sheet = vector_of_residue_begin[i]; all_resnum_in_this_sheet <= vector_of_residue_end[i]; all_resnum_in_this_sheet++)
		{
			vector_of_all_residues.push_back(all_resnum_in_this_sheet);
		}
	}
	return vector_of_all_residues;
} //get_all_residues_in_this_sheet


//get_cen_residues_in_this_sheet
vector<Size>
get_cen_residues_in_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Size sheet_id)
{
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
	return vector_of_cen_residues;
} //get_cen_residues_in_this_sheet


// check whether these sheets are too close, the closeness is checked for every possible distances
Real
get_closest_distance_between_strands(
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
} //get_closest_distance_between_strands


/// @brief Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
get_full_strands_from_sheet(
	StructureID struct_id,
	sessionOP db_session,
	core::Size sheet_id)
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
}	//get_full_strands_from_sheet



//get_num_strands_in_this_sheet
Size
get_num_strands_in_this_sheet(
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


//get_sheet_antiparallel_info
string
get_sheet_antiparallel_info(
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

// (07/19/14) but I(Doonam) confirms that 'round_to_Real' rounds as expected with 'float' argument
// 'round_to_float' didn't round as expected with 'float' argument
float
round_to_float(
	float x)
{
	return floor((x	*	10)	+	0.5)	/	10;
} //round_to_float


Real
round_to_Real(
	Real x)
{
	Real rounded = floor((x * 10) +	0.5)	/	10;
	return rounded;
} //round_to_Real


Size
round_to_Size(
	Real x)
{
	Size rounded = static_cast <Size> (floor(x+.5));
	return rounded;
} //round_to_Size



} //namespace strand_assembly
} //namespace features
} //namespace protocols
