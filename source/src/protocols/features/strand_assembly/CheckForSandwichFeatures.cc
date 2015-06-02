// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/CheckForSandwichFeatures.cc
/// @brief Check various properties for SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)

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

static thread_local basic::Tracer TR( "protocols.features.strand_assembly.CheckForSandwichFeatures" );

//absolute_vec
Real
absolute_vec (numeric::xyzVector<Real> vector)
{
	Real absolute_vec = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	return absolute_vec;
}	//absolute_vec


// calculate_dihedral_w_4_resnums
Real
calculate_dihedral_w_4_resnums(
	Pose const & pose,
	Size res1_sheet_i,
	Size res2_sheet_i,
	Size res1_sheet_j,
	Size res2_sheet_j)
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
} // calculate_dihedral_w_4_resnums


//	cal_dis_angle_to_find_sheet
//	(assign a new strand into a sheet)
vector<Real>
cal_dis_angle_to_find_sheet(
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


//cal_min_avg_dis_between_sheets_by_cen_res
pair<Real, Real>
cal_min_avg_dis_between_sheets_by_cen_res (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	utility::vector1<Size>	all_distinct_sheet_ids,
	Size	min_num_strands_in_sheet_)
{
	Real	min_dis_between_sheets_by_cen_res	=	0;
	Real	avg_dis_between_sheets_by_cen_res	=	0;
	int	appropriate_sheet_num	=	0;

	for(Size i=1; (i <= all_distinct_sheet_ids.size())	&&	(appropriate_sheet_num	!=	2); i++)
	{ // now I check all possible combinations
		if (all_distinct_sheet_ids[i] == 99999) { //all_strands[i].get_size() < min_res_in_strand_
			continue;
		}

		Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

		if (num_strands_i < min_num_strands_in_sheet_)
		{
			continue;
		}
		appropriate_sheet_num++;
		for(Size j=i+1; (j <= all_distinct_sheet_ids.size())	&&	(appropriate_sheet_num	!=	2); j++)
		{
			if (all_distinct_sheet_ids[j] == 99999) { //all_strands[j].get_size() < min_res_in_strand_
				continue;
			}
			Size num_strands_j = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]); // struct_id, db_session, sheet_id

			if (num_strands_j < min_num_strands_in_sheet_)
			{
				continue;
			}
			appropriate_sheet_num++;

			if	(appropriate_sheet_num	==	2)
			{
				pair<Real, Real>	min_dis_AND_avg_dis_between_sheets_by_cen_res	=	cal_min_avg_dis_between_two_sheets_by_cen_res	(struct_id,	db_session,	dssp_pose,	all_distinct_sheet_ids[i],	all_distinct_sheet_ids[j]);
				min_dis_between_sheets_by_cen_res	=	min_dis_AND_avg_dis_between_sheets_by_cen_res.first;
				avg_dis_between_sheets_by_cen_res	=	min_dis_AND_avg_dis_between_sheets_by_cen_res.second;
					//	avg_dis_between_sheets_by_cen_res	is the average between	min_dis_between_sheets_by_cen_res	and 2nd	min_dis_between_sheets_by_cen_res
			}
		}
	}

	if	(appropriate_sheet_num	!=	2)
		// if there are not two sheets,	both avg_dis_between_sheets_by_cen_res and	min_dis_between_sheets_by_cen_res are not meaningful, so don't calculate these
	{
		return std::make_pair(9999.0, 9999.0);
	}

	return std::make_pair(min_dis_between_sheets_by_cen_res, avg_dis_between_sheets_by_cen_res);
} //	cal_min_avg_dis_between_sheets_by_cen_res

//	cal_min_avg_dis_between_two_sheets_by_cen_res
pair<float, float>
cal_min_avg_dis_between_two_sheets_by_cen_res (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sheet_id_1,
	Size sheet_id_2)
{
	vector<Size>	vector_of_cen_residues_in_sheet_1;
	vector_of_cen_residues_in_sheet_1.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_cen_residues_in_sheet_1	=	get_central_residues_in_this_sheet(struct_id, db_session,	sheet_id_1);

	vector<Size>	vector_of_cen_residues_in_sheet_2;
	vector_of_cen_residues_in_sheet_2.clear();	// Removes all elements from the vector (which are destroyed)
	vector_of_cen_residues_in_sheet_2	=	get_central_residues_in_this_sheet(struct_id, db_session,	sheet_id_2);

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
} //	cal_min_avg_dis_between_two_sheets_by_cen_res


//	cal_min_dis_between_sheets_by_all_res
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
} //	cal_min_dis_between_sheets_by_all_res


//	cal_min_dis_between_two_sheets_by_all_res
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
} //	cal_min_dis_between_two_sheets_by_all_res


//	cal_num_of_sheets_that_surround_this_sheet
Size
cal_num_of_sheets_that_surround_this_sheet (
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	utility::vector1<Size>	all_distinct_sheet_ids,
	Size sheet_id,
	Size	min_num_strands_in_sheet_,
	Real	inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_)
{
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

		pair<Real, Real>	min_dis_AND_avg_dis_between_sheets_by_cen_res	=	cal_min_avg_dis_between_two_sheets_by_cen_res	(
			struct_id,
			db_session,
			dssp_pose,
			sheet_id,
			all_distinct_sheet_ids[i]);

		Real	min_dis_between_sheets_by_cen_res	=	min_dis_AND_avg_dis_between_sheets_by_cen_res.first;


		if (min_dis_between_sheets_by_cen_res < inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_){
			num_of_sheets_that_surround_sheet_id++;
		}
	}
		TR << "num_of_sheets_that_surround sheet_id (" << sheet_id << ") within " << inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_ << " Angstrom is " << num_of_sheets_that_surround_sheet_id << endl;
	return num_of_sheets_that_surround_sheet_id;
		// if (num_of_sheets_that_surround_sheet_id > 1) // this sheet is surrounded by more than 1 other sheets!
		// if (num_of_sheets_that_surround_sheet_id == 1) // this sheet is NOT surrounded by more than 1 other sheets, so we can use these sheets to extract sandwich
} //	cal_num_of_sheets_that_surround_this_sheet


//	cal_shortest_dis_between_facing_aro_in_sw
float
cal_shortest_dis_between_facing_aro_in_sw (
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
	utility::vector1<Size>	all_distinct_sheet_ids,
	Size	min_num_strands_in_sheet_)
{
	int	appropriate_sheet_num	=	0;
	float min_distance_between_aro = 9999.0;

	for(Size i=1; (i <= all_distinct_sheet_ids.size())	&&	(appropriate_sheet_num	!=	2); i++)
	{ // now I check all possible combinations
		if (all_distinct_sheet_ids[i] == 99999) { //all_strands[i].get_size() < min_res_in_strand_
			continue;
		}

		Size num_strands_i = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[i]); // struct_id, db_session, sheet_id

		if (num_strands_i < min_num_strands_in_sheet_)
		{
			continue;
		}
		appropriate_sheet_num++;
		for(Size j=i+1; (j <= all_distinct_sheet_ids.size())	&&	(appropriate_sheet_num	!=	2); j++)
		{
			if (all_distinct_sheet_ids[j] == 99999) { //all_strands[j].get_size() < min_res_in_strand_
				continue;
			}
			Size num_strands_j = get_num_strands_in_this_sheet(struct_id, db_session, all_distinct_sheet_ids[j]); // struct_id, db_session, sheet_id

			if (num_strands_j < min_num_strands_in_sheet_)
			{
				continue;
			}
			appropriate_sheet_num++;

			if	(appropriate_sheet_num	==	2)
			{
				vector<Size>	vector_of_aro_residues_in_sheet_1;
				vector_of_aro_residues_in_sheet_1.clear();	// Removes all elements from the vector (which are destroyed)
				vector_of_aro_residues_in_sheet_1	=	get_aro_residues_in_this_sheet(struct_id, db_session,	pose,	all_distinct_sheet_ids[i]);

				vector<Size>	vector_of_aro_residues_in_sheet_2;
				vector_of_aro_residues_in_sheet_2.clear();	// Removes all elements from the vector (which are destroyed)
				vector_of_aro_residues_in_sheet_2	=	get_aro_residues_in_this_sheet(struct_id, db_session,	pose,	all_distinct_sheet_ids[j]);

				for (Size ii=0;	ii<vector_of_aro_residues_in_sheet_1.size();	ii++)
				{
					for (Size	jj=0;	jj<vector_of_aro_residues_in_sheet_2.size();	jj++)
					{
						numeric::xyzVector< core::Real > xyz_of_centroid_of_aro_sheet_1;
						if	((pose.residue_type(vector_of_aro_residues_in_sheet_1[ii]).name3() == "PHE")	||	(pose.residue_type(vector_of_aro_residues_in_sheet_1[ii]).name3() == "TYR"))
						{
							xyz_of_centroid_of_aro_sheet_1.x()	=
							(pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CG ").xyz().x()
							+	pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CZ ").xyz().x())/2;

							xyz_of_centroid_of_aro_sheet_1.y()	=
							(pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CG ").xyz().y()
							+	pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CZ ").xyz().y())/2;

							xyz_of_centroid_of_aro_sheet_1.z()	=
							(pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CG ").xyz().z()
							+	pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CZ ").xyz().z())/2;
						}
						else	//(pose.residue_type(vector_of_aro_residues_in_sheet_1[ii]).name3() == "TRP")
						{
							xyz_of_centroid_of_aro_sheet_1.x()	=
								pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CD2").xyz().x();

							xyz_of_centroid_of_aro_sheet_1.y()	=
								pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CD2").xyz().y();

							xyz_of_centroid_of_aro_sheet_1.z()	=
								pose.residue(vector_of_aro_residues_in_sheet_1[ii]).atom(" CD2").xyz().z();
						}

						numeric::xyzVector< core::Real > xyz_of_centroid_of_aro_sheet_2;
						if	((pose.residue_type(vector_of_aro_residues_in_sheet_2[jj]).name3() == "PHE")	||	(pose.residue_type(vector_of_aro_residues_in_sheet_2[jj]).name3() == "TYR"))
						{
							xyz_of_centroid_of_aro_sheet_2.x()	=
							(pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CG ").xyz().x()
							+	pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CZ ").xyz().x())/2;

							xyz_of_centroid_of_aro_sheet_2.y()	=
							(pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CG ").xyz().y()
							+	pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CZ ").xyz().y())/2;

							xyz_of_centroid_of_aro_sheet_2.z()	=
							(pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CG ").xyz().z()
							+	pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CZ ").xyz().z())/2;
						}
						else	//(pose.residue_type(vector_of_aro_residues_in_sheet_2[ii]).name3() == "TRP")
						{
							xyz_of_centroid_of_aro_sheet_2.x()	=
								pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CD2").xyz().x();

							xyz_of_centroid_of_aro_sheet_2.y()	=
								pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CD2").xyz().y();

							xyz_of_centroid_of_aro_sheet_2.z()	=
								pose.residue(vector_of_aro_residues_in_sheet_2[jj]).atom(" CD2").xyz().z();
						}

						Real distance_between_aro = xyz_of_centroid_of_aro_sheet_1.distance(xyz_of_centroid_of_aro_sheet_2);

						if (distance_between_aro < min_distance_between_aro)
						{
							min_distance_between_aro = distance_between_aro;
						}
					}
				}
			}
		}
	}

	if	(appropriate_sheet_num	!=	2)
		// if there are not two sheets,	cal_shortest_dis_between_facing_aro_in_sw is not meaningful, so don't calculate here
	{
		return 999.0;
	}

	return min_distance_between_aro;
} //	cal_shortest_dis_between_facing_aro_in_sw

//check_canonicalness_of_LR
string
check_canonicalness_of_LR(
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

	if (loop_size == 5)
	{
		if (intra_sheet)
		{
			if (LR=="L" || LR=="BL" ) {return "F_LR";}
			else	{return "T_LR";}
		}
		else {return "U_LR";}
	}

	else // all else loop sizes
	{
		return "U_LR";
	}
} //check_canonicalness_of_LR

//check_canonicalness_of_PA
string
check_canonicalness_of_PA(
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

//check_canonicalness_of_parallel_EE
string
check_canonicalness_of_parallel_EE(
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


//check_heading_direction
string
check_heading_direction( // exclusively between preceding E and following E
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


bool
check_helix_existence(
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


// check_LR
// I tried to slightly vary the method of Nobu's recent nature paper
string
check_LR ( // check L/R chirality of sidechain
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

//check_PA
pair<string, string>
check_PA( // parallel & anti-parallel
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

//check_strand_too_closeness
// <function> check whether these sheets are too close, the closeness is checked for every possible distances
bool
check_strand_too_closeness(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	Real	min_inter_sheet_dis_CA_CA_)
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
} //check_strand_too_closeness


//check_sw_by_dis
Real
check_sw_by_dis(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell, // if false, find parallel way
	Real	min_sheet_dis_,
	Real	max_sheet_dis_
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

			Real avg_dis_CA_CA = get_avg_dis_CA_CA(pose, i_resnum,	i_resnum_1, i_resnum_2, i_resnum_3, j_resnum, j_resnum_1, j_resnum_2, j_resnum_3,	min_sheet_dis_,	max_sheet_dis_);

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
} //check_sw_by_dis


//check_whether_hairpin_connects_short_strand
bool
check_whether_hairpin_connects_short_strand(
	StructureID struct_id,
	sessionOP db_session,
	Size start_res,
	Size next_start_res)
{
	string	select_string =
	"SELECT\n"
	"	component_size	\n"
	"FROM\n"
	"	sandwich \n"
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
	"	sandwich \n"
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


//check_whether_same_direction_strands_connect_two_sheets_or_a_loop
bool
check_whether_same_direction_strands_connect_two_sheets_or_a_loop(
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
			TR.Info << "abs(torsion_between_strands): " << std::abs(torsion_between_strands) << endl;
		if ( std::abs(torsion_between_strands) >	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_)
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


//	check_whether_sheets_are_connected_with_near_bb_atoms
bool
check_whether_sheets_are_connected_with_near_bb_atoms(
	StructureID struct_id,
	sessionOP db_session,
	Pose & dssp_pose,
	Size sw_can_by_sh_id,
	Real	min_N_O_dis_between_two_sheets_,
	Real	min_N_H_O_angle_between_two_sheets_)
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
} //	check_whether_sheets_are_connected_with_near_bb_atoms


// check_whether_sw_by_sh_id_still_alive
bool
check_whether_sw_by_sh_id_still_alive(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	sandwich_PK_id \n"
	"FROM\n"
	"	sandwich\n"
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


//check_whether_sw_is_not_connected_with_continuous_atoms
string
check_whether_sw_is_not_connected_with_continuous_atoms(
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
	"	sandwich \n"
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


//check_whether_this_pdb_should_be_excluded
bool
check_whether_this_pdb_should_be_excluded (
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
} //check_whether_this_pdb_should_be_excluded


//check whether this sheet is constituted with 2 (or less) residues long only
bool
check_whether_this_sheet_is_too_short(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_i)
{
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


//check_whether_strand_i_is_in_sheet
//	<role> Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
bool
check_whether_strand_i_is_in_sheet(
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


//count_AA_w_direction
vector<Size>
count_AA_w_direction(
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


//count_AA_wo_direction
vector<Size>
count_AA_wo_direction(
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
} //count_AA_wo_direction


//determine_core_heading_surface_heading_by_distance
string
determine_core_heading_surface_heading_by_distance(
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
}	//determine_core_heading_surface_heading_by_distance


//determine_heading_direction_by_vector
string
determine_heading_direction_by_vector
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
		vector_of_cen_residues	=	get_central_residues_in_other_sheet(struct_id, db_session, sw_can_by_sh_id,	sheet_id);

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

//get_current_strands_in_sheet
//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
get_current_strands_in_sheet(
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

//get_distinct_sheet_id_from_sheet_table
utility::vector1<Size>
get_distinct_sheet_id_from_sheet_table(
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

//get_full_strands
/// @brief Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<SandwichFragment>
get_full_strands(
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
}//get_full_strands


//find_sheet (assign a new strand into a sheet)
Size
find_sheet(
	Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell, // if 'false', find a sheet in parallel way
	Real	min_CA_CA_dis_,
	Real	max_CA_CA_dis_,
	Real	min_C_O_N_angle_
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
} //find_sheet


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


//get_all_strands_in_sheet_i
//Select all strand segments reported by the secondary_structure_segments and save them in a vector
utility::vector1<SandwichFragment>
get_all_strands_in_sheet_i(
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


//get_aro_residues_in_this_sheet
vector<Size>
get_aro_residues_in_this_sheet(
	StructureID struct_id,
	sessionOP db_session,
	Pose const & pose,
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

	vector<Size> vector_of_aro_residues;
	for(Size i=0; i<vector_of_residue_begin.size(); i++ )
	{
		for(Size res	=	vector_of_residue_begin[i]; res	<=	vector_of_residue_end[i]; res++ )
		{
			if	((pose.residue_type(res).name3() == "PHE") || (pose.residue_type(res).name3() == "TRP") || (pose.residue_type(res).name3() == "TYR"))
			{
					vector_of_aro_residues.push_back(res);
			}
		}
	}
	return vector_of_aro_residues;
} //get_aro_residues_in_this_sheet


//get_avg_dis_CA_CA
Real
get_avg_dis_CA_CA(
	Pose const & pose,
	Size i_resnum,
	Size i_resnum_1,
	Size i_resnum_2,
	Size i_resnum_3,
	Size j_resnum,
	Size j_resnum_1,
	Size j_resnum_2,
	Size j_resnum_3,
	Real	min_sheet_dis_,
	Real	max_sheet_dis_)
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
} // get_avg_dis_CA_CA


//get_avg_dis_strands
// check whether these sheets are too close, the closeness is checked for every possible distances
// <where to use	now?> get_current_bs_id_and_closest_edge_bs_id_in_different_sheet
// <where to use	now?> in main, check whether strands are too distant to each other
// <where to use	now?> see_whether_sheet_is_antiparallel
Real
get_avg_dis_strands(
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
} // get_avg_dis_strands


//get_central_residues_in_other_sheet
vector<Size>
get_central_residues_in_other_sheet(
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
} //get_central_residues_in_other_sheet


//get_central_residues_in_this_sheet
vector<Size>
get_central_residues_in_this_sheet(
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
} //get_central_residues_in_this_sheet


//get_chain_B_resNum
utility::vector1<Size>
get_chain_B_resNum(
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
	"	sandwich AS sw, \n"
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


//get_central_residues_in_each_of_two_edge_strands
pair<int, int>
get_central_residues_in_each_of_two_edge_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i,
	Real	min_CA_CA_dis_,
	Real	max_CA_CA_dis_)
{
	utility::vector1<SandwichFragment> strands_from_sheet_i = get_full_strands_from_sheet(struct_id, db_session, sheet_i);

	// get central residue numbers in edge strands
	vector<int> vector_of_central_residues_in_sheet_i;
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
										strands_from_sheet_i[i].get_end(),
										min_CA_CA_dis_,
										max_CA_CA_dis_);
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

	int terminal_cen_res_pos_1 = vector_of_central_residues_in_sheet_i[index_terminal_cen_res_pos_1];	// index of the first central_residue
	int terminal_cen_res_pos_2 = vector_of_central_residues_in_sheet_i[index_terminal_cen_res_pos_2];	// index of the second central_residue

	return std::make_pair(terminal_cen_res_pos_1, terminal_cen_res_pos_2);
} //get_central_residues_in_each_of_two_edge_strands


// get_current_bs_id_and_closest_edge_bs_id_in_different_sheet
pair<Size, Size>
get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end)
{
	// <begin> retrieve current sandwich_bs_id
	string select_string =
	"SELECT\n"
	"	sandwich_bs_id \n"
	"FROM\n"
	"	sandwich \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?) \n"
	"	AND (residue_begin = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	sw_can_by_sh_id);
	select_statement.bind(3,	residue_begin);
	result res(basic::database::safely_read_from_database(select_statement));

	Size current_sandwich_bs_id;
	while(res.next())
	{
		res >> current_sandwich_bs_id;
	}
	// <end> retrieve current sandwich_bs_id


	// <begin> retrieve other edge_strands in different sheet
	string select_string_2 =
	"SELECT\n"
	"	residue_begin, \n"
	"	residue_end \n"
	"FROM\n"
	"	sandwich \n"
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
		// <begin> retrieve the closest sandwich_bs_id
		string select_string_3 =
		"SELECT\n"
		"	sandwich_bs_id \n"
		"FROM\n"
		"	sandwich \n"
		"WHERE\n"
		"	(struct_id = ?) \n"
		"	AND (sw_can_by_sh_id = ?) \n"
		"	AND (residue_begin = ?);";

		statement select_statement_3(basic::database::safely_prepare_statement(select_string_3,db_session));
		select_statement_3.bind(1,	struct_id);
		select_statement_3.bind(2,	sw_can_by_sh_id);
		select_statement_3.bind(3,	residue_begin_of_nearest_strand);
		result res_3(basic::database::safely_read_from_database(select_statement_3));

		Size closest_sandwich_bs_id;
		while(res_3.next())
		{
			res_3 >> closest_sandwich_bs_id;
		}
		// <end> retrieve the closest sandwich_bs_id

	// <end> see which other edge_strand in different sheet is closest to a current strand

	return std::make_pair(current_sandwich_bs_id, closest_sandwich_bs_id);

} // get_current_bs_id_and_closest_edge_bs_id_in_different_sheet


//get_distinct_sw_id_from_sandwich_table
utility::vector1<Size>
get_distinct_sw_id_from_sandwich_table(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	distinct sw_can_by_sh_id\n"
	"FROM\n"
	"	sandwich \n"
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
} //get_distinct_sw_id_from_sandwich_table


//get_full_strands_from_sheet
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


//get_next_starting_res_for_connecting_strands
pair<Size, Size> //Size
get_next_starting_res_for_connecting_strands(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Size former_ending_res)
{
	string select_string =
	"SELECT\n"
	"	min(residue_begin) \n"
	"FROM\n"
	"	sandwich \n"
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
	"	sandwich \n"
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

//get_num_of_distinct_sheet_id
Size
get_num_of_distinct_sheet_id(
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


//get_num_of_sheets_that_surround_this_sheet
Size
get_num_of_sheets_that_surround_this_sheet(
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


//	get_list_of_residues_in_sheet_i
utility::vector1<Size>
get_list_of_residues_in_sheet_i(
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

//get_max_sheet_id
Size
get_max_sheet_id(
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


//get_segment_id
Size
get_segment_id(
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


// get_shortest_among_4_vals (simple one with just four parameters)
Real
get_shortest_among_4_vals(
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
} // get_shortest_among_4_vals (simple one with just four parameters)


//get_size_sandwich_PK_id
Size
get_size_sandwich_PK_id(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	count(sandwich_PK_id) \n"
	"FROM\n"
	"	sandwich \n"
	"WHERE\n"
	"	(struct_id = ?) \n"
	"	AND (sw_can_by_sh_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	select_statement.bind(2,sw_can_by_sh_id);
	result res(basic::database::safely_read_from_database(select_statement));

	Size size_of_sandwich_PK_id;
	while(res.next())
	{
		res >> size_of_sandwich_PK_id;
	}
	return size_of_sandwich_PK_id;
} //get_size_sandwich_PK_id


// get_start_end_res_num_in_the_longest_strand which is used for judge_facing
utility::vector1<SandwichFragment>
get_start_end_res_num_in_the_longest_strand(
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


//get_starting_res_for_connecting_strands
pair<Size, Size>
get_starting_res_for_connecting_strands(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id,
	Size former_res_end)
{
	string select_string =
	"SELECT\n"
	"	min(residue_end) \n"
	"FROM\n"
	"	sandwich \n"
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
	"	sandwich \n"
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
} //	get_starting_res_for_connecting_strands

//get_tag
string
get_tag(
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

utility::vector1<Size>
get_vec_AA_kind (
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
	"	sandwich\n"
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


// get_vec_of_sw_can_by_sh_id
// role:	get_distinct(sw_can_by_sh_id) as a vector
// result:	mostly get_distinct(sw_can_by_sh_id) is just 1
utility::vector1<Size>
get_vec_of_sw_can_by_sh_id(
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
} //get_vec_of_sw_can_by_sh_id


utility::vector1<Size>
get_vector_of_strand_AA_distribution (
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
		"	sandwich \n"
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
		"	sandwich \n"
		"WHERE\n"
		"	strand_edge = ? \n"
		"	AND struct_id = ? ;";
	}
	statement sum_statement(basic::database::safely_prepare_statement(sum_string, db_session));
	sum_statement.bind(1,	strand_location);
	sum_statement.bind(2,	struct_id);

	result res(basic::database::safely_read_from_database(sum_statement));

	Size num_A, num_C,	num_D,	num_E, num_F, num_G, num_H, num_I, num_K, num_L, num_M, num_N, num_P, num_Q, num_R, num_S, num_T, num_V, num_W, num_Y;

	utility::vector1<Size> vector_strand_AA_distribution;

	while(res.next())
	{
		res >> num_A >>  num_C >> 	num_D >> 	num_E >>  num_F >>  num_G >>  num_H >>  num_I >>  num_K >>  num_L >>  num_M >>  num_N >>  num_P >>  num_Q >>  num_R >>  num_S >>  num_T >>  num_V >>  num_W >>  num_Y;
		vector_strand_AA_distribution.push_back(num_A);
		vector_strand_AA_distribution.push_back(num_C);
		vector_strand_AA_distribution.push_back(num_D);
		vector_strand_AA_distribution.push_back(num_E);
		vector_strand_AA_distribution.push_back(num_F);

		vector_strand_AA_distribution.push_back(num_G);
		vector_strand_AA_distribution.push_back(num_H);
		vector_strand_AA_distribution.push_back(num_I);
		vector_strand_AA_distribution.push_back(num_K);
		vector_strand_AA_distribution.push_back(num_L);

		vector_strand_AA_distribution.push_back(num_M);
		vector_strand_AA_distribution.push_back(num_N);
		vector_strand_AA_distribution.push_back(num_P);
		vector_strand_AA_distribution.push_back(num_Q);
		vector_strand_AA_distribution.push_back(num_R);

		vector_strand_AA_distribution.push_back(num_S);
		vector_strand_AA_distribution.push_back(num_T);
		vector_strand_AA_distribution.push_back(num_V);
		vector_strand_AA_distribution.push_back(num_W);
		vector_strand_AA_distribution.push_back(num_Y);
	}

	return vector_strand_AA_distribution;
} // get_vector_of_strand_AA_distribution

//get_vec_distinct_sheet_id
utility::vector1<Size>
get_vec_distinct_sheet_id(
	StructureID struct_id,
	sessionOP db_session,
	Size sw_can_by_sh_id)
{
	string select_string =
	"SELECT\n"
	"	distinct sheet_id\n"
	"FROM\n"
	"	sandwich \n"
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


// identify_sheet_id_by_residue_end
Size
identify_sheet_id_by_residue_end(
	StructureID struct_id,
	sessionOP db_session,
	Size residue_end)
{
	string select_string =
	"SELECT\n"
	"	sheet_id \n"
	"FROM\n"
	"	sandwich \n"
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


// is_this_strand_at_edge
// <role> This function is used to write edge_strand info into db
// <note> If a strand is within "short_edge", the strand is not an "edge"
string
is_this_strand_at_edge	(
	Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id,
	Size residue_begin,
	Size residue_end,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_)
{
//		TR << "is_this_strand_at_edge" << endl;
//	time_t start_time = time(NULL);

	if (residue_end - residue_begin + 1 < 3)
	{
		return "short_edge"; // Like res_num: 78-79 in [1ten], 2 residues long strand is better to be classified as 'short edge'
	}

	// <begin> see whether this sheet is consisted with two strands only
	//	get_full_strands_from_sheet is refactored
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
		// get_closest_distance_between_strands is refactored
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


//is_this_strand_at_edge_by_looking_db
string
is_this_strand_at_edge_by_looking_db(
	StructureID struct_id,
	sessionOP db_session,
	Size residue_begin)
{
	string select_string =
	"SELECT\n"
	"	strand_edge\n"
	"FROM\n"
	"	sandwich \n"
	"WHERE\n"
	"	struct_id = ? \n"
	"  AND residue_begin = ? ;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	select_statement.bind(2,	residue_begin);
	result res(basic::database::safely_read_from_database(select_statement));

	string edge;
	while(res.next())
	{
		res >> edge ;
	}
	return edge;
} //is_this_strand_at_edge_by_looking_db


//	judge_facing
int
judge_facing(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size sheet_i,
	Size sheet_j,
	Real	min_CA_CA_dis_,
	Real	max_CA_CA_dis_,
	Real	min_sheet_angle_by_four_term_cen_res_,
	Real	max_sheet_angle_by_four_term_cen_res_,
	Real	min_sheet_torsion_cen_res_,
	Real	max_sheet_torsion_cen_res_,
	Real	max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_)
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
		return 0; // these two sheets are linear or do not face each other properly!
	}


	// <begin> identify four terminal central residues
	pair<int, int>
	two_central_residues_in_two_edge_strands =	get_central_residues_in_each_of_two_edge_strands(
		struct_id,
		db_session,
		pose,
		sheet_i,
		min_CA_CA_dis_,
		max_CA_CA_dis_);

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
		sheet_j,
		min_CA_CA_dis_,
		max_CA_CA_dis_);

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
		(torsion_i_j > min_sheet_torsion_cen_res_)	&&
		(torsion_i_j < max_sheet_torsion_cen_res_))
	{
		return 1; // these two strand_pairs face each other properly, so constitute a sandwich
	}

	else
	{
			TR <<	"these two strand_pairs are linear or	do not face	each other properly!" << endl;
		return 0;
	}

} // judge_facing


void
process_decoy(
	Pose &dssp_pose,
	core::scoring::ScoreFunction const& scorefxn )
{
	scorefxn( dssp_pose );
} // process_decoy


// report_heading_directions_of_all_AA_in_a_strand
string
report_heading_directions_of_all_AA_in_a_strand	(
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
		// determine_heading_direction_by_vector is refactored well
		string	heading	=	determine_heading_direction_by_vector	(struct_id,	db_session,	pose,	sw_can_by_sh_id,	sheet_id,	residue_begin,	residue_end,	ii);
		heading_directions	= heading_directions	+	"	"	+	heading;
	}

	return heading_directions;
} // report_heading_directions_of_all_AA_in_a_strand


utility::vector1<int>
retrieve_residue_num_of_rkde(
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

		utility::vector1<int> vector_of_residue_num_of_rkde;
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

			utility::vector1<int> vector_of_residue_num_of_rkde;
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

			utility::vector1<int> vector_of_residue_num_of_rkde;
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


//	(07/19/14)
//	Doonam confirms that 'round_to_float' didn't round as expected with 'float' argument
float
round_to_float(
	float x)
{
	return floor((x	*	10)	+	0.5)	/	10;
} //round_to_float


//	(07/19/14)
//	Doonam confirms that 'round_to_Real' rounds as expected with 'float' argument
Real
round_to_Real(
	Real x)
{
	Real rounded = floor((x * 10) +	0.5)	/	10;
	return rounded;
} //round_to_Real

//round_to_Size
Size
round_to_Size(
	Real x)
{
	Size rounded = static_cast <Size> (floor(x+.5));
	return rounded;
} //round_to_Size

//	see_edge_or_core_or_loop_or_short_edge
string
see_edge_or_core_or_loop_or_short_edge (
	StructureID struct_id,
	sessionOP db_session,
	Size	residue_num)
{
	string	sum_string =
		"SELECT\n"
		"	strand_edge \n"
		"FROM\n"
		"	sandwich \n"
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

	string	edge_or_core;
	if	(strand_edge == "edge")
	{
		edge_or_core = "edge";
	}
	else if	(strand_edge == "core")
	{
		edge_or_core = "core";
	}
	else if (strand_edge == "short_edge")
	{
		edge_or_core = "short_edge";
	}
	else
	{
		edge_or_core = "loop";
	}
	return edge_or_core;
}	//	see_edge_or_core_or_loop_or_short_edge


// see_whether_sheet_is_antiparallel
// <Definition> In order to be an antiparallel sheet, all strands in the sheet should be anti-parallel to each other
// <Overall plot>
	// 1. Get central residues of each edge strands (save -999 as resnum if the strand is at core)
	// 2. Get array of sum of all distances from central residues
	// 3. Get the largest distance and its residue index
	// 4. Search the nearest strand from the current strand
string
see_whether_sheet_is_antiparallel(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size i_sheet,
	Real	min_CA_CA_dis_,
	Real	max_CA_CA_dis_,
	Real	min_C_O_N_angle_)
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
										strands_from_i[i].get_end(),
										min_CA_CA_dis_,
										max_CA_CA_dis_);
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
		Size return_of_find_sheet =
		find_sheet (
			pose,
			current_strand, nearest_strand_from_the_current_strand, true,
			min_CA_CA_dis_,
			max_CA_CA_dis_,
			min_C_O_N_angle_);

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
} //see_whether_sheet_is_antiparallel


//see_whether_sheets_can_be_combined
bool
see_whether_sheets_can_be_combined(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Pose const & pose,
	Size i_sheet,
	Size j_sheet,
	Real	min_CA_CA_dis_,
	Real	max_CA_CA_dis_,
	Real	min_C_O_N_angle_)
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

			return_of_find_sheet_antiparallel =
				find_sheet (
					pose, temp_strand_i, temp_strand_j, true,
					min_CA_CA_dis_,
					max_CA_CA_dis_,
					min_C_O_N_angle_);

			if (return_of_find_sheet_antiparallel == 999)	// since these two strands are too distant to each other, there is virtually no chance to be sheet!
			{
				break; // since these two strands are too distant to each other, there is virtually no chance to be sheet!
			}

			if (!return_of_find_sheet_antiparallel)
			{
				return_of_find_sheet_parallel =
				find_sheet (
					pose, temp_strand_i, temp_strand_j, false,
					min_CA_CA_dis_,
					max_CA_CA_dis_,
					min_C_O_N_angle_);
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
} //see_whether_sheets_can_be_combined

// see_whether_this_sw_has_SS_bond
bool
see_whether_this_sw_has_SS_bond(
	StructureID struct_id,
	sessionOP db_session)
{
	string select_string =
	"SELECT\n"
	"	resNum1 \n"
	"FROM\n"
	"	residue_scores_lr_2b\n"
	"WHERE\n"
	"	(struct_id = ?);";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,	struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	bool sw_has_SS_bond = false;
	while(res.next())
	{
		sw_has_SS_bond = true;
	}
	return sw_has_SS_bond;
} //see_whether_this_sw_has_SS_bond


} //namespace strand_assembly
} //namespace features
} //namespace protocols
