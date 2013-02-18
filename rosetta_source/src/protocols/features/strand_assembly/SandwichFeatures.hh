// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file SandwichFeatures.hh
/// @brief
/// @author Doo Nam Kim (started with Tim Jacobs' code)

#ifndef INCLUDED_protocols_features_strand_assembly_SandwichFeatures_hh
#define INCLUDED_protocols_features_strand_assembly_SandwichFeatures_hh

//Unit
#include <protocols/features/strand_assembly/SandwichFeatures.fwd.hh>

//External
#include <boost/uuid/uuid.hpp>

//Protocols
#include <protocols/features/FeaturesReporter.hh>

//Devel
#include <protocols/features/strand_assembly/SandwichFragment.hh>

//Utility
#include <utility/vector1.hh>

// for string return
#include <string>

#include <vector> // for get_sw_can_by_sh_id

namespace protocols {
namespace features {
namespace strand_assembly {

class SandwichFeatures : public protocols::features::FeaturesReporter
{

public:

	SandwichFeatures();

	void init_from_options();

	virtual
	std::string
	type_name() const
	{
		return "SandwichFeatures";
	}

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
		core::pose::Pose const & pose, //core::pose::Pose & pose, // dropped 'const' for dssp info addition
		utility::vector1<bool> const & relevant_residues,
		boost::uuids::uuid struct_id,
		utility::sql_database::sessionOP db_session);

	utility::vector1<SandwichFragment>
	get_full_strands(
					 boost::uuids::uuid struct_id,
					 utility::sql_database::sessionOP db_session);

	bool
	check_whether_strand_i_is_in_sheet(
									   boost::uuids::uuid struct_id,
									   utility::sql_database::sessionOP db_session,
									   core::Size segment_id);

	utility::vector1<SandwichFragment>
	get_current_strands_in_sheet(
								 boost::uuids::uuid struct_id,
								 utility::sql_database::sessionOP db_session);

	core::Size
	get_max_sheet_id(
					 boost::uuids::uuid struct_id,
					 utility::sql_database::sessionOP db_session);

	utility::vector1<SandwichFragment>
	get_chain_B_resNum(
					   boost::uuids::uuid struct_id,
					   utility::sql_database::sessionOP db_session);

	std::string
	get_tag(
			boost::uuids::uuid struct_id,
			utility::sql_database::sessionOP db_session);


	core::Size
	get_num_strands(
					boost::uuids::uuid struct_id,
					utility::sql_database::sessionOP db_session,
					core::Size sheet_id);


	utility::vector1<SandwichFragment>
	get_full_strands_from_sheet(
								boost::uuids::uuid struct_id,
								utility::sql_database::sessionOP db_session,
								core::Size sheet_id);

	//utility::vector1<SandwichFragment>
	utility::vector1<core::Size>
	get_distinct_sheet_id(
						  boost::uuids::uuid struct_id,
						  utility::sql_database::sessionOP db_session);

	core::Size
	find_sheet	(
				 core::pose::Pose const & pose,
				 SandwichFragment strand_i,
				 SandwichFragment strand_j,
				 bool antiparalell // if false, try to find a sheet in parallel way
				 );

	std::vector<core::Real>
	cal_dis_angle_to_find_sheet( // calculate distance and angle to find sheet
								core::pose::Pose const & pose,
								core::Size res_i_0,
								core::Size res_i_1,
								core::Size res_i_2,
								core::Size res_j_0,
								core::Size res_j_1,
								core::Size res_j_2);

	bool
	see_whether_sheets_can_be_combined(
									   boost::uuids::uuid struct_id,
									   utility::sql_database::sessionOP db_session,
									   core::pose::Pose const & pose,
									   core::Size i_sheet,
									   core::Size j_sheet);

	core::Size
	update_sheet_id(
					boost::uuids::uuid struct_id,
					utility::sql_database::sessionOP db_session,
					core::Size new_sheet_id,
					core::Size old_sheet_id);

	core::Size
	update_sheet_antiparallel(
							  boost::uuids::uuid struct_id,
							  utility::sql_database::sessionOP db_session,
							  core::Size sheet_id,
							  core::Size antiparallel);

	core::Size	
	get_sheet_antiparallel_info(
								boost::uuids::uuid struct_id,
								utility::sql_database::sessionOP db_session,
								core::Size sheet_id);
		

	bool
	change_sheet_id_if_possible(
								boost::uuids::uuid struct_id,
								utility::sql_database::sessionOP db_session,
								core::pose::Pose const & pose);

	core::Real
	check_sw_by_dis(
					core::pose::Pose const & pose,
					SandwichFragment strand_i,
					SandwichFragment strand_j,
					bool antiparalell);

	bool
	see_whether_sheet_is_antiparallel(
									  boost::uuids::uuid struct_id,
									  utility::sql_database::sessionOP db_session,
									  core::pose::Pose const & pose,
									  core::Size i_sheet);


	bool
	check_strand_too_closeness	(
								 core::pose::Pose const & pose,
								 SandwichFragment strand_i,
								 SandwichFragment strand_j);

	core::Real
	get_avg_dis_strands(
						core::pose::Pose const & pose,
						SandwichFragment strand_i,
						SandwichFragment strand_j);

	core::Real
	get_avg_dis_CA_CA(
					  core::pose::Pose const & pose,
					  core::Size i_resnum,
					  core::Size i_resnum_1,
					  core::Size i_resnum_2,
					  core::Size i_resnum_3,
					  core::Size j_resnum,
					  core::Size j_resnum_1,
					  core::Size j_resnum_2,
					  core::Size j_resnum_3);

	core::Size round(
					 core::Real x);


	bool
	can_this_strand_represent_a_terminal(
										 core::pose::Pose const & pose,
										 utility::vector1<SandwichFragment> strands_from_sheet_i,
										 core::Size current_strand_id_as_i);

	std::pair<core::Real, core::Real>
	get_two_central_residues(
							 boost::uuids::uuid struct_id,
							 utility::sql_database::sessionOP db_session,
							 core::pose::Pose const & pose,
							 core::Size sheet_i);

	core::Real
	get_shortest_among_4(
						 core::Real arr_dis_inter_sheet[]);

	bool
	judge_facing(
				 boost::uuids::uuid struct_id,
				 utility::sql_database::sessionOP db_session,
				 core::pose::Pose const & pose,
				 core::Size sheet_i,
				 core::Size sheet_j);

	core::Size
	write_to_sheet (
					boost::uuids::uuid struct_id,
					utility::sql_database::sessionOP db_session,
					core::Size sheet_PK_id_counter,
					core::Size sheet_id,
					core::Size segment_id);

	core::Size
	write_to_sw_can_by_sh	(
							 boost::uuids::uuid struct_id,
							 utility::sql_database::sessionOP db_session,
							 core::Size sw_can_by_sh_PK_id_counter,
							 std::string tag,
							 core::Size sw_can_by_sh_id_counter,
							 core::Size sheet_id,
							 core::Size num_strands_from_sheet);

	utility::vector1<SandwichFragment>
	prepare_to_fill_sw_can_by_components(
										 boost::uuids::uuid struct_id,
										 utility::sql_database::sessionOP db_session);


	bool	
	see_whether_strand_is_at_edge	(
									 core::pose::Pose const & pose,
									 boost::uuids::uuid struct_id,
									 utility::sql_database::sessionOP db_session,
									 core::Size sheet_id,
									 core::Size sheet_antiparallel,
									 core::Size residue_begin,
									 core::Size residue_end);
		
	core::Size
	fill_sw_can_by_components(
							  boost::uuids::uuid	struct_id,
							  utility::sql_database::sessionOP	db_session,
							  core::Size	sw_can_by_components_PK_id_counter,
							  std::string	tag,
							  core::Size	sw_can_by_sh_id,
							  core::Size	sheet_id,
							  core::Size	sheet_antiparellel,
							  core::Size	sw_can_by_components_bs_id,
							  core::Size	strand_is_at_edge,
							  core::Size	residue_begin,
							  core::Size	residue_end);

	utility::vector1<Size>
	get_vec_sw_can_by_sh_id(
							 boost::uuids::uuid struct_id,
							 utility::sql_database::sessionOP db_session);

	core::Size
	get_size_sw_can_by_components_PK_id(
										boost::uuids::uuid struct_id,
										utility::sql_database::sessionOP db_session,
										core::Size sw_can_by_sh_id);

	std::pair<core::Size, core::Size>
	get_starting_res_for_connecting_strands(
											boost::uuids::uuid struct_id,
											utility::sql_database::sessionOP db_session,
											core::Size sw_can_by_sh_id,
											core::Size former_res_end);

	std::pair<core::Size, core::Size>
	get_next_starting_res_for_connecting_strands(
												 boost::uuids::uuid struct_id,
												 utility::sql_database::sessionOP db_session,
												 core::Size sw_can_by_sh_id,
												 core::Size former_ending_res);

	core::Size
	update_sheet_con(
					 boost::uuids::uuid struct_id,
					 utility::sql_database::sessionOP db_session,
					 core::Size sw_can_by_components_PK_id_counter,
					 std::string tag,
					 core::Size sw_can_by_sh_id,
					 bool intra_sheet_con, // if false, then inter_sheet_con
					 core::Size intra_sheet_con_id,
					 core::Size inter_sheet_con_id,
					 core::Size start_res,
					 core::Size end_res);


	bool
	see_whether_other_strands_are_contained(
											boost::uuids::uuid struct_id,
											utility::sql_database::sessionOP db_session,
											core::Size sw_can_by_sh_id,
											core::Size res_num);

	core::Size
	delete_this_sw_can_by_sh_id(
								boost::uuids::uuid struct_id,
								utility::sql_database::sessionOP db_session,
								core::Size sw_can_by_sh_id);
		
	core::Size
	get_segment_id(
				   boost::uuids::uuid struct_id,
				   utility::sql_database::sessionOP db_session,
				   core::Size all_strands_index);
	
	core::Size
	get_num_of_distinct_sheet_id(
								boost::uuids::uuid struct_id,
								utility::sql_database::sessionOP db_session);

private:

	core::Size
	min_num_strands_to_deal_;
	
	core::Size
	max_num_strands_to_deal_;
		
	core::Size
	min_strand_size_;

	core::Real
	min_CA_CA_dis_;

	core::Real
	max_CA_CA_dis_;

	core::Real
	min_O_N_dis_;

	core::Real
	max_O_N_dis_;

	core::Real
	min_C_O_N_angle_;

	core::Real
	min_sheet_dis_;

	core::Real
	max_sheet_dis_;

	core::Real
	min_sheet_angle_;

	core::Real
	max_sheet_angle_;

	core::Real
	min_sheet_torsion_cen_res_;

	core::Real
	max_sheet_torsion_cen_res_;

	core::Real
	min_inter_sheet_dis_CA_CA_;

	core::Real
	max_inter_sheet_dis_CA_CA_;

	bool
	extract_sandwich_;
		
	bool	
	write_chain_B_resnum_;
	
	bool	
	write_phi_psi_;
		
	core::Size	
	max_num_sw_per_pdb_;
		
	

}; // class SandwichFeatures : public protocols::features::FeaturesReporter

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
