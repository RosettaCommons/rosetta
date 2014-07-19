// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_WriteToDBFromSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_WriteToDBFromSandwichFeatures_HH

//Devel
#include <protocols/features/strand_assembly/CheckForSandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>


using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {

	core::Size
	WriteToDB_avg_b_factor_CB_at_each_component (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		core::Size	sw_can_by_sh_id);



	core::Size
	WriteToDB_hydrophobic_ratio_net_charge (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size	sw_can_by_sh_id);

	core::Size
	WriteToDB_min_dis_between_sheets_by_all_res (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size	sw_can_by_sh_id,
		core::pose::Pose & dssp_pose,
		utility::vector1<core::Size>	all_distinct_sheet_ids);



	void
	WriteToDB_turn_AA(
		core::pose::Pose const & pose,
		core::Size sw_can_by_sh_id,
		core::Size i,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string turn_type);


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* WriteToDBFromSandwichFeatures_HH_ */
