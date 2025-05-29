// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/recon_util.cc
/// @brief Utility functions for the recon_design namespace
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <protocols/recon_design/recon_util.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>


namespace protocols {
namespace recon_design {

/// Based on a pose and a resfile, get the indices of all designable residues
utility::vector1< core::Size >
get_designable_residues( core::pose::Pose & pose, std::string resfile ) {
	using namespace core::pack::task;

	PackerTaskOP design_task = TaskFactory::create_packer_task( pose);
	parse_resfile( pose, *design_task, resfile );

	utility::vector1< core::Size > designable;
	utility::vector1<bool> designing = design_task->designing_residues();
	for ( core::Size i = 1; i <= designing.size(); ++i ) {
		if ( designing[ i ] ) {
			designable.push_back( i );
		}
	}

	return designable;
}

/// Based on a pose and the indices of all designable residues, get a
/// string of all designable AAs concatenated
utility::vector1< std::string >
get_designable_sequence ( core::pose::Pose & pose, utility::vector1< core::Size > designable_residues ) {
	//std::string sequence = "";
	utility::vector1< std::string > sequence;
	for ( core::Size seqpos: designable_residues ) {
		//sequence += pose.residue( seqpos ).name1();
		sequence.push_back(pose.residue_type( seqpos ).base_name());
	}
	return sequence;
}

/// Based on a list of sequences from poses, get all the AAs present at
/// the position given by position_no
utility::vector1< std::string >
get_candidate_AAs( utility::vector1< utility::vector1< std::string > > const & other_pose_sequences,
	core::Size position_no ) {

	utility::vector1< std::string > candidate_AAs;

	// Iterate through all poses in my collection and find their AA at this position
	//for ( std::string const & sequence: other_pose_sequences ) {
	// char current_AA = sequence[ position_no-1 ]; // string is zero indexed
	// std::string current_AA_3letter = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( current_AA ) );
	// if ( std::find( candidate_AAs.begin(), candidate_AAs.end(), current_AA_3letter ) == candidate_AAs.end() ) {
	//  candidate_AAs.push_back( current_AA_3letter );
	// }
	//}
	for ( utility::vector1< std::string > const & sequence: other_pose_sequences ) {
		std::string current_AA = sequence[ position_no ]; // vector1 is 1 indexed
		//std::string current_AA_3letter = core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( current_AA ) );
		if ( std::find( candidate_AAs.begin(), candidate_AAs.end(), current_AA ) == candidate_AAs.end() ) {
			candidate_AAs.push_back( current_AA );
		}
	}

	return candidate_AAs;
}

/// Given a list of poses, find the index of a particular pose
core::Size
find_pose_in_vector( core::pose::Pose const & pose,
	utility::vector1< core::pose::PoseOP > & poses ) {

	core::Size current_pose = 0;
	core::Real ref_pose_index;
	bool return_value = core::pose::getPoseExtraScore(  pose, "msd_job_dist_index", ref_pose_index );
	if ( !return_value ) {
		utility_exit_with_message("Error: poses are not indexed correctly. "
			"If you run MSDMover from RosettaScripts you need to pass the -run:msd_job_dist flag");
	}

	for ( core::Size ii = 1; ii <= poses.size(); ++ii ) {
		core::Real current_pose_index;
		core::pose::getPoseExtraScore( *poses[ ii ], "msd_job_dist_index", current_pose_index );

		if ( ref_pose_index == current_pose_index ) {
			current_pose = ii;
			break;
		}
	}

	if ( current_pose == 0 ) {
		utility_exit_with_message("Error: pose not found in vector. Strange...");
	}

	return current_pose;
}

} // namespace recon_design
} // namespace protocols

