// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/util.hh
/// @brief Utility functions for Pose Sewing.
/// @author Frank Teets (frank.teets@proteininnovation.org)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_pose_sewing_util_hh
#define INCLUDED_protocols_pose_sewing_util_hh

#include <protocols/pose_sewing/strong_types.hh>
#include <protocols/sewing/scoring/MotifScorer.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/conformation/Residue.hh>

#include <map>
#include <set>

namespace protocols {
namespace pose_sewing {

////////////////////////////////////
/////// Blocks and Elements ////////
////////////////////////////////////

///@brief Get residues in each block through a set of equal subsets.
/// Number of subsets are the number of SS blocks

void
identify_ss_blocks(std::set<core::select::residue_selector::ResidueSubset> & outset, utility::vector1< bool > const & subset );

///@brief Will eventually replace standard function.
void
identify_ss_blocks_vec (utility::vector1< core::select::residue_selector::ResidueSubset> & out_vec, utility::vector1< bool > const & subset );

///@brief Assign E/H blocks in the given pose
void
calculate_blocks(std::map< core::Size, core::Size > & outmap, core::pose::Pose const & pose );

///@brief Assign arbitrary blocks from a Residue subset.
void
calculate_blocks_from_subset(std::map< core::Size, core::Size > & outmap, utility::vector1< bool > const & subset);

void
calculate_helices(std::map< core::Size, core::Size > & outmap, core::pose::Pose const & pose, core::Size min_length);


///@brief Is the DSSP all loop?
/// A pose unfortunately starts with all L as the DSSP when creating a pose from seq and inserting residues as
/// is done for segments.
bool
all_L_dssp(core::pose::Pose const & pose);


///////////////////////////////////////////////
/////////////// Motif Scores //////////////////
///////////////////////////////////////////////

/// @brief Calculate the motif score between two residues.
core::Real
calculate_motif_score_bt_residues(
	sewing::scoring::MotifScorer const & scorer,
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	char res1_SS,
	char res2_SS
);
core::Real
calculate_distance_score_bt_residues(
	core::conformation::Residue const & N_res,
	core::conformation::Residue const & C_res,
	core::Real min_score,
	core::Real min_distance,
	core::Real dist_mult
);

///@brief Core code to calculate BlockWise/Elementwise Pose Compatable Motif Metrics.
/// Reports worst score found for each block (Number as a string for the metric)
/// If max_pair_score is something not 0, will use it (essentially used for filtering outside of metrics).
void
calculate_bw_pose_compat_motifs(
	std::map< std::string, core::Real > & outmap,
	core::pose::Pose const & pose,
	std::set<core::select::residue_selector::ResidueSubset> const & block_selections,
	bool drop_best = true,
	bool normalize_by_residues = false,
	core::Real max_pair_score = 0,
	bool use_motifs = true
);

///@brief Core code to calculate window metrics.
/// Returns worst score for each window, represented as a string for the metric.
void
calculate_bw_window_motifs(
	std::map< std::string, core::Real > & outmap,
	core::pose::Pose const & pose,
	std::set< core::select::residue_selector::ResidueSubset> const & block_selections,
	core::Size window_width = 3,
	bool use_motifs = true
);

} //pose_sewing
} //protocols


#endif //protocols/pose_sewing_util_hh

