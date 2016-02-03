// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @brief  some utilities for fragments
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
#ifndef INCLUDED_core_fragment_util_HH
#define INCLUDED_core_fragment_util_HH

// Project Headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <set>

namespace core {
namespace fragment {

core::kinematics::Stub getxform(numeric::xyzVector<core::Real> m1,
	numeric::xyzVector<core::Real> m2,
	numeric::xyzVector<core::Real> m3,
	numeric::xyzVector<core::Real> f1,
	numeric::xyzVector<core::Real> f2,
	numeric::xyzVector<core::Real> f3);

void xform_pose(core::pose::Pose& pose,
	const core::kinematics::Stub& s,
	core::Size sres = 1,
	core::Size eres = 0);

/// @brief Removes all but the top <k> fragments from <fragments>
void retain_top(core::Size k, FragSetOP fragments);

void steal_constant_length_frag_set_from_pose (
	core::pose::Pose const & pose, core::fragment::ConstantLengthFragSet & fragset
);

void steal_frag_set_from_pose(
	pose::Pose const & pose_in,
	FragSet & fragset,
	core::fragment::FragDataCOP frag_type
);

void steal_frag_set_from_pose (
	pose::Pose const& pose_in,
	Size const begin,
	Size const end,
	FragSet& fragset,
	core::fragment::FragDataCOP frag_type
);

void steal_frag_set_from_pose(
	pose::Pose const & pose_in,
	FragSet & fragset,
	core::fragment::FragDataCOP frag_type,
	std::set< core::Size > const& active_residues
);

/// @brief Function for chopping the X-sized fragments in source into fragments that
/// go into N-sized fragments in dest. Sizes are encoded by the value of
/// max_frag_length() in both source and dest.
void chop_fragments(
	core::fragment::FragSet & source, core::fragment::FragSet & dest
);

// undocumented function which passes variables around by reference for no obviou
// reason. Awesome.
void compute_per_residue_coverage(
	core::fragment::FragSet const & frags,
	utility::vector1< core::Size > & nr_frags
);

// create new FragSet ( same type as good_frags ) and fill up with frags from
// filling such that at each residue mn_nr_frags are available bRandom =>
// select fill fragments randomly
FragSetOP merge_frags( FragSet const& good_frags, FragSet const& filling, core::Size min_nr_frags, bool bRandom = true );

void flatten_list( FrameList & frames, FragID_List& frag_ids );

/// @brief goes through all frag_data in the frame and puts the
/// best scoring one into the pose, i.e. quenches the frame
void
apply_best_scoring_fragdata(
	pose::Pose & pose,
	Frame const & frame,
	scoring::ScoreFunction const & sfxn
);

/// @brief writes FragData in multimodel pdb format
/// start_tag can be used to only write a subset
/// of the contained frag data
void dump_frames_as_pdb(
	pose::Pose const & pose,
	utility::vector1< FrameOP > const & frames,
	std::string const & filename,
	core::Size start_frag = 1
);

/// @brief filling a frameset from a multimodel pdb file
/// @brief returns true if no error occured
bool fill_template_frames_from_pdb(
	pose::Pose const & pose,
	utility::vector1< FrameOP > const & template_frames,
	std::string const & filename
);

void read_std_frags_from_cmd( FragSetOP& fragset_large, FragSetOP& fragset_small );

/// @brief given a JumpFrame with Up and DownJumpSRFDs as LAST SRFDs this will make a fold-tree compatible with the
/// Frame...   this is NOT GOOD for sampling, since it introduces cut-points outside of fragments
/// later for sampling: one could probably write a routine that looks if it can move existing Jumps in Fold-tree to
/// fit the FRAME ... if not it returns failure...
void make_simple_fold_tree_from_jump_frame( Frame const&, core::Size total_residue, kinematics::FoldTree& new_fold_tree );

void
fragment_set_slice ( ConstantLengthFragSetOP & fragset, Size const & min_res, Size const & max_res );

void
fragment_set_slice( core::fragment::ConstantLengthFragSetOP & fragset,
	utility::vector1< core::Size > const & slice_res );

void
make_pose_from_frags( pose::Pose & pose, std::string sequence, utility::vector1<FragDataCOP> frags, bool chains = false );

} //fragment
} //core

#endif
