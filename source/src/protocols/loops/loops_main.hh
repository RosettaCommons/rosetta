// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Chu Wang

#ifndef INCLUDED_protocols_loops_loops_main_hh
#define INCLUDED_protocols_loops_loops_main_hh

// fwd declaration
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <vector>

#include <utility/vector1.hh>


namespace protocols {
namespace loops {

/// @brief the main function for perform loop modeling
////voi/d
///loops_main( core::pose::Pose & pose );

/// @brief construct a fold tree from loop definition
void
fold_tree_from_loops(
	core::pose::Pose const & pose,
	Loops const & loops,
	core::kinematics::FoldTree & f,
	bool terminal_cutpoint = false // should terminal loops respect the given cutpoint
);

//////////////////////////////////////////////////////////////////////////////////
/// @details Set the fold tree to contain a single chain break at the loops' position
void set_single_loop_fold_tree(
	core::pose::Pose & pose,
	Loop const & loop
);

/// @details  Remove cutpoint variants
void remove_cutpoint_variants(
	core::pose::Pose & pose,
	bool force = false
);

//////////////////////////////////////////////////////////////////////////////////
/// @brief  Add cutpoint variants to all the cutpoints in a Pose.
void add_cutpoint_variants( core::pose::Pose & pose );

/// @brief  Add cutpoint variants around a single cutpoint (defined by a Loop).
void add_single_cutpoint_variant( core::pose::Pose & pose, const Loop & loop );

/// @brief  Add cutpoint variants around a single cutpoint (defined by residue number).
void add_single_cutpoint_variant( core::pose::Pose & pose, const core::uint cutpoint );


/// @brief create a MoveMap for use of minimization based on loop definition (wrapper)
void
loops_set_move_map(
	core::pose::Pose & pose,
	Loops const & loops,
	bool const fix_template_sc,
	core::kinematics::MoveMap & mm,
	core::Real neighbor_dist,
	bool const allow_omega_move,
	bool const allow_takeoff_torsion_move
);

void
loops_set_move_map(
	core::pose::Pose & pose,
	Loops const & loops,
	bool const fix_template_sc,
	core::kinematics::MoveMap & mm,
	core::Real neighbor_dist = 10.0
);

/// @brief create a MoveMap for use of minimization based on loop definition
void
loops_set_move_map(
	Loops const & loops,
	utility::vector1<bool> const & allow_sc_move,
	core::kinematics::MoveMap & mm
);

void
loops_set_move_map(
	Loops const & loops,
	utility::vector1<bool> const & allow_sc_move,
	core::kinematics::MoveMap & mm,
	bool const allow_omega_move,
	bool const allow_takeoff_torsion_move
);

/// @brief Create a new MoveMapOP for use in minimizing the given loop.
core::kinematics::MoveMapOP
move_map_from_loops(
	core::pose::Pose & pose,
	Loops const & loops,
	bool const fix_template_sc,
	core::Real neighbor_dist = 10.0,
	bool const flanking_residues = false
);

/// @brief Create a new MoveMapOP for use in minimizing the given loop.
core::kinematics::MoveMapOP
move_map_from_loop(
	core::pose::Pose & pose,
	Loop const & loop,
	bool const fix_template_sc,
	core::Real neighbor_dist = 10.0,
	bool const flanking_residues = false
);

void
set_move_map_for_centroid_loop(
	Loop const & loop,
	core::kinematics::MoveMap & mm
);

/// @brief add flank stem residues to the loop movemap
void //made by JQX
add_loop_flank_residues_bb_to_movemap(
	Loops const & loops,
	core::kinematics::MoveMap & mm,
	core::Size flank_size=2
);

/// @brief close loops by the CCD mechanism
void
ccd_close_loops(
	core::pose::Pose & pose,
	Loops const & loops,
	core::kinematics::MoveMap const& mm
);




/// @brief mark loop residues and its neighbors as necessary in a sequence map.
/// @details
///  Uses 10A neighbor graph for neighbors
///  THEN takes distance from selection to neighbors to trim neighbors.
///  Excludes disulfide residues.
void select_loop_residues(
	core::pose::Pose const & pose,
	Loops const & loops,
	bool const include_neighbors,
	utility::vector1<bool> & map,
	core::Real neighbor_dist = 10.0
);

/// @brief mark loop residues and its neighbors as necessary in a sequence map.
/// @details
///  Uses 10A neighbor graph for neighbors
///  THEN takes distance from selection to neighbors to trim neighbors.
///  Excludes disulfide residues.
utility::vector1<bool> select_loop_residues(
	core::pose::Pose const & pose,
	Loops const & loops,
	bool const include_neighbors,
	core::Real neighbor_dist = 10.0
);

/// @brief mark loop residues and its neighbors as necessary for one loop.
/// @details
///  Uses 10A neighbor graph for neighbors
///  THEN takes distance from selection to neighbors to trim neighbors.
///  Excludes disulfide residues.
void select_loop_residues(
	core::pose::Pose const & pose,
	Loop const & loop,
	bool const include_neighbors,
	utility::vector1<bool> & map,
	core::Real neighbor_dist = 10.0
);

/// @brief mark loop residues and its neighbors as necessary for one loop.
/// @details
///  Uses 10A neighbor graph for neighbors
///  THEN takes distance from selection to neighbors to trim neighbors.
///  Excludes disulfide residues.
utility::vector1<bool> select_loop_residues(
	core::pose::Pose const & pose,
	Loop const & loop,
	bool const include_neighbors,
	core::Real neighbor_dist = 10.0
);



/// @brief filter set of loop neighbors to a certain CB distance
///  Takes distance from selection to neighbors to trim neighbors.
void filter_loop_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> & map,
	Loops const & loops,
	core::Real & dist_cutoff
);




/// @brief helper function to set secondary structure of a Pose from an external
/// file.
bool
set_secstruct_from_psipred_ss2(
	core::pose::Pose & pose,
	std::string const & filename = std::string()
);

/// @brief another helper function to set secondary structure of a Pose from an external file.
bool
set_secstruct_from_dssp(
	core::pose::Pose & pose,
	std::string const & filename
);


/// @details   set ideal BB geometry; this must occur so that loops with missing density work.
void idealize_loop(
	core::pose::Pose & pose,
	Loop const & loop
);

/// @details  Set a loop to extended torsion angles.
void set_extended_torsions(
	core::pose::Pose & pose,
	Loop const & loop
);


void read_loop_fragments(
	std::vector< core::fragment::FragSetOP > &frag_libs
);

void read_loop_fragments(
	utility::vector1< core::fragment::FragSetOP > &frag_libs
);


//////////////////////////////////////////////////////////////////////////////////
/// @details  Rebuild a loop via fragment insertion + ccd closure + minimization
void remove_missing_density(
	core::pose::Pose & pose,
	Loop const & loop
);

core::Real native_loop_core_CA_rmsd(
	const core::pose::Pose & native_pose,
	const core::pose::Pose & pose,
	loops::Loops loops,
	int &corelength
);

/// @brief calculate rmsd of loop residues with repect to native (template aligned)
core::Real
loop_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	Loops const & loops,
	bool CA_only = false,
	bool bb_only = true
);

/// @brief As above but actuall superimposes the non-loop part
core::Real
loop_rmsd_with_superimpose(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	Loops const & loops,
	bool CA_only = false,
	bool bb_only = true
);

/// @brief As above but actually superimposes only the core part (in case there are multiple loops...)
core::Real
loop_rmsd_with_superimpose_core(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	Loops const & loops,
	Loops const & core,
	bool CA_only = false,
	bool bb_only = true
);


/// @brief calculate rmsd of loop residues with repect to native (loop fit)
core::Real
loop_local_rmsd(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	Loops const & loops
);

/// @brief  Given a sequence mapping which may have simple indels, trim back around those indels so that the loops can plausibly be closed.

void
trim_back_sequence_mapping(
	core::id::SequenceMapping & mapping,
	std::string const & source_seq,
	std::string const & target_seq,
	core::Size const min_loop_size
);


void
extend_sequence_mapping(
	core::pose::Pose const & pose,
	core::id::SequenceMapping & mapping,
	std::string & source_seq,
	std::string & target_seq
);

void
set_loop_cutpoint_in_pose_fold_tree(
	core::Size const new_cutpoint,
	core::pose::Pose & pose,
	core::Size const loop_begin,
	core::Size const loop_end
);

void
apply_sequence_mapping(
	core::pose::Pose & pose,
	std::string const & target_seq,
	core::id::SequenceMapping const & start_mapping
);

} //namespace loops
} //namespace protocols

#endif
