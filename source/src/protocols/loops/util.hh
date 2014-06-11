// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loops_util_hh
#define INCLUDED_protocols_loops_util_hh

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Package headers
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <utility/vector1.hh>
#include <string>

namespace protocols {
namespace loops {

//@brief adds coord constraints for the atoms that are not in the loops structure
void
fix_with_coord_cst( loops::Loops const& rigid, core::pose::Pose& pose, bool bCstAllAtom, utility::vector1< core::Real >& );

///@brief get frags that are fully within the Loop --- shorten(=true/false) frags that are close to the end of loops.
extern void
select_loop_frags(
	Loops const& loops,
	core::fragment::FragSet& source,
	core::fragment::FragSet& loop_frags,
	core::Size min_size = 1 /* set to 0 if you don't want to shorten at all */
);

void
set_extended_torsions_and_idealize_loops( core::pose::Pose& pose, loops::Loops loops );

/// @brief Identical to set_extended_torsions_and_idealize_loops() without the irrational
/// behavior surrounding empty loops.
void
safe_set_extended_torsions_and_idealize_loops(const protocols::loops::Loops& loops,
                                                   core::pose::Pose* pose);

void
addScoresForLoopParts(
	core::pose::Pose& pose,
	Loops loops,
	const core::scoring::ScoreFunction &scorefxn,
	core::pose::Pose& native_pose, core::Size nloops
);

loops::Loops
compute_ss_regions(
	core::Real max_loop_frac,
	core::Size min_length,
	core::fragment::SecondaryStructure const & ss
);

core::scoring::ScoreFunctionOP
get_fa_scorefxn();
core::scoring::ScoreFunctionOP
get_cen_scorefxn();

void
add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_target_pose,  Loops &exclude_regions );

/// loop_str has the format: start:end:cut,start:end:cut and can use rosetta or pdb numbering. The return value is an owning pointer to a Loops object encoding that loop
LoopsOP
loops_from_string( std::string const loop_str, core::pose::Pose const & pose );

// this function will return a bunch of "loops" that refer to residues that are considered part of the core:
// not scored are loops with 4 or more residues, short helices (<=5) that terminate a loop are not scored, too
void
define_scorable_core_from_secondary_structure( core::fragment::SecondaryStructure const&, Loops& score_core );

/// @brief Extract secondary structure chunks from the pose, using multiple secondary structure types
/// this function requires that the pose object already have secstruct information
/// to get this information from structure (DSSP), call
/// protocols::jumping::Dssp dssp_obj( *pose );	dssp_obj.insert_ss_into_pose( *pose );
/// or from secondary structure prediction (psipred_ss2 file), call
///	core::pose::read_psipred_ss2_file(pose);
Loops
extract_secondary_structure_chunks(core::pose::Pose const & pose,
														   std::string extracted_ss_types = "HE",
														   core::Size gap_size = 1,
														   core::Size minimum_length_of_chunk_helix = 5,
														   core::Size minimum_length_of_chunk_strand = 3,
														   core::Real CA_CA_distance_cutoff = 4);

/// @brief Extract secondary structure chunks from the pose, using a given secondary structure type
Loops
extract_secondary_structure_chunks(core::pose::Pose const & pose,
														   char const extracted_ss_type);

Loops
split_by_resSeq(core::pose::Pose const & pose);
Loops
split_by_resSeq(core::pose::Pose const & pose, Loops const & input_chunks);

Loops
find_non_protein_chunks(core::pose::Pose const & pose);

// TODO(cmiles) deduplicate
/// @brief Split into separate chunks if CA-CA distance is over the cutoff
Loops
split_by_ca_ca_dist(core::pose::Pose const & pose,
		Loops const & input_chunks,
		core::Real const CA_CA_distance_cutoff = 4);

/// @brief If two chunks are separated by a small gap of size <= <gap_size>, combine them
Loops
remove_small_gaps(Loops const & input_chunks, core::Size gap_size = 1);

/// @brief Remove small chunks
Loops
remove_short_chunks(Loops const & input_chunks, core::Size min_length = 3);

Loops
extract_continuous_chunks(core::pose::Pose const & pose,
		core::Size const minimum_size = 3,
		core::Real const CA_CA_distance_cutoff = 4);


///@brief Measures C-N distance and N-Ca-C, C-N-CA bond angles to detect chainbreak/poory resolved loop residues if outside cutoff.
///  Returns boolean and the position of found wonky peptide bond (0 if none found).
///  Does not use AtomTree to get angles/lengths, assumes resnum is connected to resnum +1
/// @details Values are based on the CDL.  No peptide bond without severe chainbreaks or missing residues should have values
///  out of this range. -JAB
///  Berkholz DS, Shapovalov MV, Dunbrack RL Jr, Karplus PA (2009)
///  Conformation dependence of backbone geometry in proteins. Structure 17: 1316-1325.
///  deviation is +/- from the min/max values of the CDL.
///
/// @param Stringent cutoffs are: 1.5, 15, 15
///
std::pair<bool, core::Size>
has_severe_pep_bond_geom_issues(
	core::pose::Pose const & pose,
	const Loop & loop,
	bool check_bonds = true,
	bool check_angles = true,
	core::Real max_c_n_dis = 2.0,
	core::Real allowed_ca_c_n_deviation = 25.0,
	core::Real allowed_c_n_ca_deviation = 25.0);

///@brief Measures C-N distance and N-Ca-C, C-N-CA bond angles to detect chainbreak/poory resolved loop residue outside cutoff.
///  Returns boolean and the position of found wonky peptide bond (0 if none found).
///  Does not use AtomTree to get angles/lengths, assumes resnum is connected to resnum +1
/// @details Values are based on the CDL.  No peptide bond without severe chainbreaks or missing residues should have values
///  out of this range. -JAB
///  Berkholz DS, Shapovalov MV, Dunbrack RL Jr, Karplus PA (2009)
///  Conformation dependence of backbone geometry in proteins. Structure 17: 1316-1325.
///  deviation is +/- from the min/max values of the CDL.
///
/// @param Stringent cutoffs are: 1.5, 15, 15
///
std::pair<bool, core::Size>
has_severe_pep_bond_geom_issues(
	core::pose::Pose const & pose,
	core::Size resnum,
	bool check_bonds = true,
	bool check_angles = true,
	core::Real max_c_n_dis = 2.0,
	core::Real allowed_ca_c_n_deviation = 25.0,
	core::Real allowed_c_n_ca_deviation = 25.0);
	
	
} //loops
} //protocols

#endif





