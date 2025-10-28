// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/floppy_tail/FloppyTail_publication.cc
/// @brief FloppyTail extra functions from original publication - this calculates some statistics used for the first published use of the code (Kleiger G, Saha A, Lewis S, Kuhlman B, Deshaies RJ. Rapid E2-E3 assembly and disassembly enable processive ubiquitylation of cullin-RING ubiquitin ligase substrates. Cell. 2009 Nov 25;139(5):957-68. PubMed PMID: 19945379.)
/// @author Steven Lewis smlewi@gmail.com

// Unit Headers
#include <protocols/floppy_tail/FloppyTail_publication.hh>

// Project Headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>

//JD headers
#include <protocols/jd2/util.hh>

namespace protocols {
namespace floppy_tail {

/// @details This function is specific to the original system for which this code was written
void create_extra_output( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP score_fxn ){

	//score main pose
	core::Real const full_score((*score_fxn)(pose));

	//constants for hacking up pose
	core::Size const E2_chain(pose.conformation().num_chains()); //E2 is last chain
	core::Size const E2_start(pose.conformation().chain_begin(E2_chain));

	//create and score E3_RING portion
	core::pose::Pose E3_RING(pose);
	E3_RING.conformation().delete_residue_range_slow(E2_start, E3_RING.size());
	core::Real const E3_RING_score((*score_fxn)(E3_RING));

	//create and score E2 portion
	core::pose::PoseOP E2 = pose.split_by_chain(E2_chain);
	core::Real const E2_score((*score_fxn)(*E2));

	//debugging - check these pdbs
	//protocols::jd2::output_intermediate_pose( E2, "E2");
	//protocols::jd2::output_intermediate_pose( E3_RING, "E3_RING");

	//print Binding_Energy
	core::Real const BE(full_score - E3_RING_score - E2_score);
	protocols::jd2::add_string_to_current_job("Binding_Energy = Complex - E2 alone - E3/RING alone");
	protocols::jd2::add_string_real_pair_to_current_job("Binding energy", BE);
	//print crosslink distance
	//magic numbers: crosslink experiment tested these residues
	bool const NEDD8(E2_chain == 4);
	core::Size const cl_base(pose.pdb_info()->pdb2pose("A", (NEDD8 ? 1676: 679)));   //K1676??
	core::Size const cl_tail(pose.pdb_info()->pdb2pose("C", 223));

	core::Real const cl_dist(pose.residue(cl_base).atom("CA").xyz().distance( pose.residue(cl_tail).atom("CA").xyz() ));
	protocols::jd2::add_string_to_current_job("Crosslink = distance between atom CA on residues (chain A, K679/K1676; chain C, 223)");
	protocols::jd2::add_string_real_pair_to_current_job("crosslink distance", cl_dist);

	//std::cout << full_score << " " << E3_RING_score << " " << E2_score << std::endl;

	return;
}

} //floppy_tail
} //protocols
