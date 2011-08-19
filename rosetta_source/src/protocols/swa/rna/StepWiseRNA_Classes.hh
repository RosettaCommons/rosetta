// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_Classes_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_Classes_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/GreenPacker.fwd.hh>
#include <protocols/swa/MainChainTorsionSet.hh> // should make a .fwd.hh probably
#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace rna {

struct Residue_info_struct{

	std::string name;
	Size seq_num; //Full_pose_seq_number

};

struct base_struct{
	numeric::xyzVector<Real> centroid;
	Matrix base_coordiate_matrix;
	Residue_info_struct residue_info;
};


struct base_stub{
	numeric::xyzVector<Real> centroid;
	Matrix base_coordiate_matrix;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct pose_data_struct2 {
//	Real rmsd;
//	Size group_rotamer;
//	Size subgroup_rotamer;
//	std::vector< Real> diff_torsions;
 	Real O3_C5_distance; //This is only used in the connect_double_strand_pose_overlap function

	Real score;
	pose::PoseOP pose_OP;
	std::string tag;
};

struct count_struct{
	Size output_pose_count;
	Size good_rep_rotamer_count;
	Size good_atr_rotamer_count;
	Size good_angle_count;
	Size good_distance_count;
	Size C5_O3_distance_count;
	Size Near_Native_Rotamer_Count;
	Size base_pairing_count;
	Size base_stack_count;
	Size both_count;
	Size tot_rotamer_count;
	Size fine_rmsd_count;
	Size rmsd_count;
	Size native_minimize_pass_count;
	Size native_minimize_fail_count;
	Size pass_before_last_res_reb_chain_break_U;
	Size total_before_last_res_reb_chain_break_U;
	Size pass_before_last_res_reb_chain_break_M;
	Size total_before_last_res_reb_chain_break_M;
};


////////////////////////////////

enum PuckerState{ ALL, NORTH, SOUTH };

