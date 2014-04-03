// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Classes.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_Classes_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_Classes_HH


#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>
#include <utility/exit.hh> //April 29, 2011


#define O3I_C5I_PLUS_TWO_MAX_DIST 11.4226  //See ~/minirosetta/04_2010_Early_RUN/test/Apr_9_calculate_theortical_RNA_distance_O3i_C5iplus2/140_220 for detail.
#define O3I_C5I_PLUS_ONE_MAX_DIST 3.968000 //see data below the function get_C4_C3_distance_range() in StepWiseRNA_Classes.hh
#define O3I_O3I_PLUS_ONE_MAX_DIST 7.45583 //see ~/minirosetta/test/Sept_19_calculate_MAX_O3_O3_distance/trail_1_range_0_360_bin_size_5/output.txt
#define C5I_C5I_PLUS_ONE_MAX_DIST 7.71355 //see ~/minirosetta/test/Sept_19_calculate_MAX_C5_C5_distance/trail_3_0_torsion_range_360_5_degree_bin/output.txt
#define O3I_C5I_MIN_DIST 2.000
#define O3I_C5I_MAX_DIST 4.627

#define GAP_SIZE_DUMMY 999 //encodes for an 'infinite' gap between residues, e.g., if on different chains.

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {


class Jump_point{

		public:

		Jump_point():
		five_prime_seq_num( 0 ),
		five_prime_atom( "O3'" ),
		three_prime_seq_num( 0 ),
		three_prime_atom( "P" ),
		cut_point( 0 )
	{
	}

	~Jump_point(){};

	public:

	core::Size five_prime_seq_num;
	std::string five_prime_atom;
	core::Size three_prime_seq_num;
	std::string three_prime_atom;
	core::Size cut_point; //Choose the residue five_prime of the actual cutpoint position

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Torsion_Info{
	core::id::TorsionID id;
	core::Real value;
};

struct output_data_struct {

	core::Real rmsd_wrt_correct_minimized;
	core::Real rmsd_wrt_minimized;
	core::Real rmsd;
	core::Real loop_rmsd_wrt_correct_minimized;
	core::Real loop_rmsd_wrt_minimized;
	core::Real loop_rmsd;
	core::Real diff_torsions;
	core::Real current_score;
	core::Real O3_C5_distance;
};

class InternalWorkingResidueParameter{

public:

	InternalWorkingResidueParameter():
		fake_working_moving_suite( 0 ),
		possible_working_res_1( 0 ),
		possible_working_res_2( 0 )
		{
	}

	~InternalWorkingResidueParameter(){};

	public:

		core::Size fake_working_moving_suite;
		core::Size possible_working_res_1;
		core::Size possible_working_res_2;
};





class StepWiseRNA_CountStruct{

public:

  StepWiseRNA_CountStruct():
    output_pose_count( 0 ),
    full_score_count( 0 ),
    good_bin_rep_count( 0 ),
    good_rep_rotamer_count( 0 ),
    good_atr_rotamer_count( 0 ),
    good_angle_count( 0 ),
    good_distance_count( 0 ),
    chain_break_screening_count( 0 ),
    chain_closable_geometry_count( 0 ),
    chain_closable_geometry_count2( 0 ),
    chain_closable_geometry_count3( 0 ),
    Near_Native_Rotamer_Count( 0 ),
    base_pairing_count( 0 ),
    base_stack_count( 0 ),
    strict_base_pairing_count( 0 ),
    pass_base_centroid_screen( 0 ),
    both_count( 0 ),
    tot_rotamer_count( 0 ),
    fine_rmsd_count( 0 ),
    rmsd_count( 0 ),
    non_clash_sugar( 0 ),
    fast_full_atom_VDW_repulsion_screen( 0 ),
    in_range_CCD_torsion( 0 ),
    total_bin( 0 ),
    bulge_at_chain_closure_count( 0 ),
    before_chain_break_grid_index_screening( 0 ),
    chain_break_grid_index_screening( 0 ),
    residues_contact_screen( 0 ),
    test_count_one( 0 ),
    test_count_two( 0 )
    {
  }


  ~StepWiseRNA_CountStruct(){};

public:

  core::Size output_pose_count;
  core::Size full_score_count;
  core::Size good_bin_rep_count;
  core::Size good_rep_rotamer_count;
  core::Size good_atr_rotamer_count;
  core::Size good_angle_count;
  core::Size good_distance_count;
  core::Size chain_break_screening_count;
  core::Size chain_closable_geometry_count;
  core::Size chain_closable_geometry_count2;
  core::Size chain_closable_geometry_count3;
  core::Size Near_Native_Rotamer_Count;
  core::Size base_pairing_count;
  core::Size base_stack_count;
  core::Size strict_base_pairing_count;
  core::Size pass_base_centroid_screen;
  core::Size both_count;
  core::Size tot_rotamer_count;
  core::Size fine_rmsd_count;
  core::Size rmsd_count;
  core::Size non_clash_sugar;
  core::Size fast_full_atom_VDW_repulsion_screen;
  core::Size in_range_CCD_torsion;
  core::Size total_bin;
  core::Size bulge_at_chain_closure_count;
  core::Size before_chain_break_grid_index_screening;
  core::Size chain_break_grid_index_screening;
  core::Size residues_contact_screen;
  core::Size test_count_one;
  core::Size test_count_two;
};

////////////////////////////////////////////////////////////
} //rna
} //sampling
} //stepwise
} //protocols

#endif
