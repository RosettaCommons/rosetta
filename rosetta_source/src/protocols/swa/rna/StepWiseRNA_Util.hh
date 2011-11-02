// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_Util_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_Util_hh


#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.hh> /*For PuckerState*/

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>




typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {



class Jump_point{

		public:

		Jump_point():
		five_prime_seq_num( 0 ),
		five_prime_atom("O3*"),
		three_prime_seq_num( 0 ),
		three_prime_atom("P"),
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

// Why don't we just keep this in a "regular" kinematics::Stub?
struct base_stub{
  numeric::xyzVector<core::Real> centroid;
  Matrix base_coordinate_matrix;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct pose_data_struct2{
	core::Real score;
	core::pose::PoseOP pose_OP;
	std::string tag;
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




class SillyCountStruct{

public:

  SillyCountStruct():
    output_pose_count( 0 ),
    good_rep_rotamer_count( 0 ),
    good_atr_rotamer_count( 0 ),
    good_angle_count( 0 ),
    good_distance_count( 0 ),
    C5_O3_distance_count( 0 ),
    Near_Native_Rotamer_Count( 0 ),
    base_pairing_count( 0 ),
    base_stack_count( 0 ),
    both_count( 0 ),
    tot_rotamer_count( 0 ),
    fine_rmsd_count( 0 ),
    rmsd_count( 0 )
    {
  }


  ~SillyCountStruct(){};

public:

  core::Size output_pose_count;
  core::Size good_rep_rotamer_count;
  core::Size good_atr_rotamer_count;
  core::Size good_angle_count;
  core::Size good_distance_count;
  core::Size C5_O3_distance_count;
  core::Size Near_Native_Rotamer_Count;
  core::Size base_pairing_count;
  core::Size base_stack_count;
  core::Size both_count;
  core::Size tot_rotamer_count;
  core::Size fine_rmsd_count;
  core::Size rmsd_count;
};



////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void
//Output_data_parin(std::ofstream& outfile, core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, output_data_struct & output_data, core::pose::Pose const & current_pose, std::string const & tag, core::Size const reb_res, bool const prepend, bool const write_score_only);

void
output_pair_size(std::pair<core::Size, core::Size> const & pair_size);

void
output_pair_size_vector(utility::vector1 <std::pair<core::Size, core::Size> > const & pair_size_vector, std::string const & output_string);

bool
pair_sort_citeria(std::pair<core::Size, core::Size> pair_one, std::pair<core::Size, core::Size> pair_two);

void
Sort_pair_list(utility::vector1< std::pair<core::Size, core::Size> > pair_list);

bool
Is_close_chain_break(core::pose::Pose const & pose);

core::Size
Get_five_prime_chain_break(core::pose::Pose const & pose);

void
Add_harmonic_chainbreak_constraint(core::pose::Pose & pose, core::Size const five_prime_res);

void
Output_fold_tree_info(core::kinematics::FoldTree const & fold_tree, std::string const pose_name);

void
Output_fold_tree_info(core::pose::Pose const & pose, std::string pose_name);

void
output_rotamer(utility::vector1 <core::Real > & rotamer);

void
get_bulge_rotamers( utility::vector1< utility::vector1 <core::Real> >& rotamer_list, PuckerState const & pucker1, PuckerState const & pucker2 );

void
Output_pose_data_list(utility::vector1 <pose_data_struct2> const & pose_data_list, std::string const & silent_file , bool const write_score_only);

void
Output_data(core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, std::string const & tag, bool const write_score_only, core::pose::Pose const & pose, core::pose::PoseCOP native_poseOP, core::Size const moving_base_residue, bool const Is_prepend);

void
Add_virtual_O2Star_hydrogen(core::pose::Pose & pose);

bool
Remove_virtual_O2Star_hydrogen(core::pose::Pose & pose);

core::Real
suite_rmsd(core::pose::Pose const & pose1,core::pose::Pose const & pose2, core::Size const & seq_num, bool const prepend_res); //This used to be call suite_rmsd_any_seq_num

core::Real
rmsd_over_residue_list(core::pose::Pose const & pose1, core::pose::Pose const & pose2, utility::vector1 <Residue_info> const & residue_list, std::map< core::Size, core::Size > const & res_map, std::map< core::Size, bool > const & Is_prepend_map, bool const verbose);

void
Print_heavy_atoms(core::Size moving_res_num, core::pose::Pose const & pose1, core::pose::Pose const & pose2);

core::Size
Get_num_side_chain_atom_from_res_name(core::chemical::AA const & res_aa, bool const verbose);

void
suite_square_deviation(core::pose::Pose const & pose1, core::pose::Pose const & pose2, bool const & prepend_res, core::Size const & moving_res_num, core::Size& atom_count, core::Real& sum_sd, bool verbose);

void
Output_title_text(std::string const title);

bool
Check_chain_closable(core::pose::Pose const & pose, core::Size const five_prime_chain_break_res, core::Size const gap_size = 1);

void
Freeze_sugar_torsions(core::kinematics::MoveMap & mm, core::Size const total_residue);

void
Output_boolean(bool boolean);

void
Output_movemap(core::kinematics::MoveMap const & mm, core::Size const total_residue);

void
o2star_minimize(core::pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn);

utility::vector1<std::string>
Tokenize_with_input_delimiters(std::string const & str, std::string const delimiters);

void
Correctly_position_cutpoint_phosphate_torsions(core::pose::Pose & current_pose, core::Size five_prime_chainbreak,  bool verbose=false);

core::Size
make_cut_at_moving_suite( core::pose::Pose & pose, core::Size const & moving_suite );

void
apply_rotamer( core::pose::Pose & pose,
							 utility::vector1< core::id::TorsionID > const & torsion_ids,
							 utility::vector1< core::Real > const & rotamer_values );


}
}
}

#endif
