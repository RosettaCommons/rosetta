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


#ifndef INCLUDED_protocols_swa_SWA_RNAUtil_HH
#define INCLUDED_protocols_swa_SWA_RNAUtil_HH


#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh> /*For PuckerState and Torsion_Info*/

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <set>

//Auto Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>


typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {

////////////////////////////////Convert these to Rhiju version Jan 29, 2010 Parin S///////////////////////////////////////////////////////////////////

core::Size
get_matching_atom_name(std::string const & atom_name, core::conformation::Residue const & rsd);

void
setup_suite_atom_id_map(core::conformation::Residue const & rsd_1, core::conformation::Residue const & rsd_2, core::Size const res_num, std::set<std::string> const & special_atom_set, std::string mode, core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

void
setup_suite_atom_id_map(core::pose::Pose const & pose_1, core::pose::Pose const & pose_2, core::Size const base_res, bool Is_prepend, core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

void
setup_suite_atom_id_map(core::pose::Pose const & pose_1, core::pose::Pose const & pose_2, core::Size const base_res, core::Size const base_res2, bool Is_prepend, core::id::AtomID_Map< core::id::AtomID > & atom_ID_map);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
output_pair_size(std::pair<core::Size, core::Size> const & pair_size);

void
output_pair_size_vector(utility::vector1 <std::pair<core::Size, core::Size> > const & pair_size_vector, std::string const & output_string);

bool
pair_sort_citeria(std::pair<core::Size, core::Size> pair_one, std::pair<core::Size, core::Size> pair_two);

void sort_seq_num_list(utility::vector1<core::Size> & seq_num_list);

void Output_seq_num_list(std::string const tag, utility::vector1<core::Size> const & seq_num_list, core::Size const spacing=40);

// Undefinded, commenting out to fix PyRosetta build  bool seq_num_list_sort_citeria(core::Size seq_num_1, Residue_info seq_num_2);

void
Sort_pair_list(utility::vector1< std::pair<core::Size, core::Size> > pair_list);

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

//void
//get_bulge_rotamers( utility::vector1< utility::vector1 <core::Real> >& rotamer_list, PuckerState const & pucker1, PuckerState const & pucker2 );

void
Output_pose_data_list(utility::vector1 <pose_data_struct2> const & pose_data_list, std::string const & silent_file , bool const write_score_only = false );

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

core::Real
atom_square_deviation(core::conformation::Residue const & rsd_1, core::conformation::Residue const & rsd_2, core::Size const & atomno_1, core::Size const & atomno_2, bool verbose);

void
suite_square_deviation(core::pose::Pose const & pose1, core::pose::Pose const & pose2, bool const & prepend_res, core::Size const & moving_res_num, core::Size& atom_count, core::Real& sum_sd, bool verbose);

void
Output_title_text(std::string const title);

bool
Check_chain_closable(core::pose::Pose const & pose, core::Size const five_prime_chain_break_res, core::Size const gap_size = 1);

void
Freeze_sugar_torsions(core::kinematics::MoveMap & mm, core::Size const total_residue);


void
o2star_minimize(core::pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn);

utility::vector1<std::string>
Tokenize_with_input_delimiters(std::string const & str, std::string const delimiters);

void
Correctly_position_cutpoint_phosphate_torsions(core::pose::Pose & current_pose, core::Size five_prime_chainbreak,  bool verbose=false);

//void
//apply_rotamer( core::pose::Pose & pose,
//							 utility::vector1< core::id::TorsionID > const & torsion_ids,
//							 utility::vector1< core::Real > const & rotamer_values );

void
apply_rotamer( core::pose::Pose & pose, utility::vector1< Torsion_Info >  const & rotamer_list);

//This is the version currently used for floating base sampling....slightly different from the old version of Base_centroid_screening which appear to StepWiseRNA_Sampling
//Should probably integrate this with Rhiju's class Jan 28, 2010. ***ALERT***RHIJU pointed out that the screening condition is slight different in his new class.
numeric::xyzVector<core::Real>
get_base_centroid( core::conformation::Residue const & rsd , bool verbose=false); ///Stolen this from RNA_CentroidInfo.cc


numeric::xyzMatrix< core::Real >
get_base_coordinate_system( core::conformation::Residue const & rsd, numeric::xyzVector<core::Real> const & centroid ); ///Stolen this from RNA_CentroidInfo.cc

bool
Base_centroid_screening( core::pose::Pose const & pose,
												 utility::vector1< core::Size > moving_positions,
												 utility::vector1 < base_stub > const & other_residues_base_list,
												 SillyCountStruct & count_data);

bool
Base_centroid_screening(base_stub const & moving_res_base, utility::vector1 < base_stub > const & other_residues_base_list, SillyCountStruct & count_data);

PuckerState
Get_residue_pucker_state(core::pose::Pose const & pose, core::Size const seq_num);


}
}
}

#endif
