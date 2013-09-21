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
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/angle.functions.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <string>
#include <map>

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {

bool Is_OP2_atom( std::string const & atom_name );
bool Is_OP1_atom( std::string const & atom_name );
bool Is_P_atom( std::string const & atom_name );
bool Is_O2star_atom( std::string const & atom_name );
bool Is_O3star_atom( std::string const & atom_name );
bool Is_O4star_atom( std::string const & atom_name );
bool Is_O5star_atom( std::string const & atom_name );
bool Is_C2star_atom( std::string const & atom_name );
bool Is_C3star_atom( std::string const & atom_name );
bool Is_C4star_atom( std::string const & atom_name );
bool Is_C5star_atom( std::string const & atom_name );
bool Is_1H5star_atom( std::string const & atom_name );
bool Is_2H5star_atom( std::string const & atom_name );
bool Is_H3star_atom( std::string const & atom_name );
bool Is_H4star_atom( std::string const & atom_name );
bool Is_three_prime_phosphate_atom( std::string const & atom_name );
bool Is_five_prime_phosphate_atom( std::string const & atom_name );
bool Is_phosphate_atom( std::string const & atom_name );


void
minimize_with_constraints( core::pose::Pose & pose, core::kinematics::MoveMap const & mm, core::scoring::ScoreFunctionOP const & scorefxn, core::optimization::MinimizerOptions const & options );

bool
check_can_prepend( utility::vector1< core::Size > const & seq_num_list );

bool
check_can_append( utility::vector1< core::Size > const & seq_num_list );

void
apply_protonated_H1_adenosine_variant_type( core::pose::Pose & pose, core::Size const & seq_num, bool const apply_check = true );

void
apply_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num, bool const apply_check = true );

void
apply_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num, utility::vector1< core::Size > const & working_cutpoint_closed_list, bool const apply_check = true );

void
remove_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num );

bool
has_virtual_rna_residue_variant_type( core::pose::Pose & pose, core::Size const & seq_num );

void
remove_all_variant_types( core::pose::Pose & pose );

// Undefined, commenting out to fix PyRosetta build  core::Size get_matching_atom_name(std::string const & atom_name, core::conformation::Residue const & rsd);

void
setup_suite_atom_id_map( core::conformation::Residue const & rsd_1,
												core::conformation::Residue const & rsd_2,
												core::id::AtomID_Map < core::id::AtomID > & atom_ID_map,
												bool const base_only = true );


void
setup_suite_atom_id_map( core::conformation::Residue const & rsd_1,
                        core::conformation::Residue const & rsd_2,
                        core::Size const res_num_1,
											 core::Size const res_num_2, //allow for the possibility that two poses have different sizes Jun 9, 2010
                        core::id::AtomID_Map < core::id::AtomID > & atom_ID_map,
											 bool const base_only = true );

void
setup_suite_atom_id_map( core::pose::Pose const & pose_1,
											 core::pose::Pose const & pose_2,
											 core::Size const base_res,
											 core::id::AtomID_Map < core::id::AtomID > & atom_ID_map,
											 bool const base_only = true );


void
setup_suite_atom_id_map( core::pose::Pose const & pose_1,  ////Dec 23, 2011.
											 core::pose::Pose const & pose_2,
											 core::Size const base_res_1,
											 core::Size const base_res_2,
											 core::id::AtomID_Map < core::id::AtomID > & atom_ID_map,
											 bool const base_only = true );


core::id::AtomID_Map < core::id::AtomID >
create_alignment_id_map(	core::pose::Pose & mod_pose, core::pose::Pose const & ref_pose, utility::vector1< core::Size > const & rmsd_residue_list, bool const base_only = true );

void
align_poses( core::pose::Pose & moving_pose,
					  std::string const moving_tag,
						core::pose::Pose const & static_pose,
						std::string const static_tag,
						utility::vector1< core::Size > const & working_best_alignment,
						bool const base_only = true );

utility::vector1< core::Size >
apply_full_to_sub_mapping( utility::vector1< core::Size > const & res_vector, utility::vector1< core::Size > const & is_working_res, std::map< core::Size, core::Size > const & full_to_sub );

utility::vector1< core::Size >
apply_full_to_sub_mapping( utility::vector1< core::Size > const & res_vector, StepWiseRNA_JobParametersCOP job_parameters );


void
ensure_valid_full_seq_num( core::Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters );

bool
check_is_working_res( core::Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters );

core::Size
check_validity_and_get_working_res( core::Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters );


utility::vector1< core::Size >
apply_sub_to_full_mapping( utility::vector1< core::Size > const & working_res_vector, StepWiseRNA_JobParametersCOP job_parameters );

std::map< core::Size, core::Size >
create_full_to_input_res_map( utility::vector1< core::Size > const & input_res_vector );

/*
utility::vector1< core::Size >
apply_sub_to_full_mapping( utility::vector1< core::Size > const & working_res_vector, StepWiseRNA_JobParametersOP & job_parameters );
*/
/*
utility::vector1< core::Size >
apply_full_to_sub_mapping( utility::vector1< core::Size > const & res_vector, StepWiseRNA_JobParametersOP & job_parameters );
*/
/*

*/

core::Size
string_to_int( std::string const input_string );

core::Real
string_to_real( std::string const input_string );


utility::vector1< std::string >
Tokenize( std::string const str, std::string delimiters );

bool
Is_virtual_base( core::conformation::Residue const & rsd );

void
output_pair_size( std::pair < core::Size, core::Size > const & pair_size, std::ostream & outstream = std::cout );

void
output_pair_size( utility::vector1 < std::pair < core::Size, core::Size > > const & pair_size_vector, std::string const & output_string, std::ostream & outstream = std::cout, core::Size const spacing = 40 );

bool
pair_sort_citeria( std::pair < core::Size, core::Size > pair_one, std::pair < core::Size, core::Size > pair_two );

void sort_seq_num_list( utility::vector1< core::Size > & seq_num_list );

void Output_seq_num_list( std::string const tag, utility::vector1< core::Size > const & seq_num_list, std::ostream & outstream = std::cout, core::Size const spacing = 40 );

bool
Is_equivalent_vector( utility::vector1< core::Size > const & seq_num_list_1, utility::vector1< core::Size > const & seq_num_list_2 );

void Output_is_prepend_map( std::string const tag, std::map< core::Size, bool > const & my_map, core::Size const max_seq_num, std::ostream & outstream = std::cout, core::Size const tag_spacing = 40 );

void
Output_bool_list( std::string const tag, utility::vector1< bool > const & bool_list, std::ostream & outstream = std::cout, core::Size const spacing = 40 );

void
Output_bool_list( std::string const tag, utility::vector1< core::Size > const & size_list, std::ostream & outstream = std::cout, core::Size const spacing = 40 );

void
Output_size_list( std::string const tag, utility::vector1< core::Size > const & size_list, std::ostream & outstream = std::cout, core::Size const spacing = 40 );

// Undefined, commenting out to fix PyRosetta build  bool seq_num_list_sort_citeria(core::Size seq_num_1, Residue_info seq_num_2);

void
Sort_pair_list( utility::vector1< std::pair < core::Size, core::Size > > pair_list );

bool
Is_close_chain_break( core::pose::Pose const & pose );

//core::Size
//Get_five_prime_chain_break(core::pose::Pose const & pose);

void
Add_harmonic_chainbreak_constraint( core::pose::Pose & pose, core::Size const five_prime_res );

void
Output_fold_tree_info( core::kinematics::FoldTree const & fold_tree, std::string const pose_name, std::ostream & outstream = std::cout );

void
Output_fold_tree_info( core::pose::Pose const & pose, std::string pose_name, std::ostream & outstream = std::cout );

bool
file_exists( std::string const & file_name );

void
remove_file( std::string const & file_name );

void
output_rotamer( utility::vector1 < core::Real > & rotamer );

void
Add_virtual_O2Star_hydrogen( core::pose::Pose & pose );

bool
Remove_virtual_O2Star_hydrogen( core::pose::Pose & pose );

core::Real
suite_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2, core::Size const & seq_num, bool const prepend_res, bool const ignore_virtual_atom = false );

core::Real
full_length_rmsd_over_residue_list( core::pose::Pose const & pose1, core::pose::Pose const & pose2, utility::vector1 < core::Size > const & residue_list, std::string const & full_sequence, bool const verbose, bool const ignore_virtual_atom );


core::Real
rmsd_over_residue_list( core::pose::Pose const & pose1, core::pose::Pose const & pose2, utility::vector1 < core::Size > const & residue_list, std::map< core::Size, core::Size > const & full_to_sub, std::map< core::Size, bool > const & Is_prepend_map, bool const verbose, bool const ignore_virtual_atom = false );

core::Real
rmsd_over_residue_list( core::pose::Pose const & pose1, core::pose::Pose const & pose2, StepWiseRNA_JobParametersCOP job_parameters_, bool const ignore_virtual_atom = false );

void
Print_heavy_atoms( core::Size const & suite_num_1, core::Size const & suite_num_2, core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

core::Size
Get_num_side_chain_atom_from_res_name( core::chemical::AA const & res_aa, bool const verbose );


void
base_atoms_square_deviation( core::pose::Pose const & pose1, core::pose::Pose const & pose2, core::Size const & moving_res_1, core::Size const & moving_res_2, core::Size& atom_count, core::Real& sum_sd, bool verbose, bool const ignore_virtual_atom );

void
phosphate_square_deviation( core::pose::Pose const & pose1, core::pose::Pose const & pose2, core::Size const & moving_res_1, core::Size const & moving_res_2, core::Size& atom_count, core::Real& sum_sd, bool verbose, bool const ignore_virtual_atom );

core::Real
phosphate_base_phosphate_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2, core::Size const & moving_res_num, bool const ignore_virtual_atom );

void
phosphate_base_phosphate_square_deviation( core::pose::Pose const & pose1, core::pose::Pose const & pose2, core::Size const & moving_res_1, core::Size const & moving_res_2, core::Size & atom_count, core:: Real & sum_sd, bool verbose, bool const ignore_virtual_atom );

core::Real
atom_square_deviation( core::conformation::Residue const & rsd_1, core::conformation::Residue const & rsd_2, core::Size const & atomno_1, core::Size const & atomno_2, bool verbose );

void
suite_square_deviation( core::pose::Pose const & pose1, core::pose::Pose const & pose2, bool const & prepend_res, core::Size const & moving_res_1, core::Size const & moving_res_2, core::Size& atom_count, core::Real& sum_sd, bool verbose, bool const ignore_virtual_atom );

void
Output_title_text( std::string const title, std::ostream & outstream = std::cout );

bool
Check_chain_closable( numeric::xyzVector< core::Real > const & xyz_1, numeric::xyzVector< core::Real > const & xyz_2, core::Size const gap_size );

bool
Check_chain_closable( core::pose::Pose const & pose, core::Size const five_prime_chain_break_res, core::Size const gap_size );

bool
Check_chain_closable_floating_base( core::pose::Pose const & five_prime_pose,
                                   core::pose::Pose const & three_prime_pose,
																	 core::Size const five_prime_chain_break_res,
                                   core::Size const gap_size );

void
get_C4_C3_distance_range( core::conformation::Residue const & five_prime_rsd,
												 core::conformation::Residue const & three_prime_rsd,
												 core::Distance & C4_C3_dist_min,
												 core::Distance & C4_C3_dist_max );

void
Freeze_sugar_torsions( core::kinematics::MoveMap & mm, core::Size const total_residue );

void
Output_boolean( std::string const & tag, bool boolean, std::ostream & outstream = std::cout );

void
Output_boolean( bool boolean, std::ostream & outstream = std::cout );

void
Output_movemap( core::kinematics::MoveMap const & mm, core::pose::Pose const & pose, std::ostream & outstream = std::cout );

utility::vector1< core::Size >
get_surrounding_O2star_hydrogen( core::pose::Pose const & pose, utility::vector1< core::Size > const & moving_res, bool verbose = false );

void
o2star_minimize( core::pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn );

void
o2star_minimize( core::pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn, utility::vector1< core::Size > const & O2star_seq_num_list );

core::pack::task::PackerTaskOP
create_standard_o2star_pack_task( core::pose::Pose const & pose, utility::vector1< core::Size > const & O2star_pack_seq_num );

void
print_backbone_torsions( core::pose::Pose const & pose, core::Size five_prime_chainbreak );

void
Correctly_position_cutpoint_phosphate_torsions( core::pose::Pose & current_pose, core::Size const five_prime_chainbreak,  bool verbose = false );

void
copy_torsions_FROM_TO( core::id::TorsionID const start_torsion_ID, core::id::TorsionID const end_torsion_ID, core::pose::Pose const & template_pose, core::pose::Pose & pose );

core::Size
setup_chain_break_jump_point( core::pose::Pose & pose, core::Size const jump_point_one, core::Size const jump_point_two, core::Size const five_prime_cutpoint, bool const verbose );

void
remove_chain_break_jump_point( core::pose::Pose & pose, core::Size const five_prime_cutpoint, core::kinematics::FoldTree const fold_tree_without_cutpoint );

core::Size
setup_bulge_jump_point( core::pose::Pose & pose, core::Size const & moving_base, core::Size const & reference_base, bool verbose = false );

core::Size
make_cut_at_moving_suite( core::pose::Pose & pose, core::Size const & moving_suite );

core::Size
make_cut_at_moving_suite( core::kinematics::FoldTree & fold_tree, core::Size const & moving_suite );

// Undefined, commenting out to fix PyRosetta build
// void
// get_partition_definition( ObjexxFCL::FArray1D<bool> & partition_definition, core::kinematics::FoldTree const & fold_tree, core::Size const & moving_suite );

// Undefined, commenting out to fix PyRosetta build
// void
// get_partition_definition( ObjexxFCL::FArray1D<bool> & partition_definition, core::pose::Pose const & pose , core::Size const & moving_suite );

utility::vector1< bool >
get_partition_definition( core::pose::Pose const & pose, core::Size const & moving_suite );

void
apply_rotamer( core::pose::Pose & pose, utility::vector1< Torsion_Info >  const & rotamer_list );

bool
check_for_messed_up_structure( core::pose::Pose const & pose, std::string const & tag );

core::Size
Get_residue_base_state( core::pose::Pose const & pose, core::Size const seq_num );

core::Size
Get_residue_pucker_state( core::pose::Pose const & pose, core::Size const seq_num, bool verbose = false );

bool
Is_same_ribose_pucker( core::pose::Pose const & current_pose, core::pose::Pose const & cluster_center_pose, core::Size const seq_num );


void
sleep( core::Size mseconds );


void
setup_simple_fold_tree( core::pose::Pose & pose );


void
get_atom_coordinates(
    utility::vector1< std::pair < core::id::AtomID,
    numeric::xyzVector< core::Real > > > & xyz_list,
    core::Size const & seq_num,
    core::conformation::Residue const & rsd_at_origin,
    core::kinematics::Stub const & moving_res_base_stub );


void
import_pose_from_silent_file(
    core::pose::Pose & import_pose,
    std::string const & silent_file,
    std::string const & input_tag );


std::string
path_basename( std::string const full_path );

bool
Is_residues_in_contact(
    core::Size const & res_ONE,
    core::pose::Pose const & pose_ONE,
    core::Size const & res_TWO,
    core::pose::Pose const & pose_TWO,
    core::Real const atom_atom_overlap_dist_cutoff,
    core::Size const num_atom_contacts_cutoff,
    bool const verbose = false );

void
set_CCD_torsions_to_zero(
    core::pose::Pose & pose,
    core::Size const five_prime_res );

void
get_default_allowed_bulge_res(
    utility::vector1< core::Size > & allow_bulge_res_list,
    core::pose::Pose const & pose,
    bool const verbose );

core::Size
virtualize_bulges( core::pose::Pose & input_pose,
    utility::vector1< core::Size > const & in_allow_bulge_res_list,
    core::scoring::ScoreFunctionOP const & scorefxn,
    std::string const & tag,
    bool const allow_pre_virtualize,
    bool const allow_consecutive_bulges,
    bool const verbose );

std::string
get_tag_from_pdb_filename( std::string const pdb_filename );

void
move_jump_atom_to_base(
    core::kinematics::FoldTree & fold_tree,
    std::string const & working_sequence );

void
print_JobParameters_info( StepWiseRNA_JobParametersOP const & JP, std::string const JP_name, std::ostream & outstream = std::cout, bool const Is_simple_full_length_JP = false );

void
print_JobParameters_info( StepWiseRNA_JobParametersCOP const & const_JP, std::string const JP_name, std::ostream & outstream = std::cout, bool const Is_simple_full_length_JP = false );

//void
//add_native_base_pair_stats( core::io::silent::SilentStructOP & s, core::pose::Pose const & native_pose, utility::vector1 < core::Size > const & input_rmsd_res_list);

//void
//add_base_pair_stats( core::io::silent::SilentStructOP &s, core::pose::Pose const &  pose, core::pose::Pose const & native_pose, utility::vector1 < core::Size > const & input_rmsd_res_list);

void
set_nucleotide_to_A_form( core::pose::Pose & pose, core::Size const seq_num );

void
print_atom_info( core::pose::Pose const & pose, core::Size const seq_num, std::string const pose_name );

void
print_individual_atom_info( core::conformation::Residue const & rsd, core::Size const atomno, std::string const rsd_name );

void
print_base_state( std::string const tag, core::Size const base_state, std::ostream & outstream = std::cout );

void
print_ribose_pucker_state( std::string const tag, core::Size const pucker_state, std::ostream & outstream = std::cout );


void
initialize_common_scorefxns(
    core::scoring::ScoreFunctionOP const & scorefxn,
    core::scoring::ScoreFunctionOP & sampling_scorefxn,
    core::scoring::ScoreFunctionOP & atr_rep_screening_scorefxn,
    core::scoring::ScoreFunctionOP & chainbreak_scorefxn,
    core::scoring::ScoreFunctionOP & o2star_pack_scorefxn );


void
copy_all_o2star_torsions( core::pose::Pose & mod_pose, core::pose::Pose const & template_pose );

core::scoring::ScoreFunctionOP
rescale_scorefxn( core::scoring::ScoreFunctionOP const & starting_scorefxn, core::Real const scaling_factor );

void
show_scorefxn_weight_lines( core::scoring::ScoreFunctionOP const & scorefxn, std::string const title );
//Doesn't work on MAC!!
//void
//process_mem_usage(double& vm_usage, double& resident_set, core::Size count);

void
figure_out_swa_rna_movemap( core::kinematics::MoveMap & mm, core::pose::Pose const & pose, ObjexxFCL::FArray1D < bool > const & allow_insert );

void
figure_out_swa_rna_movemap( core::kinematics::MoveMap & mm, core::pose::Pose const & pose, utility::vector1< core::Size > const & minimize_res );

void
choose_random_if_unspecified_nucleotide( char & newrestype );

bool
mutate_res_if_allowed( core::pose::Pose & pose, core::Size const mutate_res, core::Real const mutation_frequency = 0.5 );

}
}
}

#endif
