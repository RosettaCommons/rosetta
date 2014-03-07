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


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_Util_hh
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_Util_hh


#include <protocols/stepwise/sampling/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
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
#include <protocols/toolbox/AllowInsert.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <string>
#include <map>

typedef  numeric::xyzMatrix< core::Real > Matrix;
using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {

bool is_OP2_atom( std::string const & atom_name );
bool is_OP1_atom( std::string const & atom_name );
bool is_P_atom( std::string const & atom_name );
bool is_O2prime_atom( std::string const & atom_name );
bool is_O3prime_atom( std::string const & atom_name );
bool is_O4prime_atom( std::string const & atom_name );
bool is_O5prime_atom( std::string const & atom_name );
bool is_C2prime_atom( std::string const & atom_name );
bool is_C3prime_atom( std::string const & atom_name );
bool is_C4prime_atom( std::string const & atom_name );
bool is_C5prime_atom( std::string const & atom_name );
bool is_1H5prime_atom( std::string const & atom_name );
bool is_2H5prime_atom( std::string const & atom_name );
bool is_H3prime_atom( std::string const & atom_name );
bool is_H4prime_atom( std::string const & atom_name );
bool is_three_prime_phosphate_atom( std::string const & atom_name );
bool is_five_prime_phosphate_atom( std::string const & atom_name );
bool is_phosphate_atom( std::string const & atom_name );


void
minimize_with_constraints( pose::Pose & pose, kinematics::MoveMap const & mm, scoring::ScoreFunctionOP const & scorefxn, optimization::MinimizerOptions const & options );

bool
check_can_prepend( utility::vector1< Size > const & seq_num_list );

bool
check_can_append( utility::vector1< Size > const & seq_num_list );

void
apply_protonated_H1_adenosine_variant_type( pose::Pose & pose, Size const & seq_num, bool const apply_check = true );

void
apply_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num, bool const apply_check = true );

void
apply_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num, utility::vector1< Size > const & working_cutpoint_closed_list, bool const apply_check = true );

void
remove_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num );

bool
has_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num );

void
remove_all_variant_types( pose::Pose & pose );

// Undefined, commenting out to fix PyRosetta build  Size get_matching_atom_name(std::string const & atom_name, conformation::Residue const & rsd);

void
setup_suite_atom_id_map( conformation::Residue const & rsd_1,
												conformation::Residue const & rsd_2,
												id::AtomID_Map < id::AtomID > & atom_ID_map,
												bool const base_only = true );


void
setup_suite_atom_id_map( conformation::Residue const & rsd_1,
                        conformation::Residue const & rsd_2,
                        Size const res_num_1,
											 Size const res_num_2, //allow for the possibility that two poses have different sizes Jun 9, 2010
                        id::AtomID_Map < id::AtomID > & atom_ID_map,
											 bool const base_only = true );

void
setup_suite_atom_id_map( pose::Pose const & pose_1,
											 pose::Pose const & pose_2,
											 Size const base_res,
											 id::AtomID_Map < id::AtomID > & atom_ID_map,
											 bool const base_only = true );


void
setup_suite_atom_id_map( pose::Pose const & pose_1,  ////Dec 23, 2011.
											 pose::Pose const & pose_2,
											 Size const base_res_1,
											 Size const base_res_2,
											 id::AtomID_Map < id::AtomID > & atom_ID_map,
											 bool const base_only = true );


id::AtomID_Map < id::AtomID >
create_alignment_id_map(	pose::Pose & mod_pose, pose::Pose const & ref_pose, utility::vector1< Size > const & rmsd_residue_list, bool const base_only = true );

void
align_poses( pose::Pose & moving_pose,
					  std::string const moving_tag,
						pose::Pose const & static_pose,
						std::string const static_tag,
						utility::vector1< Size > const & working_best_alignment,
						bool const base_only = true );

utility::vector1< Size >
apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, utility::vector1< Size > const & is_working_res, std::map< Size, Size > const & full_to_sub );

utility::vector1< Size >
apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, StepWiseRNA_JobParametersCOP job_parameters );


void
ensure_valid_full_seq_num( Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters );

bool
check_is_working_res( Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters );

Size
check_validity_and_get_working_res( Size const full_seq_num, StepWiseRNA_JobParametersCOP const & job_parameters );


utility::vector1< Size >
apply_sub_to_full_mapping( utility::vector1< Size > const & working_res_vector, StepWiseRNA_JobParametersCOP job_parameters );

std::map< Size, Size >
create_full_to_input_res_map( utility::vector1< Size > const & input_res_vector );

Size
string_to_int( std::string const input_string );

Real
string_to_real( std::string const input_string );


utility::vector1< std::string >
tokenize( std::string const str, std::string delimiters );

bool
is_virtual_base( conformation::Residue const & rsd );

void
output_pair_size( std::pair < Size, Size > const & pair_size, std::ostream & outstream = std::cout );

void
output_pair_size( utility::vector1 < std::pair < Size, Size > > const & pair_size_vector, std::string const & output_string, std::ostream & outstream = std::cout, Size const spacing = 40 );

bool
pair_sort_criterion( std::pair < Size, Size > pair_one, std::pair < Size, Size > pair_two );

void sort_seq_num_list( utility::vector1< Size > & seq_num_list );

void output_seq_num_list( std::string const tag, utility::vector1< Size > const & seq_num_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

bool
is_equivalent_vector( utility::vector1< Size > const & seq_num_list_1, utility::vector1< Size > const & seq_num_list_2 );

void output_is_prepend_map( std::string const tag, std::map< Size, bool > const & my_map, Size const max_seq_num, std::ostream & outstream = std::cout, Size const tag_spacing = 40 );

void
output_bool_list( std::string const tag, utility::vector1< bool > const & bool_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

void
output_bool_list( std::string const tag, utility::vector1< Size > const & size_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

void
output_size_list( std::string const tag, utility::vector1< Size > const & size_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

// Undefined, commenting out to fix PyRosetta build  bool seq_num_list_sort_criterion(Size seq_num_1, Residue_info seq_num_2);

void
sort_pair_list( utility::vector1< std::pair < Size, Size > > pair_list );

void
output_fold_tree_info( kinematics::FoldTree const & fold_tree, std::string const pose_name, std::ostream & outstream = std::cout );

void
output_fold_tree_info( pose::Pose const & pose, std::string pose_name, std::ostream & outstream = std::cout );

bool
file_exists( std::string const & file_name );

void
remove_file( std::string const & file_name );

void
output_rotamer( utility::vector1 < Real > & rotamer );

void
add_virtual_O2Prime_hydrogen( pose::Pose & pose );

bool
remove_virtual_O2Prime_hydrogen( pose::Pose & pose );

Real
suite_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & seq_num, bool const prepend_res, bool const ignore_virtual_atom = false );

Real
full_length_rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::string const & full_sequence, bool const verbose, bool const ignore_virtual_atom );


Real
rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::map< Size, Size > const & full_to_sub, std::map< Size, bool > const & is_prepend_map, bool const verbose, bool const ignore_virtual_atom = false );

Real
rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, StepWiseRNA_JobParametersCOP job_parameters_, bool const ignore_virtual_atom = false );

void
print_heavy_atoms( Size const & suite_num_1, Size const & suite_num_2, pose::Pose const & pose1, pose::Pose const & pose2 );

Size
get_num_side_chain_atom_from_res_name( chemical::AA const & res_aa, bool const verbose );


void
base_atoms_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom );

void
phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom );

Real
phosphate_base_phosphate_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_num, bool const ignore_virtual_atom );

void
phosphate_base_phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size & atom_count,  Real & sum_sd, bool verbose, bool const ignore_virtual_atom );

Real
atom_square_deviation( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const & atomno_1, Size const & atomno_2, bool verbose );

void
suite_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, bool const & prepend_res, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom );

void
output_title_text( std::string const title, std::ostream & outstream = std::cout );


void
Freeze_sugar_torsions( kinematics::MoveMap & mm, Size const total_residue );

void
output_boolean( std::string const & tag, bool boolean, std::ostream & outstream = std::cout );

void
output_boolean( bool boolean, std::ostream & outstream = std::cout );

void
output_movemap( kinematics::MoveMap const & mm, pose::Pose const & pose, std::ostream & outstream = std::cout );

utility::vector1< Size >
get_surrounding_O2prime_hydrogen( pose::Pose const & pose, utility::vector1< Size > const & moving_res, bool verbose = false );

void
o2prime_trials( pose::Pose& pose, scoring::ScoreFunctionOP const & packer_scorefxn );

void
o2prime_trials( pose::Pose& pose, scoring::ScoreFunctionOP const & packer_scorefxn, utility::vector1< Size > const & O2prime_seq_num_list );

pack::task::PackerTaskOP
create_standard_o2prime_pack_task( pose::Pose const & pose, utility::vector1< Size > const & O2prime_pack_seq_num );

void
print_backbone_torsions( pose::Pose const & pose, Size five_prime_chainbreak );

void
copy_torsions_FROM_TO( id::TorsionID const start_torsion_ID, id::TorsionID const end_torsion_ID, pose::Pose const & template_pose, pose::Pose & pose );

core::Size
setup_chain_break_jump_point( core::pose::Pose & pose,
															core::Size const moving_res,
															core::Size const reference_res );

void
remove_chain_break_jump_point( core::pose::Pose & pose,
															 core::Size const moving_res,
															 core::Size const reference_res );

void
setup_chain_break_variants( core::pose::Pose & pose,  Size const cutpoint );

void
remove_chain_break_variants( core::pose::Pose & pose,  Size const & cutpoint );

Size
setup_bulge_jump_point( pose::Pose & pose, Size const & moving_base, Size const & reference_base, bool verbose = false );


utility::vector1< bool >
get_partition_definition_floating_base( pose::Pose const & pose, Size const & moving_res );

Size
get_anchor_res( Size const rebuild_res, pose::Pose const & pose );

void
apply_rotamer( pose::Pose & pose, utility::vector1< Torsion_Info >  const & rotamer_list );

bool
check_for_messed_up_structure( pose::Pose const & pose, std::string const & tag );

Size
get_residue_base_state( pose::Pose const & pose, Size const seq_num );

Size
get_residue_pucker_state( pose::Pose const & pose, Size const seq_num, bool verbose = false );

bool
is_same_sugar_pucker( pose::Pose const & current_pose, pose::Pose const & cluster_center_pose, Size const seq_num );

void
sleep( Size mseconds );

void
setup_simple_fold_tree( pose::Pose & pose );

void
import_pose_from_silent_file(
    pose::Pose & import_pose,
    std::string const & silent_file,
    std::string const & input_tag );

std::string
path_basename( std::string const full_path );

bool
is_residues_in_contact(
    Size const & res_ONE,
    pose::Pose const & pose_ONE,
    Size const & res_TWO,
    pose::Pose const & pose_TWO,
    Real const atom_atom_overlap_dist_cutoff,
    Size const num_atom_contacts_cutoff,
    bool const verbose = false );

void
set_CCD_torsions_to_zero(
    pose::Pose & pose,
    Size const five_prime_res );

void
get_default_allowed_bulge_res(
    utility::vector1< Size > & allow_bulge_res_list,
    pose::Pose const & pose,
    bool const verbose );

Size
virtualize_bulges( pose::Pose & input_pose,
    utility::vector1< Size > const & in_allow_bulge_res_list,
    scoring::ScoreFunctionOP const & scorefxn,
    std::string const & tag,
    bool const allow_pre_virtualize,
    bool const allow_consecutive_bulges,
    bool const verbose );

std::string
get_tag_from_pdb_filename( std::string const pdb_filename );

void
move_jump_atom_to_base(
    kinematics::FoldTree & fold_tree,
    std::string const & working_sequence );

void
print_JobParameters_info( StepWiseRNA_JobParametersOP const & JP, std::string const JP_name, std::ostream & outstream = std::cout, bool const is_simple_full_length_JP = false );

void
print_JobParameters_info( StepWiseRNA_JobParametersCOP const & const_JP, std::string const JP_name, std::ostream & outstream = std::cout, bool const is_simple_full_length_JP = false );

void
set_nucleotide_to_A_form( pose::Pose & pose, Size const seq_num );

void
print_atom_info( pose::Pose const & pose, Size const seq_num, std::string const pose_name );

void
print_individual_atom_info( conformation::Residue const & rsd, Size const atomno, std::string const rsd_name );

void
print_base_state( std::string const tag, Size const base_state, std::ostream & outstream = std::cout );

void
print_sugar_pucker_state( std::string const tag, Size const pucker_state, std::ostream & outstream = std::cout );


scoring::ScoreFunctionOP
get_sampling_scorefxn( scoring::ScoreFunctionCOP scorefxn_ );

void
initialize_common_scorefxns(
    scoring::ScoreFunctionOP const & scorefxn,
    scoring::ScoreFunctionOP & sampling_scorefxn,
    scoring::ScoreFunctionOP & atr_rep_screening_scorefxn,
    scoring::ScoreFunctionOP & chainbreak_scorefxn,
    scoring::ScoreFunctionOP & o2prime_pack_scorefxn );


void
copy_all_o2prime_torsions( pose::Pose & mod_pose, pose::Pose const & template_pose );

scoring::ScoreFunctionOP
rescale_scorefxn( scoring::ScoreFunctionOP const & starting_scorefxn, Real const scaling_factor );

void
show_scorefxn_weight_lines( scoring::ScoreFunctionOP const & scorefxn, std::string const title );

void
figure_out_stepwise_rna_movemap( kinematics::MoveMap & mm, pose::Pose const & pose, utility::vector1< Size > const & minimize_res );

void
figure_out_stepwise_rna_movemap( core::kinematics::MoveMap & mm, core::pose::Pose const & pose, toolbox::AllowInsertOP const & allow_insert );

void
update_allow_insert_with_extra_minimize_res( pose::Pose const & pose, toolbox::AllowInsertOP & allow_insert, utility::vector1< core::Size > const & extra_minimize_res );

void
choose_random_if_unspecified_nucleotide( char & newrestype );

bool
mutate_res_if_allowed( pose::Pose & pose, Size const mutate_res, Real const mutation_frequency = 0.5 );

std::string
create_tag( std::string const & prestring, Size const i );

std::string //silly function to convert to real to string
create_torsion_value_string( core::Real const & torsion_value );

std::string //silly function used for appending the rotamer value to the tag
create_rotamer_string( core::pose::Pose const & pose, Size const moving_res, bool const is_prepend );

std::string //silly function used for appending the rotamer value to the tag
create_rotamer_string( core::pose::Pose const & pose, Size const moving_res, Size const reference_res );

void
add_fade_chain_break_constraint_across_gap( pose::Pose & pose,
																						Size const five_prime_res,
																						Size const three_prime_res,
																						Size const gap_size );

void
add_harmonic_chain_break_constraint( pose::Pose & pose, Size const five_prime_res );

void
get_possible_O3prime_C5prime_distance_range( Size const gap_size_, Distance & min_dist, Distance & max_dist );


} //rna
} //sampling
} //stepwise
} //protocols

#endif
