// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
/// @brief
/// @details
///
/// @author Rhiju Das

#ifndef INCLUDED_protocols_stepwise_rna_util_hh
#define INCLUDED_protocols_stepwise_rna_util_hh

#include <protocols/stepwise/modeler/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh>
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
using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {

void
minimize_with_constraints( pose::Pose & pose, kinematics::MoveMap const & mm, scoring::ScoreFunctionOP const & scorefxn, optimization::MinimizerOptions const & options );

void
apply_protonated_H1_adenosine_variant_type( pose::Pose & pose, Size const & seq_num, bool const apply_check = true );

void
remove_all_variant_types( pose::Pose & pose );

// may be deprecatable soon -- in use by superimpose_pose, which may be replaced by StepWisePoseAligner,
//  and and in use by CombineLongLoopFilterer, which again could be fixed soon.
void
setup_suite_atom_id_map( conformation::Residue const & rsd_1,
	conformation::Residue const & rsd_2,
	id::AtomID_Map < id::AtomID > & atom_ID_map,
	bool const base_only = true );

// may be deprecatable soon -- in use by superimpose_pose, which may be replaced by StepWisePoseAligner,
//  and and in use by CombineLongLoopFilterer, which again could be fixed soon.
void
setup_suite_atom_id_map( conformation::Residue const & rsd_1,
	conformation::Residue const & rsd_2,
	Size const res_num_1,
	Size const res_num_2, //allow for the possibility that two poses have different sizes Jun 9, 2010
	id::AtomID_Map < id::AtomID > & atom_ID_map,
	bool const base_only = true );

// may be deprecatable soon -- in use by superimpose_pose, which may be replaced by StepWisePoseAligner,
//  and and in use by CombineLongLoopFilterer, which again could be fixed soon.
void
setup_suite_atom_id_map( pose::Pose const & pose_1,
	pose::Pose const & pose_2,
	Size const base_res,
	id::AtomID_Map < id::AtomID > & atom_ID_map,
	bool const base_only = true );

// may be deprecatable soon -- in use by superimpose_pose, which may be replaced by StepWisePoseAligner,
//  and and in use by CombineLongLoopFilterer, which again could be fixed soon.
void
setup_suite_atom_id_map( pose::Pose const & pose_1,  ////Dec 23, 2011.
	pose::Pose const & pose_2,
	Size const base_res_1,
	Size const base_res_2,
	id::AtomID_Map < id::AtomID > & atom_ID_map,
	bool const base_only = true );


// may be deprecatable soon -- in use by superimpose_pose, which may be replaced by StepWisePoseAligner,
void
align_poses( pose::Pose & moving_pose,
	std::string const & moving_tag,
	pose::Pose const & static_pose,
	std::string const & static_tag,
	utility::vector1< Size > const & working_best_alignment,
	bool const base_only = true );

utility::vector1< Size >
apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, utility::vector1< Size > const & is_working_res, std::map< Size, Size > const & full_to_sub );

utility::vector1< Size >
apply_full_to_sub_mapping( utility::vector1< Size > const & res_vector, working_parameters::StepWiseWorkingParametersCOP working_parameters );

void
ensure_valid_full_seq_num( Size const full_seq_num, working_parameters::StepWiseWorkingParametersCOP const & working_parameters );

bool
check_is_working_res( Size const full_seq_num, working_parameters::StepWiseWorkingParametersCOP const & working_parameters );

Size
check_validity_and_get_working_res( Size const full_seq_num, working_parameters::StepWiseWorkingParametersCOP const & working_parameters );

utility::vector1< Size >
apply_sub_to_full_mapping( utility::vector1< Size > const & working_res_vector, working_parameters::StepWiseWorkingParametersCOP working_parameters );

std::map< Size, Size >
create_full_to_input_res_map( utility::vector1< Size > const & input_res_vector );

Size
string_to_int( std::string const & input_string );

Real
string_to_real( std::string const & input_string );

utility::vector1< std::string >
tokenize( std::string const str, std::string delimiters );

bool
is_virtual_base( conformation::Residue const & rsd );

void sort_seq_num_list( utility::vector1< Size > & seq_num_list );

void output_seq_num_list( std::string const & tag, utility::vector1< Size > const & seq_num_list, std::ostream & outstream = std::cout, Size const spacing = 40 );

bool
file_exists( std::string const & file_name );

void
remove_file( std::string const & file_name );

//void
//output_rotamer( utility::vector1 < Real > & rotamer );

void
add_virtual_O2Prime_hydrogen( pose::Pose & pose );

void
check_instantiated_O2Prime_hydrogen( core::pose::Pose const & pose );

bool
remove_virtual_O2Prime_hydrogen( pose::Pose & pose );

Real
suite_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & seq_num, bool const prepend_res, bool const ignore_virtual_atom = false );

Real
rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::map< Size, Size > const & full_to_sub, std::map< Size, bool > const & is_prepend_map, bool const verbose, bool const ignore_virtual_atom = false );

Real
rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, working_parameters::StepWiseWorkingParametersCOP working_parameters_, bool const ignore_virtual_atom = false );

//void
//print_heavy_atoms( Size const & suite_num_1, Size const & suite_num_2, pose::Pose const & pose1, pose::Pose const & pose2 );

//Size
//get_num_side_chain_atom_from_res_name( chemical::AA const & res_aa, bool const verbose );


Real
phosphate_base_phosphate_rmsd( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_num, bool const ignore_virtual_atom );

void
phosphate_base_phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size & atom_count,  Real & sum_sd, bool verbose, bool const ignore_virtual_atom );

void
base_atoms_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom );


void
phosphate_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom );


Real
atom_square_deviation( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const & atomno_1, Size const & atomno_2, bool verbose );

void
suite_square_deviation( pose::Pose const & pose1, pose::Pose const & pose2, bool const & prepend_res, Size const & moving_res_1, Size const & moving_res_2, Size& atom_count, Real& sum_sd, bool verbose, bool const ignore_virtual_atom );


void
output_title_text( std::string const & title, std::ostream & outstream = std::cout );

void
freeze_sugar_torsions( kinematics::MoveMap & mm, Size const total_residue );

utility::vector1< Size >
get_surrounding_O2prime_hydrogen( pose::Pose const & pose, utility::vector1< Size > const & moving_res, bool verbose = false );

void
o2prime_trials( pose::Pose& pose, scoring::ScoreFunctionCOP const & packer_scorefxn,
	bool const pack_virtual_o2prime_hydrogen = false );

void
o2prime_trials( pose::Pose& pose, scoring::ScoreFunctionCOP const & packer_scorefxn, utility::vector1< Size > const & O2prime_seq_num_list,
	bool const pack_virtual_o2prime_hydrogen = false );

pack::task::PackerTaskOP
create_standard_o2prime_pack_task( pose::Pose const & pose, utility::vector1< Size > const & O2prime_pack_seq_num,
	bool const pack_virtual_o2prime_hydrogen = false );

void
setup_chain_break_variants( core::pose::Pose & pose,  Size const cutpoint );

void
remove_chain_break_variants( core::pose::Pose & pose,  Size const & cutpoint );

utility::vector1< bool >
get_partition_definition_floating_base( pose::Pose const & pose, Size const & moving_res );

Size
get_anchor_res( Size const rebuild_res, pose::Pose const & pose );

bool
check_for_messed_up_structure( pose::Pose const & pose, std::string const & tag );

//void
//sleep( Size mseconds );

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
print_base_state( std::string const & tag, Size const base_state, std::ostream & outstream = std::cout );

void
print_sugar_pucker_state( std::string const & tag, Size const pucker_state, std::ostream & outstream = std::cout );

scoring::ScoreFunctionOP
get_modeler_scorefxn( scoring::ScoreFunctionCOP scorefxn_ );

void
copy_all_o2prime_torsions( pose::Pose & mod_pose, pose::Pose const & template_pose );

scoring::ScoreFunctionOP
rescale_scorefxn( scoring::ScoreFunctionOP const & starting_scorefxn, Real const scaling_factor );

void
show_scorefxn_weight_lines( scoring::ScoreFunctionOP const & scorefxn, std::string const & title );

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

void
remove_all_virtual_phosphates( core::pose::Pose & pose );

utility::vector1< core::Size >
just_rna( utility::vector1< core::Size > const & res_list, core::pose::Pose const & pose );

void
figure_out_moving_rna_chain_breaks( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_partition_res,
	utility::vector1< core::Size > & rna_cutpoints_closed,
	utility::vector1< core::Size > & rna_five_prime_chain_breaks,
	utility::vector1< core::Size > & rna_three_prime_chain_breaks,
	utility::vector1< core::Size > & rna_chain_break_gap_sizes );
void
virtualize_free_rna_moieties( pose::Pose & pose );

bool
just_modeling_RNA( std::string const & sequence );


} //rna
} //modeler
} //stepwise
} //protocols

#endif
