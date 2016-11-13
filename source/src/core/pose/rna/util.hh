// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/util.hh
/// @brief
/// @author Rhiju

#ifndef INCLUDED_core_pose_rna_util_hh
#define INCLUDED_core_pose_rna_util_hh

#include <core/types.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/types.hh>
#include <utility/vector1.fwd.hh>
#include <core/pose/rna/VDW_Grid.hh>

namespace core {
namespace pose {
namespace rna {

bool
mutate_position( core::pose::Pose & pose, core::Size const i, char const & new_seq );

bool
mutate_position( pose::Pose & pose, Size const i, std::string const & name3 );

bool
mutate_position( pose::Pose & pose, Size const i, core::chemical::ResidueType const & rt );

void
figure_out_reasonable_rna_fold_tree( core::pose::Pose & pose,
																		 bool const force_cut_at_rna_chainbreak = false );

void
virtualize_5prime_phosphates( core::pose::Pose & pose );

bool
is_cutpoint_open( Pose const & pose, Size const i );

bool
is_rna_chainbreak( Pose const & pose, Size const i );

void
fix_sugar_coords_WORKS_BUT_SLOW(
	utility::vector1< std::string> atoms_for_which_we_need_new_dofs,
	utility::vector1< utility::vector1< id::DOF_Type > > which_dofs,
	utility::vector1< Vector > const & non_main_chain_sugar_coords,
	Pose & pose,
	Size const i);

void
prepare_scratch_residue(
	conformation::ResidueOP & scratch_rsd,
	conformation::Residue const & start_rsd,
	utility::vector1< Vector > const & non_main_chain_sugar_coords,
	Pose const & pose);

void
fix_sugar_coords(
	utility::vector1< std::string> atoms_for_which_we_need_new_dofs,
	utility::vector1< Vector > const & non_main_chain_sugar_coords,
	Pose & pose,
	Pose const & reference_pose,
	Size const i);

void
initialize_atoms_for_which_we_need_new_dofs(
	utility::vector1< std::string > & atoms_for_which_we_need_new_dofs,
	Pose const & pose,
	Size const i );

void
apply_non_main_chain_sugar_coords(
	utility::vector1< Vector > const & non_main_chain_sugar_coords,
	Pose & pose,
	Pose const & reference_pose,
	Size const i);

void
apply_ideal_c2endo_sugar_coords(
	Pose & pose,
	Size const i);

core::chemical::rna::PuckerState
assign_pucker(
	Pose const & pose,
	Size const rsd_id );

void
apply_pucker(
	Pose & pose,
	Size const i,
	core::chemical::rna::PuckerState const pucker_state = core::chemical::rna::ANY_PUCKER,
	bool const skip_same_state = false,
	bool const idealize_coord = true );

// DEPRECATED
void
apply_ideal_c2endo_sugar_coords(
	pose::Pose & pose,
	Size const i );

void
position_cutpoint_phosphate_torsions( pose::Pose & current_pose,
																			Size const five_prime_chainbreak,
																			Size three_prime_chainbreak = 0 );

bool is_cutpoint_closed_torsion( pose::Pose const & pose, id::TorsionID const & torsion_id );

bool is_cutpoint_closed_atom( core::chemical::ResidueType const & rsd, id::AtomID const & id );

void output_boolean( std::string const & tag, bool boolean );

void print_torsion_info( pose::Pose const & pose, id::TorsionID const & torsion_id );

bool is_torsion_valid( pose::Pose const & pose, id::TorsionID const & torsion_id,
	bool verbose = false, bool skip_chainbreak_torsions = false );


// might be deprecatable soon -- legacy of VirtualSugarSampler, which could get away with bulge variant.
void
apply_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num, bool const apply_check = true );

// might be deprecatable soon -- legacy of VirtualSugarSampler, which could get away with bulge variant.
void
apply_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num, utility::vector1< Size > const & working_cutpoint_closed_list, bool const apply_check = true );

// might be deprecatable soon -- legacy of VirtualSugarSampler, which could get away with bulge variant.
void
remove_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num );

// might be deprecatable soon -- legacy of VirtualSugarSampler, which could get away with bulge variant.
bool
has_virtual_rna_residue_variant_type( pose::Pose & pose, Size const & seq_num );

void
apply_Aform_torsions( pose::Pose & pose, Size const n );

core::chemical::rna::ChiState
get_residue_base_state( core::pose::Pose const & pose, Size const seq_num );

core::chemical::rna::PuckerState
get_residue_pucker_state( core::pose::Pose const & pose, Size const seq_num );

Real
get_op2_op1_sign( pose::Pose const & pose );

Real
get_op2_op1_sign( pose::Pose const & pose , Size res_num);

void
make_phosphate_nomenclature_matches_mini( pose::Pose & pose);

void
add_virtual_O2Prime_hydrogen( pose::Pose & pose );

Atom_Bin
get_atom_bin( numeric::xyzVector< core::Real > const & atom_pos, numeric::xyzVector< core::Real > const & ref_xyz,
	core::Real const atom_bin_size, int const bin_offset );

bool
is_atom_bin_in_range( Atom_Bin const & atom_pos_bin, int const bin_max );

utility::vector1< std::string >
tokenize( std::string const & str, std::string const & delimiters );

core::Size
string_to_int( std::string const & input_string );

/// @brief 'suite' backbone torsion -- useful in setting up cutpoint with OVL1, OVL2, OVU atoms
utility::vector1< std::pair< id::TorsionID, Real > >
get_suite_torsion_info( core::pose::Pose const & pose, Size const i );

/// @brief useful in setting up cutpoint with OVL1, OVL2, OVU atoms
void
apply_suite_torsion_info( core::pose::Pose & pose,
	utility::vector1< std::pair< id::TorsionID, Real > > const & suite_torsion_info );

std::string
remove_bracketed( std::string const & sequence );

void
remove_and_store_bracketed(
	std::string const & working_sequence,
	std::string & working_sequence_clean,
	std::map< Size, std::string > & special_res );

} //ns rna
} //ns pose
} //ns core

#endif
