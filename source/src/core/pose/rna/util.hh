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
/// @author Rhiju Das

#ifndef INCLUDED_core_pose_rna_util_hh
#define INCLUDED_core_pose_rna_util_hh

#include <core/types.hh>
#include <core/chemical/rna/util.hh> // for some enums declared here, no fwd
#include <core/pose/PDBInfo.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/scoring/rna/RNA_Motif.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/BasePair.fwd.hh>
#include <core/pose/rna/BasePairStep.fwd.hh>
#include <core/pose/rna/StubStubType.fwd.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/rna/BaseStack.fwd.hh>
#include <core/pose/rna/StubStubType.fwd.hh>
#include <core/pose/rna/VDW_Grid.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <tuple>
#include <map>

namespace core {
namespace pose {
namespace rna {

bool
mutate_position( core::pose::Pose & pose, core::Size const i, char const new_seq );

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

bool
check_in_base_pair_list( core::pose::rna::BasePair const & base_pair /*from native*/,
	utility::vector1< core::pose::rna::BasePair > const & base_pair_list /*for decoy*/);

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
get_number_base_pairs(pose::Pose const & pose,
	Size & N_WC,
	Size & N_NWC,
	Size & N_BP
);

void
add_number_base_pairs( pose::Pose const & pose, io::silent::SilentStruct & s );

void
add_number_base_pairs( pose::Pose & pose );

void
get_number_native_base_pairs(pose::Pose & pose, pose::Pose const & native_pose,
	Size & pN_WC,
	Size & pN_NWC,
	Size & pN_BP,
	Size & pnatWC,
	Size & pnatNWC,
	Size & pnatBP,
	Real & pf_natWC,
	Real & pf_natNWC,
	Real & pf_natBP
);

void
add_number_native_base_pairs(pose::Pose & pose, pose::Pose const & native_pose, io::silent::SilentStruct & s );

void
add_number_native_base_pairs(pose::Pose & pose, pose::Pose const & native_pose );

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

void
add_chi_constraints( pose::Pose & pose,
	core::scoring::func::FuncOP chi_potential_restraint,
	utility::vector1< Size > const & rna_chi_res );

void
add_syn_chi_constraints( core::pose::Pose & pose );

void
add_anti_chi_constraints( core::pose::Pose & pose );

void
add_syn_anti_chi_constraints( core::pose::Pose & pose );

utility::vector1< core::id::TorsionID >
get_suite_torsion_ids( Size const i );

void
get_stub_stub( core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::kinematics::Stub & stub1,
	core::kinematics::Stub & stub2,
	StubStubType const & stub_stub_type );

void
setup_base_pair_constraints(
	core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size > > const &  pairings,
	core::Real const scale_factor = 1.0,
	bool const use_flat_harmonic = false );

void
get_base_pairing_list(
	core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size> > & base_pairing_list );

void
assert_phosphate_nomenclature_matches_mini( core::pose::Pose const & pose);

void
set_output_res_and_chain( core::pose::Pose & extended_pose,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & output_resnum_and_chain_and_segid );

void
virtualize_free_rna_moieties( core::pose::Pose & pose );

utility::vector1< bool >
detect_base_contacts( core::pose::Pose const & pose );

bool
check_phosphate_contacts_donor( utility::vector1< core::Vector > const & op_xyz_list,
	utility::vector1< core::Vector > const & donor_atom_xyz_list,
	utility::vector1< core::Vector > const & donor_base_atom_xyz_list );

bool
check_phosphate_contacts_donor( core::pose::Pose const & pose, core::Size const n );

utility::vector1< bool >
detect_phosphate_contacts( core::pose::Pose const & pose );


utility::vector1< bool >
detect_sugar_contacts( core::pose::Pose const & pose );

bool
detect_sugar_contacts( core::pose::Pose const & pose, core::Size const moving_res,
	core::Distance const o2prime_contact_distance_cutoff_ = 3.2 /*hydrogen bond*/ );


void
setup_three_prime_phosphate_based_on_next_residue( core::pose::Pose & pose, core::Size const n );


enum Terminus {
	FIVE_PRIME,
	THREE_PRIME
};

void
get_phosphate_atom_and_neighbor_list( core::pose::Pose const & pose,
	Size const n,
	Terminus const t,
	utility::vector1< core::Vector > & donor_atom_xyz_list,
	utility::vector1< core::Vector > & donor_base_atom_xyz_list,
	utility::vector1< core::Size > & neighbor_copy_dofs );




bool
moveable_jump( core::id::AtomID const & jump_atom_id1,
	core::id::AtomID const & jump_atom_id2,
	core::pose::toolbox::AtomLevelDomainMap const & atom_level_domain_map);

bool
moveable_jump( core::Size const jump_pos1,
	core::Size const jump_pos2,
	core::pose::toolbox::AtomLevelDomainMap const & atom_level_domain_map);

bool
base_pair_step_moving( core::pose::rna::BasePairStep const & base_pair_step,
	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::pose::Pose const & pose );

bool
base_pair_moving( core::pose::rna::BasePair const & base_pair,
	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::pose::Pose const & pose );




utility::vector1< core::Size >
get_rigid_body_jumps( core::pose::Pose const & pose );

void
fill_in_default_jump_atoms( core::kinematics::FoldTree & f, core::pose::Pose const & pose );

void
fill_in_default_jump_atoms( core::pose::Pose & pose );



Vector
get_sugar_centroid( core::conformation::Residue const & rsd );


void
figure_out_secstruct( core::pose::Pose & pose );

void
get_base_pairing_info( core::pose::Pose const & pose,
	core::Size const & seqpos,
	char & secstruct,
	ObjexxFCL::FArray1D <bool> & edge_is_base_pairing );

/// @brief Output base pairs detected for RNA, including noncanonicals. Must previously score pose with RNA_LORES.
void
output_base_pairs( std::ostream & out, core::pose::rna::RNA_BasePairList const & base_pair_list, core::pose::Pose const & pose  );

/// @brief Output base stacks detected for RNA. Must previously score pose with RNA_LORES.
void
output_base_stacks( std::ostream & out, core::pose::rna::RNA_BaseStackList const & base_stack_list, core::pose::Pose const & pose  );

/// @brief Output stems (>=2 base-pair helices) detected for RNA, including noncanonicals. Must previously score pose with RNA_LORES.
void
output_stems( std::ostream & out, core::scoring::rna::RNA_Motifs const & rna_motifs, core::pose::Pose const & pose );

/// @brief Figure out chains that have RNA in them. Can supply chains from command-line to focus on particular RNA chains.
utility::vector1< std::pair< char, std::string > >
figure_out_rna_chains( pose::Pose const & pose, utility::vector1< std::string > const & chains = utility::vector1< std::string >() );

/// @brief Pull out RNA from pose, in chains specified by chain_segids.
pose::Pose
extract_rna_chains( pose::Pose const & full_pose, utility::vector1< pose::ChainSegID > const & chain_segids );

/// @brief Output contacts of RNA chains (specified in chain_segids) to any non-RNA chains ("ligands").
void
output_ligands( std::ostream & out, pose::Pose const & pose,  utility::vector1< pose::ChainSegID > const & chain_segids );

/// @brief Output residue-residue interactions that are not base pairs or base stacks;
void
output_other_contacts( std::ostream & out, pose::Pose const & pose );

/// @brief get rid of Upper and Lower from RNA; useful for cleaner output of annotated_sequence.
void
remove_upper_lower_variants_from_RNA( pose::Pose & pose );

} //ns rna
} //ns pose
} //ns core

#endif
