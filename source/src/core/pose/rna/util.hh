// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/rna/util.hh
/// @brief
/// @author Rhiju

#ifndef INCLUDED_core_pose_rna_util_hh
#define INCLUDED_core_pose_rna_util_hh

#include <core/types.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/types.hh>
#include <utility/vector1.fwd.hh>

using namespace core::chemical::rna;

namespace core {
namespace pose {
namespace rna{

bool
mutate_position( core::pose::Pose & pose, core::Size const i, char const & new_seq );

void
figure_out_reasonable_rna_fold_tree( core::pose::Pose & pose );

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
		PuckerState const pucker_state = ANY_PUCKER,
		bool const skip_same_state = false,
		bool const idealize_coord = true );

// DEPRECATED
void
apply_ideal_c2endo_sugar_coords(
																pose::Pose & pose,
																Size const i );

void
correctly_position_cutpoint_phosphate_torsions( pose::Pose & current_pose, Size const five_prime_chainbreak );

} //ns rna
} //ns pose
} //ns core

#endif
