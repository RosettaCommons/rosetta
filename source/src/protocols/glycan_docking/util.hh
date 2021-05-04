// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/glycan_docking/util.hh
/// @brief Utility methods for protocols/glycan_docking/GlycanDockProtocol
/// @author Morgan Nance (morganlnance@gmail.com)


#ifndef INCLUDED_protocols_glycan_docking_util_hh
#define INCLUDED_protocols_glycan_docking_util_hh


// Project Headers
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
//#include <protocols/moves/PyMOLMover.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ChainSelector.hh>

// Basic Headers

// Numeric Headers

// Utility Headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>


namespace protocols {
namespace glycan_docking {


////////////////////////////////////////////////////////////////////////////////
/// @brief Creates a ChainSelector based on the downstream partner (glycoligand)
/// @brief of the -docking:partners (e.g. AB_X) after some quality assurance
core::select::residue_selector::ChainSelectorCOP
setup_glycoligand_selector
( std::string const & docking_partners /* e.g. A_X */ );


/// @brief Return a residue subset where carbohydrates of the glycoligand
/// @brief that have (movable) glycosidic torsion angles
/// @note Reducing end carbohydrate does not have dihedrals
utility::vector1< bool >
get_glycolig_subset_with_torsions
( core::pose::Pose const & pose,
	core::select::residue_selector::ChainSelectorCOP glycolig_selector );


/// @brief Create a TaskFactory appropriate for the GlycanDockProtocol
/// @note If -prepack_only is true, TaskFactory is not restricted to interface
core::pack::task::TaskFactoryOP
setup_GlycanDock_taskfactory
( core::Size const docking_jump_num,
	core::Real const interface_packing_distance,
	bool const prepack_only );


/// @brief Create a MoveMap appropriate for structure minimization
/// @note Not affected by -prepack_only because minimization is not used
core::kinematics::MoveMapOP
setup_GlycanDock_movemap
( core::pose::Pose const & pose,
	core::select::residue_selector::ChainSelectorCOP glycolig_selector,
	bool const lock_rings,
	core::Size const dock_jump_num );


/// @brief Create SequenceMover holding Stage 1 rigid-body perturbation Movers
protocols::moves::SequenceMoverOP
setup_GlycanDock_stage1_rb_seqmover
( core::pose::Pose const & pose,
	core::Size const docking_jump_num,
	core::Real const stage1_rotation_magnitude,
	core::Real const stage1_translation_magnitude,
	bool const stage1_rotate_glycolig_about_com );


/// @brief Create RandomMover holding Stage 2 rigid-body perturbation Movers
protocols::moves::RandomMoverOP
setup_GlycanDock_stage2_rb_randmover
( core::Size const docking_jump_num,
	core::Real const stage2_rotation_magnitude,
	core::Real const stage2_translaation_magnitude,
	bool const slide_glycolig_into_contact );


/// @brief Create RandomMover holding Stage 2 torsion angle perturbation Movers
protocols::moves::RandomMoverOP
setup_GlycanDock_stage2_tor_randmover
( core::pose::Pose const & pose,
	utility::vector1< bool > const & glycolig_subset_with_torsions,
	core::scoring::ScoreFunctionOP sf,
	core::Real const mc_kt,
	core::Size const n_shear_moves,
	bool const refine_only );


/// @brief Apply small uniform perturbation on each glycosidic torsion angle
/// @note Stage 1 serves to increase the conformational diversity of
/// @note a docking trajectory by injecting a little initial randomness
void
do_stage1_tor_uniform_perturbation
( core::pose::Pose & pose,
	utility::vector1< bool > const & glycolig_subset,
	core::Real const stage1_tor_perturb_mag );


/// @brief Apply Stage 1 initial conformation perturbation to input Pose
/// @note Stage 1 serves to increase the conformational diversity of
/// @note a docking trajectory by injecting a little initial randomness
void
do_stage1_conformation_initialization
( core::pose::Pose & pose,
	protocols::moves::SequenceMoverOP stage1_rb_seqmover,
	/* protocols::moves::RandomMoverOP stage1_tor_mover TODO */
	utility::vector1< bool > const & glycolig_subset,
	bool const glycolig_has_dihedrals,
	core::Real const stage1_tor_perturb_mag );


/// @brief Apply an inner cycle of Stage 2 sampling and optimization
/// @note Stage 2 is the main docking refinement portion of GlycanDock
core::Size /* n moves accepted */
do_stage2_sample_pack_min_cycle
( core::pose::Pose & pose,
	protocols::moves::RandomMoverOP stage2_mover,
	core::Size const n_rounds,
	minimization_packing::PackRotamersMoverOP full_packer,
	minimization_packing::EnergyCutRotamerTrialsMoverOP ecut_packer,
	core::Size const full_pack_every_x_rounds,
	minimization_packing::MinMoverOP min_mover,
	moves::MonteCarloOP mc );


/// @brief Determine if rsd_ii is in contact with rsd_jj
/// @details If any non-hydrogen atom of rsd_ii is within dist_cutoff
/// @details (generally 5 Ang) of any non-hydrogen atom of rsd_jj,
/// @details then the two residues are considered contacting
/// @note This does not check if an atom is virtual, but virtual
/// @note carbohydrate atoms should track their corresponding real atom
// Essentially the same as protocols/docking/metrics::calc_res_contact
bool
calc_res_contact
( core::conformation::ResidueCOP rsd_ii,
	core::conformation::ResidueCOP rsd_jj,
	core::Real const dist_cutoff );


/// @brief Generate a map of Pose interface contacts between the
/// @brief receptor_nbrhood_subset and the glycolig_subset
/// @details receptor_nbrhood_subset: A pre-filtered residue subset
/// @details of protein receptor residues at the interface based on
/// @details a cutoff of 10 Ang (using the tenA energy map)
/// @details glycolig_subset: Subset indicating the carbohydrate
/// @details residues of the glycoligand
ObjexxFCL::FArray2D_bool
gen_intf_contact_map
( core::pose::Pose const & pose,
	utility::vector1< bool > const & receptor_nbrhood_subset,
	utility::vector1< bool > const & glycolig_subset,
	core::Real const dist_cutoff );


/// @brief Calculate various metrics describing the protein-glycoligand interface
/// @details The return vector1 looks like so:
/// @details < num intf res, num intf res in the ref pose,
/// @details   the fraction of ref pose intf res recovered in the decoy,
/// @details   num intf contacts, num intf contacts in the ref pose,
/// @details   the fraction of ref pose intf contacts recovered in the decoy >
utility::vector1< core::Real >
calc_GlycanDock_intf_metrics
( core::pose::Pose & pose,
	core::pose::Pose const & ref_pose,
	utility::vector1< bool > const & glycolig_subset,
	core::Real const dist_cutoff = 5 /* Ang */);


/// @brief Create an InterfaceAnalyzerMover appropriate for GlycanDock output
protocols::analysis::InterfaceAnalyzerMoverOP
get_GlycanDock_IAM
( std::string const & docking_partners,
	core::scoring::ScoreFunctionOP sf );


} // glycan_docking
} // protocols

#endif //protocols_glycan_docking_util_hh
