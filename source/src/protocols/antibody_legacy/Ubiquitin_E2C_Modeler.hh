// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Ubiquitin_E2C_Modeler.hh
/// @brief Build model for di-ubiquitin and e2c enzyme complex
/// @details
///
///
/// @author Aroop Sircar


#ifndef INCLUDED_protocols_antibody_legacy_Ubiquitin_E2C_Modeler_hh
#define INCLUDED_protocols_antibody_legacy_Ubiquitin_E2C_Modeler_hh
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/antibody_legacy/Ubiquitin_E2C_Modeler.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace ub_e2c {

class ubi_e2c_modeler: public moves::Mover {
public:

	// default constructor
	ubi_e2c_modeler();

	// default destructor
	~ubi_e2c_modeler() override;

	void set_default();

	void apply( core::pose::Pose & pose_in ) override;
	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;

	// find out domain ends and center of masses
	void
	setup_key_residues( const core::pose::Pose & pose_in );

	void setup_move_maps();

	void setup_complex_fold_tree(
		core::pose::Pose & pose_in,
		bool trim = false );

	void initial_cter_perturbation( core::pose::Pose & pose_in );

	void setup_simple_fold_tree(
		core::Size jumppoint1,
		core::Size cutpoint,
		core::Size jumppoint2,
		core::Size nres,
		core::pose::Pose & pose_in );

	void trim_cter( core::pose::Pose & pose_in );

	void restore_cter(
		core::pose::Pose & pose_in,
		core::pose::Pose without_cter );

	void init_k48r_perturbation( core::pose::Pose & pose_in );

	void init_d77_perturbation( core::pose::Pose & pose_in );

	core::Real initial_perturbation( core::pose::Pose & pose_in );

	core::Real centroid_mode_perturbation( core::pose::Pose & pose_in );

	core::Real fullatom_mode_perturbation( core::pose::Pose & pose_in );

	void initial_repack( core::pose::Pose & pose_in );

	void setup_packer_task( core::pose::Pose & pose_in );

	void restrict_to_interfacial_loop_packing( core::pose::Pose & pose_in );

	void set_e2g2_diubi_fold_tree( core::pose::Pose & pose_in );

	core::Real calc_interaction_energy(
		const core::pose::Pose & pose_in,
		bool dimer = true );

	core::Real CSP_fraction(
		const core::pose::Pose & pose_in,
		bool non_CSP = false,
		bool trim = false,
		bool swap = false );

	bool centroid_filter( core::pose::Pose & pose_in );

	bool fullatom_filter( core::pose::Pose & pose_in );

	core::Real
	calc_Lrmsd(
		const core::pose::Pose & pose_in,
		const core::pose::Pose & native_pose,
		core::Size ubiquitin );

	void evaluate_native( core::pose::Pose & pose_in );

	void optimize_cov_bond( core::pose::Pose & pose_in );

private:

	void compute_trim_CSPs();

	void compute_swap_trim_CSPs();

	void assign_CSPs( const core::pose::Pose & pose_in );

	void assign_non_CSPs( const core::pose::Pose & pose_in );

	core::Size e2_ctr_of_mass_;
	core::Size k48r_ctr_of_mass_;
	core::Size d77_ctr_of_mass_;
	core::Size e2_end_;
	core::Size k48r_end_;
	core::Size d77_end_;
	core::Size k48r_48_lys_;
	core::Size d77_48_lys_;
	core::Size k48r_trim_end_;
	core::Size d77_trim_end_;
	core::Size d77_trim_ctr_mass_;

	bool k48r_swap_;

	// move maps
	core::kinematics::MoveMapOP k48r_docking_map_;
	core::kinematics::MoveMapOP d77_docking_map_;
	core::kinematics::MoveMapOP docking_map_;
	core::kinematics::MoveMapOP flex_cter_map_;
	core::kinematics::MoveMapOP all_dof_map_;
	core::kinematics::MoveMapOP init_k48r_docking_map_;
	core::kinematics::MoveMapOP init_d77_docking_map_;
	core::kinematics::MoveMapOP init_docking_map_;
	core::kinematics::MoveMapOP init_flex_cter_map_;
	core::kinematics::MoveMapOP init_all_dof_map_;


	// jump between E2 G2 and UbK48R
	core::Size const e2_k48r_jump_;
	// jump between E2 G2 and UbD77
	core::Size const e2_d77_jump_;
	core::Real const max_k48_cter_dist_;
	core::Real const temperature_;
	std::map < std::string, core::Real > score_map_;
	// number of flexible c-terminal residues in Ubiquitin
	// actually its one more than flex_cter residues
	// so if flex_cter equals 3, then 4 residues are flexible
	core::Size const flex_cter_;
	// CA atom index
	core::Size const CA;
	// If a decoy does not pass the filters in fullatom perturbation
	// how many attempts to make before failure
	core::Size const max_repeats_;
	// strict constraint cutoffs for fullatom
	// 0.25 per constraint
	core::Real fullatom_constraint_cutoff_; // depends on mode
	// CSP Satisfaction Factor
	core::Real const centroid_allowed_CSP_fraction_;
	core::Real const fullatom_allowed_CSP_fraction_;
	// CSP score weights
	core::Real const centroid_CSP_weight_;
	core::Real const centroid_non_CSP_weight_;
	core::Real const fullatom_CSP_weight_;
	core::Real const fullatom_non_CSP_weight_;
	// weighting terms
	core::Real const cen_vdw_;
	core::Real const cen_constraint_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real const full_vdw_;
	core::Real full_constraint_; // depends on which mode being run
	// mode
	bool const refinement_mode_;
	bool const cov_bond_only_flag_;
	bool const monoub_mode_;
	bool const higher_d77_pert_mode_;
	// MTSL PRE Tags
	core::Size k48r_48k_mtsl_ ;
	core::Size d77_75g_mtsl_;
	// filters
	bool passed_centroid_filter_;
	bool passed_fullatom_filter_;
	// Applied Full Atom Perturbation
	bool applied_fullatom_pert_;
	//packer task
	core::pack::task::TaskFactoryOP tf_;
	core::pack::task::TaskFactoryOP init_task_factory_;

	//sets up minimization parameters
	core::Real min_tolerance_;
	std::string min_type_;
	bool nb_list_;

	// CSPs
	utility::vector1<core::Size> CSP_;
	utility::vector1<core::Size> CSP_trim_;
	utility::vector1<core::Size> CSP_swap_trim_;

	// non_CSPs
	utility::vector1<core::Size> non_CSP_;
	utility::vector1<core::Size> non_CSP_trim_;
	utility::vector1<core::Size> non_CSP_swap_trim_;

	// score functions
	core::scoring::ScoreFunctionOP dock_lowres_scorefxn_;
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP dock_lowres_cst_scorefxn_;
	core::scoring::ScoreFunctionOP lowres_cst_scorefxn_;
	core::scoring::ScoreFunctionOP pack_scorefxn_;
	core::scoring::ScoreFunctionOP dockfa_scorefxn_;
	core::scoring::ScoreFunctionOP dockfa_cst_scorefxn_;
	core::scoring::ScoreFunctionOP dockfa_min_scorefxn_;
	core::scoring::ScoreFunctionOP dockfa_cst_min_scorefxn_;
	core::scoring::ScoreFunctionOP pack_cst_scorefxn_;
	core::scoring::ScoreFunctionOP output_cen_scorefxn_;
	core::scoring::ScoreFunctionOP output_full_scorefxn_;

	// Mono Ubiquitin Mode Variables
	void monoub_assign_CSPs(const core::pose::Pose & pose_in );

	void monoub_apply( core::pose::Pose & pose_in );

	void monoub_setup_key_residues( const core::pose::Pose & pose_in );

	void monoub_setup_move_maps();

	void monoub_fold_tree( core::pose::Pose & pose_in );

	void monoub_initial_cter_perturbation( core::pose::Pose & pose_in );

	void monoub_first_perturbation( core::pose::Pose & pose_in );

	void monoub_initial_perturbation( core::pose::Pose & pose_in );

	void monoub_centroid_mode_perturbation( core::pose::Pose & pose_in );

	void monoub_fullatom_mode_perturbation( core::pose::Pose & pose_in );

	core::Real monoub_calc_interaction_energy(
		const core::pose::Pose & pose_in );

	core::Real monoub_CSP_fraction( const core::pose::Pose & pose_in );

	bool monoub_centroid_filter( core::pose::Pose & pose_in );

	bool monoub_fullatom_filter( core::pose::Pose & pose_in );

	core::Real monoub_calc_Lrmsd (
		const core::pose::Pose & pose_in,
		const core::pose::Pose & native_pose );


	core::Size monoub_end_;
	core::Size monoub_ctr_of_mass_;

	core::kinematics::MoveMapOP monoub_all_dof_map_;
	core::kinematics::MoveMapOP monoub_docking_map_;
	core::kinematics::MoveMapOP monoub_flex_cter_map_;
	core::kinematics::MoveMapOP init_monoub_all_dof_map_;
	core::kinematics::MoveMapOP init_monoub_docking_map_;
	core::kinematics::MoveMapOP init_monoub_flex_cter_map_;

}; // class ubi_e2c_modeler


} // namespace ub_e2c
} // namespace protocols
#endif
