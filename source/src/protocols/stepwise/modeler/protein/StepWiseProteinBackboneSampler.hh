// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SWA_Screener.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinBackboneSampler_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinBackboneSampler_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/protein/MainChainTorsionSet.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/RamaPrePro.fwd.hh>
#include <string>
#include <map>

//Auto Headers
#include <core/id/TorsionID.fwd.hh>
namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

typedef std::map< core::Size, core::Size > ResMap;
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseProteinBackboneSampler: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseProteinBackboneSampler( protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters );

	//destructor -- necessary?
	~StepWiseProteinBackboneSampler();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;


	void
	set_silent_file( std::string const & setting );

	void
	set_rmsd_cutoff( core::Real const setting );

	void
	set_n_sample( core::Size const setting );

	void
	set_nstruct_centroid( core::Size const setting );

	void
	set_filter_native_big_bins( bool const setting );

	void
	set_centroid_screen( bool const setting );

	void
	set_ghost_loops( bool const setting );

	void
	set_apply_vdw_cut( bool const setting );

	void
	set_centroid_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	void
	set_centroid_score_diff_cut( core::Real const setting );

	utility::vector1< MainChainTorsionSetList > const & main_chain_torsion_set_lists() const;

	utility::vector1< utility::vector1< core::Real > >
	main_chain_torsion_set_lists_real( core::pose::Pose const & pose ) const;

	utility::vector1< core::id::TorsionID > which_torsions( core::pose::Pose const & pose );


	void
	setup_centroid_screen(core::Real const centroid_score_diff_cut,
		std::string const & centroid_weights,
		core::Size const nstruct_centroid,
		bool const ghost_loops);

	core::Size
	get_big_bin( core::Real const phi, core::Real const psi ) const;

	void
	set_moving_residues( utility::vector1< core::Size > const & moving_res );

	void
	set_fixed_residues( utility::vector1< core::Size > const & fixed_res );

	void set_expand_loop_takeoff( bool const setting ){ expand_loop_takeoff_ = setting; }
	bool expand_loop_takeoff() const { return expand_loop_takeoff_; }

private:

	void
	define_moving_res( core::pose::Pose const & pose );

	void
	setup_torsion_sets( core::pose::Pose const & pose );

	void
	filter_main_chain_torsion_sets();

	void
	prepare_ghost_pose( core::pose::Pose const & pose );

	void
	initialize_ghost_pose(
		core::pose::PoseOP & ghost_pose,
		std::string const & desired_sequence,
		core::pose::Pose const & template_pose,
		ResMap const & ghost_map,
		core::kinematics::FoldTree f );

	void
	copy_coords( core::pose::Pose & pose, core::pose::Pose const & template_pose, ResMap const & ghost_map ) const;

	core::kinematics::FoldTree
	figure_out_fold_tree( ResMap const & ghost_map ) const;

	void
	sample_residues_recursively(
		Size const which_res,
		Size & count,
		core::pose::Pose & pose );

	void
	filter_and_save( core::pose::Pose & pose,
		std::string const & tag  );

	void
	get_main_chain_torsion_set_list(
		core::Size const n,
		core::pose::Pose const & pose,
		MainChainTorsionSetList & main_chain_torsion_set_list );

	void
	get_main_chain_torsion_set_list_coarse(
		core::Size const n,
		core::pose::Pose const & pose,
		MainChainTorsionSetList & main_chain_torsion_set_list );

	void
	get_main_chain_torsion_set_list_full( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list );


	void
	get_main_chain_torsion_set_list_n_terminus( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list );


	void
	get_main_chain_torsion_set_list_c_terminus( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list );


	void
	get_main_chain_torsion_set_list_sample_phi_only( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list );


	void
	get_main_chain_torsion_set_list_sample_psi_only( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list );

	void
	filter_native_BIG_BINS(
		core::Size const n,
		MainChainTorsionSetList & main_chain_torsion_set_list );

	void
	filter_based_on_desired_secstruct(
		char const secstruct,
		MainChainTorsionSetList & main_chain_torsion_set_list );

	void
	filter_big_bin( Size const big_bin,
		MainChainTorsionSetList & main_chain_torsion_set_list );


private:

	void
	output_centroid_silent_struct(
		core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
		std::string const & silent_file, std::string const & tag );

	void
	convert_to_centroid( core::pose::Pose & pose );

	void
	sample_cis_omega( MainChainTorsionSetList & main_chain_torsion_set_list );

	void
	initialize_is_fixed_res();

	core::Real
	get_rotamer_angle( core::Size const i, core::Size const N_SAMPLE );

	core::Real
	rama_energy(
		core::Size const n,
		core::pose::Pose const & pose,
		utility::vector1< core::Real > const & mainchain_torsions );

private:

	protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters_;

	utility::vector1< Size > const moving_residues_input_;
	utility::vector1< Size > moving_residues_;
	core::Size n_sample_;
	core::Size n_sample_beta_;
	core::Real rmsd_cutoff_;

	core::scoring::ScoreFunctionOP centroid_scorefxn_;
	std::string silent_file_;
	bool filter_native_big_bins_;
	bool centroid_screen_;
	core::Real centroid_score_ref_;
	core::Real centroid_score_diff_cut_;

	bool apply_vdw_cut_;
	core::Real centroid_vdw_ref_;

	Size nstruct_centroid_;

	core::scoring::Ramachandran const & ramachandran_;
	core::scoring::RamaPrePro const & rama_pre_pro_;

	MainChainTorsionSetList  main_chain_torsion_set_for_moving_residues_;
	utility::vector1< MainChainTorsionSetList > main_chain_torsion_sets_for_moving_residues_;

	utility::vector1< core::Real > centroid_scores_;

	bool ghost_loops_;
	ResMap ghost_map_;
	core::pose::PoseOP ghost_pose_, ghost_native_pose_;

	utility::vector1< bool > is_pre_proline_;

	utility::vector1< bool > is_fixed_res_;
	utility::vector1< bool > is_fixed_res_input_;

	bool expand_loop_takeoff_;

	bool use_rama_pre_pro_;

};

} //protein
} //modeler
} //stepwise
} //protocols

#endif
