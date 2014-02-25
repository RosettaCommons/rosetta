// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_Screener.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_StepWiseScreener_hh
#define INCLUDED_protocols_stepwise_StepWiseScreener_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/stepwise/MainChainTorsionSet.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/Ramachandran.hh>
#include <string>
#include <map>

//Auto Headers


namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

typedef std::map< std::string, core::pose::PoseOP > PoseList;
typedef std::map< core::Size, core::Size > ResMap;
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseScreener: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseScreener(
		utility::vector1< Size > const & moving_residues
	);

	//destructor -- necessary? -- YES destructors are necessary.
	~StepWiseScreener();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void
	set_silent_file( std::string const & setting );

	void
	set_rmsd_cutoff( core::Real const & setting );

	void
	set_n_sample( core::Size const & setting );

	void
	set_nstruct_centroid( core::Size const & setting );

	void
	set_filter_native_big_bins( bool const & setting );

	void
	set_centroid_screen( bool const & setting );

	void
	set_ghost_loops( bool const & setting );

	void
	set_apply_vdw_cut( bool const & setting );

	void
	set_centroid_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	void
	set_centroid_score_diff_cut( core::Real const & setting );

	utility::vector1< MainChainTorsionSetList > const & main_chain_torsion_set_lists() const;

private:

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
		core::kinematics::FoldTree f
	);

	void
	copy_coords( core::pose::Pose & pose, core::pose::Pose const & template_pose, ResMap const & ghost_map ) const;

	core::kinematics::FoldTree
	figure_out_fold_tree( ResMap const & ghost_map ) const;

	void
	sample_residues_recursively(
		Size const which_res,
		Size & count,
		core::pose::Pose & pose
	);

	void
	filter_and_output(
		core::pose::Pose & pose,
		std::string const & tag
	);

	void
	get_main_chain_torsion_set_list(
		core::Size const & n,
		core::pose::Pose const & pose,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);

	void
	get_main_chain_torsion_set_list_coarse(
		core::Size const & n,
		core::pose::Pose const & pose,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);

	void
	get_main_chain_torsion_set_list_full(
		core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);


	void
	get_main_chain_torsion_set_list_n_terminus(
		core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);


	void
	get_main_chain_torsion_set_list_c_terminus(
		core::Size const & n,
		core::pose::Pose const & pose,
		core::Real const & best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);


	void
	get_main_chain_torsion_set_list_sample_phi_only(
		core::Size const & n,
		core::pose::Pose const & pose,
		core::Real const & best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);


	void
	get_main_chain_torsion_set_list_sample_psi_only(
		core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);

	void
	filter_native_BIG_BINS(
		core::Size const & n,
		MainChainTorsionSetList & main_chain_torsion_set_list
	);

private:

	core::Size
	get_big_bin( core::Real const phi, core::Real const psi ) const;

	void
	output_centroid_silent_struct(
																core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
																std::string const silent_file, std::string const & tag );

	void
	convert_to_centroid( core::pose::Pose & pose );


private:

	utility::vector1< Size > moving_residues_;
	core::Size n_sample_;
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

	core::scoring::Ramachandran ramachandran_;

	MainChainTorsionSetList main_chain_torsion_set_for_moving_residues_;
	utility::vector1< MainChainTorsionSetList > main_chain_torsion_sets_for_moving_residues_;
	utility::vector1< core::Real > centroid_scores_;

	bool ghost_loops_;
	ResMap ghost_map_;
	core::pose::PoseOP ghost_pose_, ghost_native_pose_;

};

} //protein
} //sampling
} //stepwise
} //protocols

#endif
