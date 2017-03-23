// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/RNA_FragmentMonteCarlo.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_FragmentMonteCarlo_HH
#define INCLUDED_protocols_farna_RNA_FragmentMonteCarlo_HH

#include <protocols/moves/Mover.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.fwd.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_FragmentMover.fwd.hh>
#include <protocols/rna/denovo/fragments/RNA_Fragments.fwd.hh>
#include <protocols/rna/denovo/base_pairs/RNA_BasePairHandler.fwd.hh>
#include <protocols/rna/denovo/setup/RNA_DeNovoPoseInitializer.fwd.hh>
#include <protocols/rna/denovo/libraries/RNA_ChunkLibrary.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_JumpMover.fwd.hh>
#include <protocols/rna/movers/RNA_LoopCloser.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_Relaxer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/rna/StubStubType.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>
#include <numeric/MathNTensor.fwd.hh>

namespace protocols {
namespace rna {
namespace denovo {

class RNA_FragmentMonteCarlo: public protocols::moves::Mover {

public:

	//constructor
	RNA_FragmentMonteCarlo( options::RNA_FragmentMonteCarloOptionsCOP = nullptr );

	//destructor
	~RNA_FragmentMonteCarlo() override;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override { return "RNA_FragmentMonteCarlo"; }

	void
	set_refine_pose( bool const setting ){ refine_pose_ = setting; }

	void
	set_chunk_coverage( core::Real const & setting ){ chunk_coverage_ = setting; }

	void
	set_all_lores_score_final( std::list< core::Real > const & setting ){ all_lores_score_final_ = setting; }

	void
	set_denovo_scorefxn( core::scoring::ScoreFunctionCOP setting ){ denovo_scorefxn_ = setting; }

	void
	set_hires_scorefxn( core::scoring::ScoreFunctionCOP setting ){ hires_scorefxn_ = setting; }

	void
	set_vdw_grid( protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP setting ){ vdw_grid_ = setting; }

	void
	set_is_rna_and_protein( bool const & setting ){ is_rna_and_protein_ = setting; }

	void
	set_out_file_tag( std::string const & setting ){ out_file_tag_ = setting; }

	std::list< core::Real > const & all_lores_score_final() const { return all_lores_score_final_; }

	core::pose::PoseOP lores_pose() { return lores_pose_; }

	core::Real lores_score_early() const { return lores_score_early_; }

	core::Real lores_score_final() const { return lores_score_final_; }

	void
	set_rna_base_pair_handler( base_pairs::RNA_BasePairHandlerCOP setting ) { rna_base_pair_handler_ = setting; }

	base_pairs::RNA_BasePairHandlerCOP rna_base_pair_handler() const { return rna_base_pair_handler_; }

	void
	set_rna_de_novo_pose_initializer( setup::RNA_DeNovoPoseInitializerCOP setting ) { rna_de_novo_pose_initializer_ = setting; }

	void
	set_user_input_chunk_library( libraries::RNA_ChunkLibraryCOP setting ) { user_input_rna_chunk_library_ = setting; }

	libraries::RNA_ChunkLibraryCOP rna_chunk_library() const { return rna_chunk_library_; }

	void show(std::ostream & output) const override;

	void
	align_pose( core::pose::Pose & pose, bool const verbose = false ) const;

	core::Real
	get_rmsd_no_superimpose( core::pose::Pose const & pose ) const;

	core::Real
	get_rmsd_stems_no_superimpose ( core::pose::Pose const & pose ) const;

	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map() const { return atom_level_domain_map_; }

	protocols::rna::movers::RNA_LoopCloserCOP rna_loop_closer() const { return rna_loop_closer_; }

	bool
	loop_modeling() const;

	numeric::MathNTensorOP< core::Size, 6 > jump_histogram() const{ return jump_histogram_; }
	void set_jump_histogram( numeric::MathNTensorOP< core::Size, 6 > setting ) { jump_histogram_ = setting; }

private:

	void
	initialize( core::pose::Pose & pose );

	void
	initialize_libraries( core::pose::Pose & pose );

	void
	initialize_movers();

	void
	initialize_score_functions();

	void
	initialize_parameters();

	void
	do_random_moves( core::pose::Pose & pose );

	void
	randomize_rigid_body_orientations( core::pose::Pose & pose );

	void
	randomize_rnp_rigid_body_orientations( core::pose::Pose & pose );

	void
	update_denovo_scorefxn_weights( core::Size const r );

	core::Size
	figure_out_constraint_separation_cutoff( core::Size const r, core::Size const  max_dist );

	void
	update_pose_constraints( core::Size const r, core::pose::Pose & pose );

	void
	update_frag_size( core::Size const r );

	void
	do_move_trial( core::Size const & i, core::pose::Pose & pose );

	void
	RNA_move_trial( core::pose::Pose & pose );

	void
	random_fragment_trial( core::pose::Pose & pose );

	bool
	random_chunk_trial( core::pose::Pose & pose );

	void
	random_jump_trial( core::pose::Pose & pose );

	void
	setup_rnp_fold_tree( core::pose::Pose & pose );

	void
	rnp_docking_trial( core::pose::Pose & pose );

	core::kinematics::FoldTree
	get_rnp_docking_fold_tree( core::pose::Pose const & pose );

	void
	setup_rna_protein_docking_mover( core::pose::Pose const & pose, core::Size const round );

	void
	setup_rigid_body_mover( core::pose::Pose const & pose, core::Size const r );

	void
	setup_rigid_body_mover2( core::pose::Pose const & pose, core::Size const r );

	bool
	check_score_filter( core::Real const lores_score_, std::list< core::Real > & all_lores_score_ );

	void
	apply_chem_shift_data( core::pose::Pose & pose );

	void
	setup_monte_carlo_cycles( core::pose::Pose const & pose );

	void
	final_score( core::pose::Pose & pose );

	utility::vector1< Size >
	reroot_pose_before_align_and_return_moving_res( core::pose::Pose & pose ) const;

	core::Real
	get_rmsd( core::pose::Pose const & const_pose ) const;

	void
	check_for_loop_modeling_case( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ) const;

	void
	initialize_output_score();

	void
	output_score_if_desired( core::Size const & r,
		core::Size const & i,
		core::pose::Pose & pose );

	void
	finish_output_score();

	core::kinematics::RT
	get_output_jump_RT( core::pose::PoseCOP pose ) const;

	void
	get_output_jump_stub_stub( core::pose::Pose const & pose,
		core::kinematics::Stub & stub1,
		core::kinematics::Stub & stub2 ) const;

	void
	output_jump_information( core::pose::Pose const & pose );

private:

	// The parameters in this OptionsCOP should not change:
	options::RNA_FragmentMonteCarloOptionsCOP options_;
	std::string out_file_tag_;

	// Movers (currently must be set up outside, but should write auto-setup code)
	base_pairs::RNA_BasePairHandlerCOP rna_base_pair_handler_;
	libraries::RNA_ChunkLibraryCOP user_input_rna_chunk_library_;
	libraries::RNA_ChunkLibraryOP rna_chunk_library_;
	setup::RNA_DeNovoPoseInitializerCOP rna_de_novo_pose_initializer_;
	movers::RNA_FragmentMoverOP rna_fragment_mover_;
	movers::RNA_JumpMoverOP rna_jump_mover_;
	protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer_;
	movers::RNA_MinimizerOP rna_minimizer_;
	movers::RNA_RelaxerOP rna_relaxer_;
	protocols::rigid::RigidBodyPerturbMoverOP rnp_docking_mover_;
	protocols::rigid::RigidBodyPerturbMoverOP rigid_body_mover_;

	protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map_;

	core::scoring::ScoreFunctionCOP denovo_scorefxn_;
	core::scoring::ScoreFunctionCOP hires_scorefxn_;
	core::scoring::ScoreFunctionCOP chem_shift_scorefxn_;
	core::scoring::ScoreFunctionCOP final_scorefxn_;
	core::scoring::ScoreFunctionOP working_denovo_scorefxn_;

	protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_grid_;

	// Parameters that change during run:
	protocols::moves::MonteCarloOP monte_carlo_;
	core::Size const monte_carlo_cycles_max_default_;
	core::Size monte_carlo_cycles_;
	core::Size rounds_;
	core::Size frag_size_;
	bool do_close_loops_;
	bool refine_pose_;
	core::Real jump_change_frequency_;
	core::Real lores_score_early_;
	core::Real lores_score_final_;
	core::Real chunk_coverage_;

	std::list< core::Real > all_lores_score_final_;
	core::scoring::constraints::ConstraintSetOP constraint_set_;
	core::pose::PoseCOP align_pose_;

	core::pose::PoseOP lores_pose_;
	bool is_rna_and_protein_;
	bool do_rnp_docking_;
	core::kinematics::FoldTree rnp_docking_ft_;

	utility::io::ozstream running_score_output_;
	numeric::MathNTensorOP< core::Size, 6 > jump_histogram_;
	utility::vector1< core::Real > jump_histogram_min_, jump_histogram_max_, jump_histogram_bin_width_;
	core::pose::rna::StubStubType output_stub_stub_type_;
	core::kinematics::RTOP reference_RT_;

};

} //denovo
} //rna
} //protocols

#endif
