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
#include <protocols/rna/denovo/movers/RNA_DeNovoMasterMover.fwd.hh>
#include <core/import_pose/RNA_BasePairHandler.fwd.hh>
#include <protocols/rna/denovo/RNA_DeNovoPoseInitializer.fwd.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.fwd.hh>
#include <protocols/rna/movers/RNA_LoopCloser.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_Relaxer.fwd.hh>
#include <protocols/rna/denovo/movers/RNP_HighResMover.fwd.hh>
#include <protocols/rna/denovo/movers/RNA_HelixMover.fwd.hh>
#include <protocols/rna/denovo/output/RNA_FragmentMonteCarloOutputter.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.fwd.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/format.hh>

namespace protocols {
namespace rna {
namespace denovo {

class RNA_FragmentMonteCarlo: public protocols::moves::Mover {

public:

	//constructor
	RNA_FragmentMonteCarlo( core::import_pose::options::RNA_FragmentMonteCarloOptionsCOP = nullptr );

	//destructor
	~RNA_FragmentMonteCarlo() override;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override { return "RNA_FragmentMonteCarlo"; }

	void
	set_refine_pose( bool const setting ){ refine_pose_ = setting; }

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
	core::pose::PoseOP get_additional_output() override;

	core::Real lores_score_early() const { return lores_score_early_; }

	core::Real lores_score_final() const { return lores_score_final_; }

	void
	set_rna_base_pair_handler( core::import_pose::RNA_BasePairHandlerCOP setting ) { rna_base_pair_handler_ = setting; }

	core::import_pose::RNA_BasePairHandlerCOP rna_base_pair_handler() const { return rna_base_pair_handler_; }

	void
	set_rna_de_novo_pose_initializer( protocols::rna::denovo::RNA_DeNovoPoseInitializerCOP setting ) { rna_de_novo_pose_initializer_ = setting; }

	void
	set_user_input_chunk_library( core::import_pose::libraries::RNA_ChunkLibraryCOP setting ) { user_input_rna_chunk_library_ = setting; }

	void
	set_user_input_chunk_initialization_library( core::import_pose::libraries::RNA_ChunkLibraryCOP setting ) { user_input_rna_chunk_initialization_library_ = setting; }

	core::import_pose::libraries::RNA_ChunkLibraryCOP rna_chunk_library() const { return rna_chunk_library_; }

	void show(std::ostream & output) const override;

	void
	align_pose( core::pose::Pose & pose, bool const verbose = false ) const;

	core::Real
	get_rmsd_no_superimpose( core::pose::Pose const & pose ) const;

	core::Real
	get_rmsd_stems_no_superimpose ( core::pose::Pose const & pose ) const;

	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map() const { return atom_level_domain_map_; }

	protocols::rna::movers::RNA_LoopCloserCOP rna_loop_closer() const { return rna_loop_closer_; }

	bool
	loop_modeling() const;

	void set_outputter( output::RNA_FragmentMonteCarloOutputterOP const & setting ){ outputter_ = setting; }
	output::RNA_FragmentMonteCarloOutputterOP outputter() const { return outputter_; }

private:

	void
	initialize( core::pose::Pose & pose );

	void
	initialize_libraries( core::pose::Pose & pose );

	void
	initialize_movers( core::pose::Pose const & pose );
	// initialize_movers();

	void
	initialize_score_functions();

	void
	initialize_parameters();

	void
	update_denovo_scorefxn_weights( core::Size const r );

	core::Size
	figure_out_constraint_separation_cutoff( core::Size const r, core::Size const  max_dist );

	void
	update_pose_constraints( core::Size const r, core::pose::Pose & pose );

	void
	update_rna_denovo_master_mover( core::Size const & r,
		core::pose::Pose const & pose );

	core::Size
	update_frag_size( core::Size const & r ) const;

	void
	get_rigid_body_move_mags( core::Size const & r,
		core::Real & rot_mag,
		core::Real & trans_mag,
		core::Real const & rot_mag_init,
		core::Real const & trans_mag_init,
		core::Real const & rot_mag_final,
		core::Real const & trans_mag_final ) const;

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
	setup_full_initial_structure( core::pose::Pose & pose ) const;

	core::Real
	randomize_and_close_all_chains( core::pose::Pose & pose ) const;

	void
	copy_structure_keep_fold_tree( core::pose::Pose & pose,
		core::pose::Pose const & pose_to_copy ) const;

private:

	// The parameters in this OptionsCOP should not change:
	core::import_pose::options::RNA_FragmentMonteCarloOptionsCOP options_;
	std::string out_file_tag_;

	// Movers (currently must be set up outside, but should write auto-setup code)
	core::import_pose::RNA_BasePairHandlerCOP rna_base_pair_handler_;
	core::import_pose::libraries::RNA_ChunkLibraryCOP user_input_rna_chunk_library_;
	core::import_pose::libraries::RNA_ChunkLibraryCOP user_input_rna_chunk_initialization_library_;
	core::import_pose::libraries::RNA_ChunkLibraryOP rna_chunk_library_;
	core::import_pose::libraries::RNA_ChunkLibraryOP rna_chunk_initialization_library_;
	protocols::rna::denovo::RNA_DeNovoPoseInitializerCOP rna_de_novo_pose_initializer_;
	protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer_;
	protocols::rna::movers::RNA_LoopCloserOP rna_loop_closer_init_;
	movers::RNA_DeNovoMasterMoverOP rna_denovo_master_mover_;
	movers::RNA_DeNovoMasterMoverOP rna_denovo_master_mover_init_;
	movers::RNA_HelixMoverOP rna_helix_mover_;
	movers::RNA_MinimizerOP rna_minimizer_;
	movers::RNA_RelaxerOP rna_relaxer_;
	movers::RNP_HighResMoverOP rnp_high_res_mover_;

	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map_;

	core::scoring::ScoreFunctionCOP denovo_scorefxn_;
	core::scoring::ScoreFunctionCOP hires_scorefxn_;
	core::scoring::ScoreFunctionCOP chem_shift_scorefxn_;
	core::scoring::ScoreFunctionOP chainbreak_sfxn_;
	core::scoring::ScoreFunctionCOP final_scorefxn_;
	core::scoring::ScoreFunctionOP working_denovo_scorefxn_;

	protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_grid_;

	// Parameters that change during run:
	core::Size const monte_carlo_cycles_max_default_;
	core::Size monte_carlo_cycles_;
	core::Size rounds_;
	bool refine_pose_;
	core::Real lores_score_early_;
	core::Real lores_score_final_;

	std::list< core::Real > all_lores_score_final_;
	core::scoring::constraints::ConstraintSetOP constraint_set_;
	core::pose::PoseCOP align_pose_;

	core::pose::PoseOP lores_pose_;
	bool is_rna_and_protein_;

	output::RNA_FragmentMonteCarloOutputterOP outputter_;

};

} //denovo
} //rna
} //protocols

#endif
