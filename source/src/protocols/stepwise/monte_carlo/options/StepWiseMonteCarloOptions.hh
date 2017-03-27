// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseMonteCarloOptions_HH
#define INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseMonteCarloOptions_HH

#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.fwd.hh>
#include <protocols/stepwise/options/StepWiseBasicOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <core/types.hh>

#include <utility/tag/Tag.fwd.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace options {

class StepWiseMonteCarloOptions: virtual
	public protocols::stepwise::options::StepWiseBasicOptions,
	public protocols::stepwise::monte_carlo::mover::options::StepWiseMoveSelectorOptions {

public:

	//constructor
	StepWiseMonteCarloOptions();

	StepWiseMonteCarloOptions( StepWiseMonteCarloOptions const & src );

	//destructor
	~StepWiseMonteCarloOptions();

public:

	/// @brief Describe this instance to a given output stream
	virtual
	void
	show(
		std::ostream & ) const{}

	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP
	){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "StepWiseMonteCarloOptions";}

public:

	StepWiseMonteCarloOptionsOP clone() const;

	void
	initialize_from_command_line();

	protocols::stepwise::modeler::options::StepWiseModelerOptionsOP setup_modeler_options() const;

	bool const & verbose_scores() const { return verbose_scores_; }
	void set_verbose_scores( bool const & setting ){ verbose_scores_ = setting; };

	bool const & integration_test_mode() const { return integration_test_mode_; }
	void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; };

	bool const & force_centroid_interaction() const { return force_centroid_interaction_; }
	void set_force_centroid_interaction( bool const & setting ){ force_centroid_interaction_ = setting; }

	core::Real const & sampler_max_centroid_distance() const { return sampler_max_centroid_distance_; }
	void set_sampler_max_centroid_distance( core::Real const & setting ){ sampler_max_centroid_distance_ = setting; }

	bool const & use_phenix_geo() const { return use_phenix_geo_; }
	void set_use_phenix_geo( bool const & setting ){ use_phenix_geo_ = setting; }

	bool const & skip_deletions() const { return skip_deletions_; }
	void set_skip_deletions( bool const & setting ){ skip_deletions_ = setting; }

	bool const & erraser() const { return erraser_; }
	void set_erraser(  bool const & setting ){ erraser_ = setting; }

	core::Size const & cycles() const { return cycles_; }
	void set_cycles( core::Size const & setting ){ cycles_ = setting; }

	core::Real const & add_proposal_density_factor() const { return add_proposal_density_factor_; }
	void set_add_proposal_density_factor( core::Real const & setting ){ add_proposal_density_factor_ = setting; }

	core::Real const & minimize_single_res_frequency() const { return minimize_single_res_frequency_; }
	void set_minimize_single_res_frequency(  core::Real const & setting ){ minimize_single_res_frequency_ = setting; }

	core::Real const & just_min_after_mutation_frequency() const { return just_min_after_mutation_frequency_; }
	void set_just_min_after_mutation_frequency(  core::Real const & setting ){ just_min_after_mutation_frequency_ = setting; }

	core::Real const & temperature() const { return temperature_; }
	void set_temperature(  core::Real const & setting ){ temperature_ = setting; }

	core::Real const & max_missing_weight() const { return max_missing_weight_; }
	void set_max_missing_weight( core::Real const & setting ){ max_missing_weight_ = setting; }

	core::Real const & chainbreak_weight() const { return chainbreak_weight_; }
	void set_chainbreak_weight(  core::Real const & setting ){ chainbreak_weight_ = setting; }

	bool const & allow_split_off() const { return allow_split_off_; }
	void set_allow_split_off( bool const & setting ){ allow_split_off_ = setting; }

	bool const & virtual_sugar_keep_base_fixed() const { return virtual_sugar_keep_base_fixed_; }
	void set_virtual_sugar_keep_base_fixed(  bool const & setting ){ virtual_sugar_keep_base_fixed_ = setting; }

	bool const & virtual_sugar_do_minimize() const { return virtual_sugar_do_minimize_; }
	void set_virtual_sugar_do_minimize(  bool const & setting ){ virtual_sugar_do_minimize_ = setting; }

	bool const & make_movie() const { return make_movie_; }
	void set_make_movie( bool const & setting ){ make_movie_ = setting; }

	bool const & sampler_perform_phosphate_pack() const { return sampler_perform_phosphate_pack_; }
	void set_sampler_perform_phosphate_pack( bool const & setting ){ sampler_perform_phosphate_pack_ = setting; }

	bool const & force_phosphate_instantiation() const { return force_phosphate_instantiation_; }
	void set_force_phosphate_instantiation( bool const & setting ){ force_phosphate_instantiation_ = setting; }

	bool const & virtualize_packable_moieties_in_screening_pose() const { return virtualize_packable_moieties_in_screening_pose_; }
	void set_virtualize_packable_moieties_in_screening_pose( bool const & setting ){ virtualize_packable_moieties_in_screening_pose_ = setting; }

	bool const & rebuild_bulge_mode() const { return rebuild_bulge_mode_; }
	void set_rebuild_bulge_mode( bool const & setting ){ rebuild_bulge_mode_ = setting; }

	bool const & tether_jump() const { return tether_jump_; }
	void set_tether_jump( bool const & setting ){ tether_jump_ = setting; }

	bool const & local_redock_only() const { return local_redock_only_; }
	void set_local_redock_only( bool const & setting ){ local_redock_only_ = setting; }

	bool const & skip_coord_constraints() const { return skip_coord_constraints_; }
	void set_skip_coord_constraints( bool const & setting ){ skip_coord_constraints_ = setting; }

	bool const & filter_native_big_bins() const { return filter_native_big_bins_; }
	void set_filter_native_big_bins( bool const & setting ){ filter_native_big_bins_ = setting; }

	bool const & allow_virtual_o2prime_hydrogens() const { return allow_virtual_o2prime_hydrogens_; }
	void set_allow_virtual_o2prime_hydrogens( bool const & setting ){ allow_virtual_o2prime_hydrogens_ = setting; }

	bool const & allow_virtual_side_chains() const { return allow_virtual_side_chains_; }
	void set_allow_virtual_side_chains( bool const & setting ){ allow_virtual_side_chains_ = setting; }

	Size const & n_sample() const { return n_sample_; }
	void set_n_sample( Size const & setting ){ n_sample_ = setting; }

	bool const & protein_prepack() const { return protein_prepack_; }
	void set_protein_prepack( bool const & setting ){ protein_prepack_ = setting; }

	bool const & o2prime_legacy_mode() const { return o2prime_legacy_mode_; }
	void set_o2prime_legacy_mode( bool const & setting ){ o2prime_legacy_mode_ = setting; }

	bool const & recover_low() const { return recover_low_; }
	void set_recover_low( bool const & setting ){ recover_low_ = setting; }

	bool const & enumerate() const { return enumerate_; }
	void set_enumerate( bool const & setting ){ enumerate_ = setting; }

	bool const & preminimize() const { return preminimize_; }
	void set_preminimize( bool const & setting ){ preminimize_ = setting; }

	bool const & skip_preminimize() const { return skip_preminimize_; }
	void set_skip_preminimize( bool const & setting ){ skip_preminimize_ = setting; }

	bool const & new_move_selector() const { return new_move_selector_; }
	void set_new_move_selector( bool const & setting ){ new_move_selector_ = setting; }

	bool const & test_all_moves() const { return test_all_moves_; }
	void set_test_all_moves( bool const & setting ){ test_all_moves_ = setting; }

	bool const & save_times() const { return save_times_; }
	void set_save_times( bool const & setting ){ save_times_ = setting; }

	bool const & use_precomputed_library() const { return use_precomputed_library_; }
	void set_use_precomputed_library( bool const & setting ){ use_precomputed_library_ = setting; }

	bool const & minimize_after_delete() const { return minimize_after_delete_; }
	void set_minimize_after_delete( bool const & setting ){ minimize_after_delete_ = setting; }

	bool const & use_first_jump_for_submotif() const { return use_first_jump_for_submotif_; }
	void set_use_first_jump_for_submotif( bool const & setting ){ use_first_jump_for_submotif_ = setting; }

	bool checkpoint() const { return checkpoint_; }
	void set_checkpoint( bool const setting ){ checkpoint_ = setting; }

	Size checkpointing_frequency() const { return checkpointing_frequency_; }
	void set_checkpointing_frequency( Size const setting ) { checkpointing_frequency_ = setting; }
	
	bool continue_until_none_missing() const { return continue_until_none_missing_; }
	void set_continue_until_none_missing( bool const setting ){ continue_until_none_missing_ = setting; }


private:

	bool verbose_scores_;
	bool integration_test_mode_;
	bool force_centroid_interaction_;
	core::Real sampler_max_centroid_distance_;
	bool use_phenix_geo_;
	bool skip_deletions_;
	bool erraser_;
	core::Size cycles_;
	core::Real add_proposal_density_factor_;
	core::Real minimize_single_res_frequency_;
	core::Real just_min_after_mutation_frequency_;
	core::Real temperature_;
	core::Real max_missing_weight_;
	core::Real chainbreak_weight_;
	bool allow_split_off_;
	bool virtual_sugar_keep_base_fixed_;
	bool virtual_sugar_do_minimize_;
	bool make_movie_;
	bool sampler_perform_phosphate_pack_;
	bool force_phosphate_instantiation_;
	bool virtualize_packable_moieties_in_screening_pose_;
	bool rebuild_bulge_mode_;
	bool tether_jump_;
	bool local_redock_only_;
	bool skip_coord_constraints_;
	bool filter_native_big_bins_;
	bool allow_virtual_o2prime_hydrogens_;
	bool allow_virtual_side_chains_;
	Size n_sample_;
	bool protein_prepack_;
	bool o2prime_legacy_mode_;
	bool recover_low_;
	bool enumerate_;
	bool preminimize_;
	bool skip_preminimize_;
	bool new_move_selector_;
	bool test_all_moves_;
	bool save_times_;
	bool use_precomputed_library_;
	bool minimize_after_delete_;
	bool use_first_jump_for_submotif_;
	bool checkpoint_;
	core::Size checkpointing_frequency_;
	bool continue_until_none_missing_;

};

} //options
} //monte_carlo
} //stepwise
} //protocols

#endif
