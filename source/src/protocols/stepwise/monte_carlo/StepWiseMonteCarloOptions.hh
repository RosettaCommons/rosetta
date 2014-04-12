// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarloOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseMonteCarloOptions_HH
#define INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseMonteCarloOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarloOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#include <utility/tag/Tag.fwd.hh>

#ifdef WIN32
	#include <utility/tag/Tag.hh>
#endif


// using namespace core;
// Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
// by our coding convention due to problems it create on modern compilers and because of the name clashing.
// For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	class StepWiseMonteCarloOptions: public basic::resource_manager::ResourceOptions {

	public:

		//constructor
		StepWiseMonteCarloOptions();

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
		type() const{ return "StepWiseRNA_ModelerOptions";}

	public:

		void
		initialize_from_command_line();

		protocols::stepwise::sampling::rna::StepWiseRNA_ModelerOptionsOP setup_rna_modeler_options() const;

		protocols::stepwise::sampling::protein::StepWiseProteinModelerOptionsOP setup_protein_modeler_options() const;

		std::string const & silent_file() const { return silent_file_; }
		void set_silent_file( std::string const & setting ){ silent_file_ = setting; };

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

		bool const & allow_internal_hinge_moves() const { return allow_internal_hinge_moves_; }
		void set_allow_internal_hinge_moves( bool const & setting ){ allow_internal_hinge_moves_ = setting; }

		bool const & allow_internal_local_moves() const { return allow_internal_local_moves_; }
		void set_allow_internal_local_moves( bool const & setting ){ allow_internal_local_moves_ = setting; }

		core::Size const & num_random_samples() const { return num_random_samples_; }
		void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }

		core::Size const & num_pose_minimize() const { return num_pose_minimize_; }
		void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

		core::Size const & cycles() const { return cycles_; }
		void set_cycles( core::Size const & setting ){ cycles_ = setting; }

		core::Real const & add_delete_frequency() const { return add_delete_frequency_; }
		void set_add_delete_frequency( core::Real const & setting ){ add_delete_frequency_ = setting; }

		core::Real const & intermolecular_frequency() const { return intermolecular_frequency_; }
		void set_intermolecular_frequency( core::Real const & setting ){ intermolecular_frequency_ = setting; }

		core::Real const & minimize_single_res_frequency() const { return minimize_single_res_frequency_; }
		void set_minimize_single_res_frequency(  core::Real const & setting ){ minimize_single_res_frequency_ = setting; }

		bool const & minimizer_allow_variable_bond_geometry() const { return minimizer_allow_variable_bond_geometry_; }
		void set_minimizer_allow_variable_bond_geometry( bool const & setting ){ minimizer_allow_variable_bond_geometry_ = setting; }

		core::Real const & minimizer_vary_bond_geometry_frequency() const { return minimizer_vary_bond_geometry_frequency_; }
		void set_minimizer_vary_bond_geometry_frequency( core::Real const & setting ){ minimizer_vary_bond_geometry_frequency_ = setting; }

		core::Real const & switch_focus_frequency() const { return switch_focus_frequency_; }
		void set_switch_focus_frequency( core::Real const & setting ){ switch_focus_frequency_ = setting; }

		core::Real const & just_min_after_mutation_frequency() const { return just_min_after_mutation_frequency_; }
		void set_just_min_after_mutation_frequency(  core::Real const & setting ){ just_min_after_mutation_frequency_ = setting; }

		core::Real const & temperature() const { return temperature_; }
		void set_temperature(  core::Real const & setting ){ temperature_ = setting; }

		utility::vector1< Size > const & bulge_res() const { return bulge_res_; }
		void set_bulge_res( utility::vector1< Size > const & setting ){ bulge_res_ = setting; }

		core::Real const & max_missing_weight() const { return max_missing_weight_; }
		void set_max_missing_weight( core::Real const & setting ){ max_missing_weight_ = setting; }

		core::Real const & chainbreak_weight() const { return chainbreak_weight_; }
		void set_chainbreak_weight(  core::Real const & setting ){ chainbreak_weight_ = setting; }

		bool const & allow_skip_bulge() const { return allow_skip_bulge_; }
		void set_allow_skip_bulge( bool const & setting ){ allow_skip_bulge_ = setting; }

		core::Real const & from_scratch_frequency() const { return from_scratch_frequency_; }
		void set_from_scratch_frequency( core::Real const & setting ){ from_scratch_frequency_ = setting; }

		bool const & allow_split_off() const { return allow_split_off_; }
		void set_allow_split_off( bool const & setting ){ allow_split_off_ = setting; }

		bool const & virtual_sugar_keep_base_fixed() const { return virtual_sugar_keep_base_fixed_; }
		void set_virtual_sugar_keep_base_fixed(  bool const & setting ){ virtual_sugar_keep_base_fixed_ = setting; }

		bool const & virtual_sugar_do_minimize() const { return virtual_sugar_do_minimize_; }
		void set_virtual_sugar_do_minimize(  bool const & setting ){ virtual_sugar_do_minimize_ = setting; }

		bool const & make_movie() const { return make_movie_; }
		void set_make_movie( bool const & setting ){ make_movie_ = setting; }

		utility::vector1< core::Size > const & extra_minimize_res() const { return extra_minimize_res_; }
		void set_extra_minimize_res(  utility::vector1< core::Size > const & setting ){ extra_minimize_res_ = setting; }

		utility::vector1< core::Size > const & syn_chi_res_list() const { return syn_chi_res_list_; }
		void set_syn_chi_res_list(  utility::vector1< core::Size > const & setting ){  syn_chi_res_list_ = setting; }

		utility::vector1< core::Size > const & terminal_res() const { return terminal_res_; }
		void set_terminal_res(  utility::vector1< core::Size > const & setting ){  terminal_res_ = setting; }

		bool const & sampler_perform_phosphate_pack() const { return sampler_perform_phosphate_pack_; }
		void set_sampler_perform_phosphate_pack( bool const & setting ){ sampler_perform_phosphate_pack_ = setting; }

		bool const & rebuild_bulge_mode() const { return rebuild_bulge_mode_; }
		void set_rebuild_bulge_mode( bool const & setting ){ rebuild_bulge_mode_ = setting; }

		bool const & unified_framework() const { return unified_framework_; }
		void set_unified_framework( bool const & setting ){ unified_framework_ = setting; }

		bool const & tether_jump() const { return tether_jump_; }
		void set_tether_jump( bool const & setting ){ tether_jump_ = setting; }

		bool const & local_redock_only() const { return local_redock_only_; }
		void set_local_redock_only( bool const & setting ){ local_redock_only_ = setting; }

		bool const & skip_coord_constraints() const { return skip_coord_constraints_; }
		void set_skip_coord_constraints( bool const & setting ){ skip_coord_constraints_ = setting; }

		bool const & output_minimized_pose_list() const { return output_minimized_pose_list_; }
		void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

		bool const & filter_native_big_bins() const { return filter_native_big_bins_; }
		void set_filter_native_big_bins( bool const & setting ){ filter_native_big_bins_ = setting; }

		bool const & allow_virtual_side_chains() const { return allow_virtual_side_chains_; }
		void set_allow_virtual_side_chains( bool const & setting ){ allow_virtual_side_chains_ = setting; }

		bool const & protein_prepack() const { return protein_prepack_; }
		void set_protein_prepack( bool const & setting ){ protein_prepack_ = setting; }

		bool const & protein_atr_rep_screen() const { return protein_atr_rep_screen_; }
		void set_protein_atr_rep_screen( bool const & setting ){ protein_atr_rep_screen_ = setting; }

		Real const & rmsd_screen() const { return rmsd_screen_; }
		void set_rmsd_screen( Real const & setting ){ rmsd_screen_ = setting; }

	private:

		std::string silent_file_;
		bool verbose_scores_;
		bool integration_test_mode_;
		bool force_centroid_interaction_;
		core::Real sampler_max_centroid_distance_;
		bool use_phenix_geo_;
		bool skip_deletions_;
		bool erraser_;
		bool allow_internal_hinge_moves_;
		bool allow_internal_local_moves_;
		core::Size num_random_samples_;
		core::Size num_pose_minimize_;
		core::Size cycles_;
		core::Real add_delete_frequency_;
		core::Real intermolecular_frequency_;
		core::Real minimize_single_res_frequency_;
		bool minimizer_allow_variable_bond_geometry_;
		core::Real minimizer_vary_bond_geometry_frequency_;
		core::Real switch_focus_frequency_;
		core::Real just_min_after_mutation_frequency_;
		core::Real temperature_;
		utility::vector1< Size > bulge_res_; // disallow addition of these.
		Real max_missing_weight_;
		Real chainbreak_weight_;
		bool allow_skip_bulge_;
		core::Real from_scratch_frequency_;
		bool allow_split_off_;
		bool virtual_sugar_keep_base_fixed_;
		bool virtual_sugar_do_minimize_;
		bool make_movie_;
		bool sampler_perform_phosphate_pack_;
		bool rebuild_bulge_mode_;
		bool unified_framework_;
		bool tether_jump_;
		bool local_redock_only_;
		bool skip_coord_constraints_;
		bool output_minimized_pose_list_;
		bool filter_native_big_bins_;
		bool allow_virtual_side_chains_;
		bool protein_prepack_;
		bool protein_atr_rep_screen_;
		Real rmsd_screen_;

		utility::vector1< core::Size > extra_minimize_res_;
		utility::vector1< core::Size > syn_chi_res_list_;
		utility::vector1< core::Size > terminal_res_;

	};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
