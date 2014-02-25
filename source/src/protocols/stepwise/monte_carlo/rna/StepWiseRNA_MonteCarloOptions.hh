// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseRNA_MonteCarloOptions_HH
#define INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseRNA_MonteCarloOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	class StepWiseRNA_MonteCarloOptions: public basic::resource_manager::ResourceOptions {

	public:

		//constructor
		StepWiseRNA_MonteCarloOptions();

		//destructor
		~StepWiseRNA_MonteCarloOptions();

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
		type() const{ return "StepWiseRNA_ModelOptions";}

	public:

		void
		initialize_from_command_line();

		protocols::stepwise::sampling::rna::StepWiseRNA_ModelerOptionsOP setup_modeler_options() const;

		bool const & verbose_scores() const { return verbose_scores_; }
		void set_verbose_scores( bool const & setting ){ verbose_scores_ = setting; };

		bool const & force_centroid_interaction() const { return force_centroid_interaction_; }
		void set_force_centroid_interaction( bool const & setting ){ force_centroid_interaction_ = setting; }

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

		core::Size const & cycles() const { return cycles_; }
		void set_cycles( core::Size const & setting ){ cycles_ = setting; }

		core::Real const & add_delete_frequency() const { return add_delete_frequency_; }
		void set_add_delete_frequency( core::Real const & setting ){ add_delete_frequency_ = setting; }

		core::Real const & minimize_single_res_frequency() const { return minimize_single_res_frequency_; }
		void set_minimize_single_res_frequency(  core::Real const & setting ){ minimize_single_res_frequency_ = setting; }

		bool const & minimizer_allow_variable_bond_geometry() const { return minimizer_allow_variable_bond_geometry_; }
		void set_minimizer_allow_variable_bond_geometry( bool const & setting ){ minimizer_allow_variable_bond_geometry_ = setting; }

		Real const & minimizer_vary_bond_geometry_frequency() const { return minimizer_vary_bond_geometry_frequency_; }
		void set_minimizer_vary_bond_geometry_frequency( Real const & setting ){ minimizer_vary_bond_geometry_frequency_ = setting; }

		core::Real const & switch_focus_frequency() const { return switch_focus_frequency_; }
		void set_switch_focus_frequency( core::Real const & setting ){ switch_focus_frequency_ = setting; }

		core::Real const & just_min_after_mutation_frequency() const { return just_min_after_mutation_frequency_; }
		void set_just_min_after_mutation_frequency(  core::Real const & setting ){ just_min_after_mutation_frequency_ = setting; }

		core::Real const & temperature() const { return temperature_; }
		void set_temperature(  core::Real const & setting ){ temperature_ = setting; }

		utility::vector1< Size > const & sample_res() const { return sample_res_; }
		void set_sample_res(  utility::vector1< Size > const & setting ){  sample_res_ = setting; }

		utility::vector1< Size > const & bulge_res() const { return bulge_res_; }
		void set_bulge_res( utility::vector1< Size > const & setting ){ bulge_res_ = setting; }

		Real const & max_missing_weight() const { return max_missing_weight_; }
		void set_max_missing_weight( Real const & setting ){ max_missing_weight_ = setting; }

		Real const & chainbreak_weight() const { return chainbreak_weight_; }
		void set_chainbreak_weight(  Real const & setting ){ chainbreak_weight_ = setting; }

		bool const & allow_skip_bulge() const { return allow_skip_bulge_; }
		void set_allow_skip_bulge( bool const & setting ){ allow_skip_bulge_ = setting; }

		bool const & allow_from_scratch() const { return allow_from_scratch_; }
		void set_allow_from_scratch( bool const & setting ){ allow_from_scratch_ = setting; }

		bool const & allow_split_off() const { return allow_split_off_; }
		void set_allow_split_off( bool const & setting ){ allow_split_off_ = setting; }

		bool const & virtual_sugar_keep_base_fixed() const { return virtual_sugar_keep_base_fixed_; }
		void set_virtual_sugar_keep_base_fixed(  bool const & setting ){ virtual_sugar_keep_base_fixed_ = setting; }

		core::Real const & constraint_x0() const { return constraint_x0_; }
		void set_constraint_x0(  core::Real const & setting ){ constraint_x0_ = setting; }

		core::Real const & constraint_tol() const { return constraint_tol_; }
		void set_constraint_tol( core::Real const & setting ){ constraint_tol_ = setting; }

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

	private:

		bool verbose_scores_;
		bool force_centroid_interaction_;
		bool use_phenix_geo_;
		bool skip_deletions_;
		bool erraser_;
		bool allow_internal_hinge_moves_;
		bool allow_internal_local_moves_;
		core::Size num_random_samples_;
		core::Size cycles_;
		core::Real add_delete_frequency_;
		core::Real minimize_single_res_frequency_;
		bool minimizer_allow_variable_bond_geometry_;
		Real minimizer_vary_bond_geometry_frequency_;
		core::Real switch_focus_frequency_;
		core::Real just_min_after_mutation_frequency_;
		core::Real temperature_;
		utility::vector1< Size > sample_res_;
		utility::vector1< Size > bulge_res_; // disallow addition of these.
		Real max_missing_weight_;
		Real chainbreak_weight_;
		bool allow_skip_bulge_;
		bool allow_from_scratch_;
		bool allow_split_off_;
		bool virtual_sugar_keep_base_fixed_;
		core::Real constraint_x0_;
		core::Real constraint_tol_;
		bool make_movie_;
		bool sampler_perform_phosphate_pack_;
		bool rebuild_bulge_mode_;
		utility::vector1< core::Size > extra_minimize_res_;
		utility::vector1< core::Size > syn_chi_res_list_;
		utility::vector1< core::Size > terminal_res_;

	};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
