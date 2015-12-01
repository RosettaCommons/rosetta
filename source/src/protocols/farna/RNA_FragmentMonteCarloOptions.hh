// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/RNA_FragmentMonteCarloOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_FragmentMonteCarloOptions_HH
#define INCLUDED_protocols_farna_RNA_FragmentMonteCarloOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/farna/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace farna {

	class RNA_FragmentMonteCarloOptions: public virtual basic::resource_manager::ResourceOptions {

	public:

		//constructor
		RNA_FragmentMonteCarloOptions();

		RNA_FragmentMonteCarloOptions( RNA_FragmentMonteCarloOptions const & src );

		//destructor
		~RNA_FragmentMonteCarloOptions();

	public:

		RNA_FragmentMonteCarloOptionsOP clone() const;

		void
		initialize_from_command_line();

		/// @brief Initialize from the recursive "tag" structure.
		virtual
		void
		parse_my_tag( utility::tag::TagCOP ){}

		/// @brief The class name (its type) for a particular ResourceOptions instance.
		/// This function allows for better error message delivery.
		virtual
		std::string
		type() const{ return "RNA_FragmentMonteCarloOptions";}

		void set_rounds( core::Size const & setting ){ rounds_ = setting; }
		core::Size rounds() const { return rounds_; }

		void
		set_dump_pdb( bool const setting ){ dump_pdb_ = setting; };
		bool dump_pdb() const { return dump_pdb_; }

		void set_close_loops_after_each_move( bool const & setting ){ close_loops_after_each_move_ = setting; }
		bool close_loops_after_each_move() const { return close_loops_after_each_move_; }

		void set_refine_pose( bool const & setting ){ refine_pose_ = setting; }
		bool refine_pose() const { return refine_pose_; }

		void set_temperature( core::Real const setting ){ temperature_ = setting; };
		core::Real temperature() const { return temperature_; }

		void set_autofilter( bool const & setting ){ autofilter_ = setting; }
		bool autofilter() const { return autofilter_; }

		void set_autofilter_score_quantile( core::Real const & setting ){ autofilter_score_quantile_ = setting; }
		core::Real autofilter_score_quantile() const { return autofilter_score_quantile_; }

		void set_minimize_structure( bool const & setting ){ minimize_structure_ = setting; }
		bool minimize_structure() const { return minimize_structure_; }

		void set_relax_structure( bool const & setting ){ relax_structure_ = setting; }
		bool relax_structure() const { return relax_structure_; }

		//void set_use_chem_shift_data( bool const & setting ){ use_chem_shift_data_ = setting; }
		bool use_chem_shift_data() const { return use_chem_shift_data_; }

		void set_jump_change_frequency( core::Real const & setting ){ jump_change_frequency_ = setting; }
		core::Real jump_change_frequency() const { return jump_change_frequency_; }

		void set_close_loops( bool const setting ){	close_loops_ = setting; }
		bool close_loops() const { return close_loops_; }

		void set_allow_bulge( bool const & setting ){ allow_bulge_ = setting; }
		bool allow_bulge() const { return allow_bulge_; }

		void set_allow_consecutive_bulges( bool const & setting ){ allow_consecutive_bulges_ = setting; }
		bool allow_consecutive_bulges() const { return allow_consecutive_bulges_; }

		void set_simple_rmsd_cutoff_relax( bool const & setting ){ simple_rmsd_cutoff_relax_ = setting; }
		bool simple_rmsd_cutoff_relax() const { return simple_rmsd_cutoff_relax_; }

		void
		set_filter_lores_base_pairs( bool const setting ){ filter_lores_base_pairs_ = setting; }
		bool filter_lores_base_pairs() const { return filter_lores_base_pairs_; }

		void
		set_filter_lores_base_pairs_early( bool const setting ){ filter_lores_base_pairs_early_ = setting;	}
		bool filter_lores_base_pairs_early() const { return filter_lores_base_pairs_early_; }

		void
		set_filter_chain_closure( bool const setting ){ filter_chain_closure_ = setting; }
		bool filter_chain_closure() const { return filter_chain_closure_; }

		void
		set_filter_chain_closure_distance( core::Real const setting ){ filter_chain_closure_distance_ = setting; }
		core::Real filter_chain_closure_distance() const { return filter_chain_closure_distance_; }

		void
		set_filter_chain_closure_halfway( bool const setting ){ filter_chain_closure_halfway_ = setting; }
		bool filter_chain_closure_halfway() const { return filter_chain_closure_halfway_; }


		void set_chainbreak_weight( core::Real const & setting ){ chainbreak_weight_ = setting; }
		core::Real chainbreak_weight() const { return chainbreak_weight_; }

		void set_linear_chainbreak_weight( core::Real const & setting ){ linear_chainbreak_weight_ = setting; }
		core::Real linear_chainbreak_weight() const { return linear_chainbreak_weight_; }

		void set_allowed_bulge_res( utility::vector1< core::Size > const & setting ){ allowed_bulge_res_ = setting; };
		utility::vector1< core::Size > const & allowed_bulge_res() const { return  allowed_bulge_res_; };

		void set_minimizer_use_coordinate_constraints( bool const & setting ){ minimizer_use_coordinate_constraints_ = setting; }
		bool minimizer_use_coordinate_constraints() const { return minimizer_use_coordinate_constraints_; }

		void set_monte_carlo_cycles( core::Size const setting ){ monte_carlo_cycles_ = setting; }
		core::Size monte_carlo_cycles() const { return monte_carlo_cycles_; }

		void set_user_defined_cycles( bool const & setting ){ user_defined_cycles_ = setting; }
		bool user_defined_cycles() const { return user_defined_cycles_; }

		void set_ignore_secstruct( bool const & setting ){ ignore_secstruct_ = setting; }
		bool ignore_secstruct() const { return ignore_secstruct_; }

		void set_vary_bond_geometry( bool const & setting ){ vary_bond_geometry_ = setting; }
		bool vary_bond_geometry() const { return vary_bond_geometry_; }

		void set_titrate_stack_bonus( bool const & setting ){ titrate_stack_bonus_ = setting; }
		bool titrate_stack_bonus() const { return titrate_stack_bonus_; }

		void set_staged_constraints( bool const & setting ){ staged_constraints_ = setting; }
		bool staged_constraints() const { return staged_constraints_; }

		void set_move_first_rigid_body( bool const & setting ){ move_first_rigid_body_ = setting; }
		bool move_first_rigid_body() const { return move_first_rigid_body_; }

		void set_root_at_first_rigid_body( bool const & setting ){ root_at_first_rigid_body_ = setting; }
		bool root_at_first_rigid_body() const { return root_at_first_rigid_body_; }

		void set_suppress_bp_constraint( core::Real const & setting ){ suppress_bp_constraint_ = setting; }
		core::Real suppress_bp_constraint() const { return suppress_bp_constraint_; }

		void set_bps_moves( bool const & setting ){ bps_moves_ = setting; }
		bool bps_moves() const { return bps_moves_; }

		void set_extra_minimize_res( utility::vector1< Size > const & extra_minimize_res ) {	extra_minimize_res_ = extra_minimize_res;	}
		utility::vector1< Size > const & extra_minimize_res() const { return extra_minimize_res_; };

		void set_extra_minimize_chi_res( utility::vector1< Size > const & extra_minimize_chi_res ) {	extra_minimize_chi_res_ = extra_minimize_chi_res;	}
		utility::vector1< Size > const & extra_minimize_chi_res() const { return extra_minimize_chi_res_; };

		void set_jump_library_file( std::string const & setting ){ jump_library_file_ = setting; }
		std::string jump_library_file() const { return jump_library_file_; }

		void set_all_rna_fragments_file( std::string const & setting ){ all_rna_fragments_file_ = setting; }
		std::string all_rna_fragments_file() const { return all_rna_fragments_file_; }

		void set_vall_torsions_file( std::string const & setting ){ all_rna_fragments_file_ = setting; }
		std::string vall_torsions_file() const { return all_rna_fragments_file_; }

		void set_rna_params_file( std::string const & setting ){ rna_params_file_ = setting; }
		std::string rna_params_file() const { return rna_params_file_; }

		void
		set_chunk_pdb_files( utility::vector1< std::string > const & chunk_pdb_files ) {
			chunk_pdb_files_ = chunk_pdb_files;
		}
		utility::vector1< std::string > const & chunk_pdb_files() const { return chunk_pdb_files_; };

		void
		set_chunk_silent_files( utility::vector1< std::string > const & chunk_silent_files ) {
			chunk_silent_files_ = chunk_silent_files;
		}
		utility::vector1< std::string > const & chunk_silent_files() const { return chunk_silent_files_; };

		void set_input_res( utility::vector1< Size > const & input_res ) {	input_res_ = input_res;	}
		utility::vector1< Size > const & input_res() const { return input_res_; };

		void set_verbose( bool const & setting ){ verbose_ = setting; }
		bool verbose() const { return verbose_; }

		void set_minimize_bps( bool const & setting ){ minimize_bps_ = setting; }
		bool minimize_bps() const { return minimize_bps_; }

	private:

		Size rounds_;
		Size monte_carlo_cycles_;
		bool user_defined_cycles_;
		bool dump_pdb_;
		bool minimize_structure_;
		bool relax_structure_;

		core::Real temperature_; // default temperature for monte carlo

		bool ignore_secstruct_;

		core::Real chainbreak_weight_;
		core::Real linear_chainbreak_weight_;
		bool close_loops_;
		bool close_loops_in_last_round_;
		bool close_loops_after_each_move_;

		bool allow_bulge_, allow_consecutive_bulges_;

		core::Real jump_change_frequency_;

		bool autofilter_;
		core::Real autofilter_score_quantile_;
		bool titrate_stack_bonus_;
		bool move_first_rigid_body_;
		bool root_at_first_rigid_body_;
		core::Real suppress_bp_constraint_;

		bool filter_lores_base_pairs_;
		bool filter_lores_base_pairs_early_;

		bool filter_chain_closure_;
		core::Real filter_chain_closure_distance_;
		bool filter_chain_closure_halfway_;

		bool staged_constraints_;

		utility::vector1< core::Size > allowed_bulge_res_;
		utility::vector1< core::Size > extra_minimize_res_;
		utility::vector1< core::Size > extra_minimize_chi_res_;

		bool vary_bond_geometry_;
		bool simple_rmsd_cutoff_relax_;

		bool refine_from_silent_;
		bool refine_pose_;
		bool bps_moves_;
		bool minimizer_use_coordinate_constraints_;
		bool minimize_bps_;
		bool use_chem_shift_data_;
		bool verbose_;

		std::string all_rna_fragments_file_;
		std::string rna_params_file_;
		std::string jump_library_file_;

		utility::vector1< std::string > chunk_pdb_files_;
		utility::vector1< std::string > chunk_silent_files_;
		utility::vector1< core::Size > input_res_;

	};

} //farna
} //protocols

#endif
