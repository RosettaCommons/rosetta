// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinModelerOptions_HH
#define INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinModelerOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class StepWiseProteinModelerOptions: public basic::resource_manager::ResourceOptions {

	public:

		//constructor
		StepWiseProteinModelerOptions();

		StepWiseProteinModelerOptions( StepWiseProteinModelerOptions const & src );

		//destructor
		~StepWiseProteinModelerOptions();

	public:

		StepWiseProteinModelerOptionsOP clone() const;

		StepWiseProteinModelerOptions &
		operator = ( StepWiseProteinModelerOptions const & src );

		/// @brief Describe this instance to a given output stream
		virtual
		void
		show( std::ostream & ) const{}

		/// @brief Initialize from the recursive "tag" structure.
		virtual
		void
		parse_my_tag( utility::tag::TagCOP ){}

		/// @brief The class name (its type) for a particular ResourceOptions instance.
		/// This function allows for better error message delivery.
		virtual
		std::string
		type() const{ return "StepWiseProteinModelerOptions";}

	public:

		void
		initialize_from_command_line();

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }
		std::string silent_file() const{ return silent_file_; }

		void set_global_optimize( bool const & setting ){ global_optimize_ = setting; }
		bool global_optimize() const{ return global_optimize_; }

		void set_mapfile_activated( bool const & setting ){ mapfile_activated_ = setting; }
		bool mapfile_activated() const{ return mapfile_activated_; }

		void set_disallow_backbone_sampling( bool const & setting ){ disallow_backbone_sampling_ = setting; }
		bool disallow_backbone_sampling() const{ return disallow_backbone_sampling_; }

		void set_dump( bool const & setting ){ dump_ = setting; }
		bool dump() const{ return dump_; }

		void set_sample_beta( bool const & setting ){ sample_beta_ = setting; }
		bool sample_beta() const{ return sample_beta_; }

		void set_move_jumps_between_chains( bool const & setting ){ move_jumps_between_chains_ = setting; }
		bool move_jumps_between_chains() const{ return move_jumps_between_chains_; }

		void set_disable_sampling_of_loop_takeoff( bool const & setting ){ disable_sampling_of_loop_takeoff_ = setting; }
		bool disable_sampling_of_loop_takeoff() const{ return disable_sampling_of_loop_takeoff_; }

		void set_rescore_only( bool const & setting ){ rescore_only_ = setting; }
		bool rescore_only() const{ return rescore_only_; }

		void set_cart_min( bool const & setting ){ cart_min_ = setting; }
		bool cart_min() const{ return cart_min_; }

		void set_n_sample( Size const & setting ){ n_sample_ = setting; }
		Size n_sample() const{ return n_sample_; }

		void set_filter_native_big_bins( bool const & setting ){ filter_native_big_bins_ = setting; }
		bool filter_native_big_bins() const{ return filter_native_big_bins_; }

		void set_allow_virtual_side_chains( bool const & setting ){ allow_virtual_side_chains_ = setting; }
		bool allow_virtual_side_chains() const{ return allow_virtual_side_chains_; }

		void set_prepack( bool const & setting ){ prepack_ = setting; }
		bool prepack() const{ return prepack_; }

		void set_atr_rep_screen( bool const & setting ){ atr_rep_screen_ = setting; }
		bool atr_rep_screen() const{ return atr_rep_screen_; }

		void set_centroid_output( bool const & setting ){ centroid_output_ = setting; }
		bool centroid_output() const{ return centroid_output_; }

		void set_centroid_screen( bool const & setting ){ centroid_screen_ = setting; }
		bool centroid_screen() const{ return centroid_screen_; }

		void set_ghost_loops( bool const & setting ){ ghost_loops_ = setting; }
		bool ghost_loops() const{ return ghost_loops_; }

		void set_ccd_close( bool const & setting ){ ccd_close_ = setting; }
		bool ccd_close() const{ return ccd_close_; }

		void set_nstruct_centroid( Size const & setting ){ nstruct_centroid_ = setting; }
		Size nstruct_centroid() const{ return nstruct_centroid_; }

		void set_centroid_weights( std::string const & setting ){ centroid_weights_ = setting; }
		std::string centroid_weights() const{ return centroid_weights_; }

		void set_centroid_score_diff_cut( Real const & setting ){ centroid_score_diff_cut_ = setting; }
		Real centroid_score_diff_cut() const{ return centroid_score_diff_cut_; }

		void set_rmsd_screen( Real const & setting ){ rmsd_screen_ = setting; }
		Real rmsd_screen() const{ return rmsd_screen_; }

		void set_cluster_radius( Real const & setting ){ cluster_radius_ = setting; }
		Real cluster_radius() const{ return cluster_radius_; }

		void set_cluster_by_all_atom_rmsd( bool const & setting ){ cluster_by_all_atom_rmsd_ = setting; }
		bool cluster_by_all_atom_rmsd() const{ return cluster_by_all_atom_rmsd_; }

		void set_pack_weights( std::string const & setting ){ pack_weights_ = setting; }
		std::string pack_weights() const{ return pack_weights_; }

		void set_use_green_packer( bool const & setting ){ use_green_packer_ = setting; }
		bool use_green_packer() const{ return use_green_packer_; }

		void set_use_packer_instead_of_rotamer_trials( bool const & setting ){ use_packer_instead_of_rotamer_trials_ = setting; }
		bool use_packer_instead_of_rotamer_trials() const{ return use_packer_instead_of_rotamer_trials_; }

		void set_min_tolerance( Real const & setting ){ min_tolerance_ = setting; }
		Real min_tolerance() const{ return min_tolerance_; }

		void set_min_type( std::string const & setting ){ min_type_ = setting; }
		std::string min_type() const{ return min_type_; }

		void set_skip_minimize( bool const & setting ){ skip_minimize_ = setting; }
		bool skip_minimize() const{ return skip_minimize_; }

		bool const & choose_random() const { return choose_random_; }
		void set_choose_random( bool const & setting ){ choose_random_ = setting; }

		Size const & num_random_samples() const { return num_random_samples_; };
		void set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }

		Size const & num_pose_minimize() const { return num_pose_minimize_; };
		void set_num_pose_minimize( Size const & setting ){ num_pose_minimize_ = setting; }

		void set_max_decoys( Size const & setting ){ max_decoys_ = setting; }
		Size max_decoys() const{ return max_decoys_; }

		bool const & expand_loop_takeoff() const { return expand_loop_takeoff_; }
		void set_expand_loop_takeoff( bool const & setting ){ expand_loop_takeoff_ = setting; }

		bool const & skip_coord_constraints() const { return skip_coord_constraints_; }
		void set_skip_coord_constraints( bool const & setting ){ skip_coord_constraints_ = setting; }

		bool const & output_minimized_pose_list() const { return output_minimized_pose_list_; }
		void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

	private:
		void
		initialize_variables();

	private:

		std::string silent_file_;
		bool global_optimize_;
		bool mapfile_activated_;
		bool disallow_backbone_sampling_; //a.k.a., skip_sampling
		bool dump_;
		bool sample_beta_;
		bool move_jumps_between_chains_;
		bool disable_sampling_of_loop_takeoff_;
		bool rescore_only_;
		bool cart_min_;
		Size n_sample_;
		Size max_decoys_;

		bool filter_native_big_bins_;
		bool allow_virtual_side_chains_;
		bool prepack_;
		bool atr_rep_screen_;

		bool centroid_output_;
		bool centroid_screen_;
		Real centroid_score_diff_cut_;
		std::string centroid_weights_;
		Size nstruct_centroid_;
		bool ghost_loops_;

		bool ccd_close_; // should this be in here, or in job parameters?

		Real rmsd_screen_;
		Real cluster_radius_;
		bool cluster_by_all_atom_rmsd_;

		std::string pack_weights_;
		bool use_green_packer_;
		bool use_packer_instead_of_rotamer_trials_;

		std::string min_type_;
		Real min_tolerance_;
		bool skip_minimize_;
		bool choose_random_;
		Size num_random_samples_;
		Size num_pose_minimize_;
		bool expand_loop_takeoff_;
		bool skip_coord_constraints_;
		bool output_minimized_pose_list_;

	};

} //protein
} //sampling
} //stepwise
} //protocols

#endif
