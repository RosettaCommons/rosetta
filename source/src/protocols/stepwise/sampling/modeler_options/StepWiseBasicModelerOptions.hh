// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/modeler_options/StepWiseBasicModelerOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_StepWiseBasicModelerOptions_HH
#define INCLUDED_protocols_stepwise_sampling_StepWiseBasicModelerOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseBasicModelerOptions.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace sampling {

	class StepWiseBasicModelerOptions: public virtual basic::resource_manager::ResourceOptions {

	public:

		//constructor
		StepWiseBasicModelerOptions();

		StepWiseBasicModelerOptions( StepWiseBasicModelerOptions const & src );

		//destructor
		~StepWiseBasicModelerOptions();

	public:

		StepWiseBasicModelerOptionsOP clone() const;

		StepWiseBasicModelerOptions &
		operator = ( StepWiseBasicModelerOptions const & src );

		/// @brief Initialize from the recursive "tag" structure.
		virtual
		void
		parse_my_tag( utility::tag::TagCOP ){}

		/// @brief The class name (its type) for a particular ResourceOptions instance.
		/// This function allows for better error message delivery.
		virtual
		std::string
		type() const{ return "StepWiseBasicModelerOptions";}

		std::string const & silent_file() const { return silent_file_; }
		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		core::Size const & sampler_num_pose_kept() const { return sampler_num_pose_kept_; }
		void set_sampler_num_pose_kept( core::Size const & setting ){ sampler_num_pose_kept_ = setting; }

		core::Size const & num_pose_minimize() const { return num_pose_minimize_; }
		void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

		bool const & use_green_packer() const { return use_green_packer_; }
		void set_use_green_packer( bool const & setting ){ use_green_packer_ = setting; }

		bool const & verbose() const { return verbose_; }
		void set_verbose( bool const & setting ){ verbose_ = setting; }

		bool const & choose_random() const { return choose_random_; }
		void set_choose_random( bool const & setting ){ choose_random_ = setting; }

		core::Size const & num_random_samples() const { return num_random_samples_; };
		void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }

		void set_dump( bool const & setting ){ dump_ = setting; }
		bool dump() const{ return dump_; }

		void set_atr_rep_screen( bool const & setting ){ atr_rep_screen_ = setting; }
		bool atr_rep_screen() const{ return atr_rep_screen_; }

		void set_rmsd_screen( core::Real const & setting ){ rmsd_screen_ = setting; }
		core::Real rmsd_screen() const{ return rmsd_screen_; }

		core::Real const & cluster_rmsd() const { return cluster_rmsd_; }
		void set_cluster_rmsd( core::Real const & setting ){ cluster_rmsd_ = setting; }

		void set_skip_minimize( bool const & setting ){ skip_minimize_ = setting; }
		bool skip_minimize() const{ return skip_minimize_; }

		bool const & output_minimized_pose_list() const { return output_minimized_pose_list_; }
		void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

		bool const & disallow_realign() const { return disallow_realign_; }
		void set_disallow_realign( bool const & setting ){ disallow_realign_ = setting; }

		bool const & unified_sampler() const { return unified_sampler_; }
		void set_unified_sampler( bool const & setting ){ unified_sampler_ = setting; }

protected:

		void
		initialize_from_command_line();

    void
		initialize_variables();

protected:

		std::string silent_file_;
		core::Size sampler_num_pose_kept_;
		core::Size num_pose_minimize_;
		bool use_green_packer_;
		bool verbose_;
		bool choose_random_;
		core::Size num_random_samples_;
		bool dump_;
		bool atr_rep_screen_;
		core::Real rmsd_screen_;
		core::Real cluster_rmsd_;
		bool skip_minimize_;
		bool output_minimized_pose_list_;
		bool disallow_realign_;
		bool unified_sampler_;

	};

} //sampling
} //stepwise
} //protocols

#endif
