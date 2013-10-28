// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/RNA_StepWiseMonteCarlo.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_monte_carlo_RNA_StepWiseMonteCarlo_HH
#define INCLUDED_protocols_swa_monte_carlo_RNA_StepWiseMonteCarlo_HH

#include <protocols/moves/Mover.hh>
#include <protocols/swa/monte_carlo/RNA_StepWiseMonteCarlo.fwd.hh>
#include <protocols/swa/monte_carlo/RNA_AddMover.fwd.hh>
#include <protocols/swa/monte_carlo/RNA_DeleteMover.fwd.hh>
#include <protocols/swa/monte_carlo/RNA_AddOrDeleteMover.fwd.hh>
#include <protocols/swa/monte_carlo/RNA_ResampleMover.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace monte_carlo {

	class RNA_StepWiseMonteCarlo: public protocols::moves::Mover {

	public:

	//constructor
		RNA_StepWiseMonteCarlo( core::scoring::ScoreFunctionOP scorefxn );

		//destructor
		~RNA_StepWiseMonteCarlo();

		virtual std::string get_name() const {
			return "RNA_StepWiseMonteCarlo";
		}

		/// @brief Apply the loop-rebuild protocol to the input pose
		virtual
		void apply ( core::pose::Pose & pose );

	public:

		void set_verbose_scores( bool const & setting ){ verbose_scores_ = setting; }
		bool verbose_scores() const{ return verbose_scores_; }

		void set_use_phenix_geo( bool const & setting ){ use_phenix_geo_ = setting; }
		bool use_phenix_geo() const{ return use_phenix_geo_; }

		void set_skip_deletions( bool const & setting ){ skip_deletions_ = setting; }
		bool skip_deletions() const{ return skip_deletions_; }

		void set_erraser( bool const & setting ){ erraser_ = setting; }
		bool erraser() const{ return erraser_; }

		void set_allow_internal_moves( bool const & setting ){ allow_internal_moves_ = setting; }
		bool allow_internal_moves() const{ return allow_internal_moves_; }

		void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }
		core::Size num_random_samples() const{ return num_random_samples_; }

		void set_cycles( core::Size const & setting ){ cycles_ = setting; }
		core::Size cycles() const{ return cycles_; }

		void set_add_delete_frequency( core::Real const & setting ){ add_delete_frequency_ = setting; }
		core::Real add_delete_frequency() const{ return add_delete_frequency_; }

		void set_minimize_single_res_frequency( core::Real const & setting ){ minimize_single_res_frequency_ = setting; }
		core::Real minimize_single_res_frequency() const{ return minimize_single_res_frequency_; }

		void set_switch_focus_frequency( core::Real const & setting ){ switch_focus_frequency_ = setting; }
		core::Real switch_focus_frequency() const{ return switch_focus_frequency_; }

		void set_just_min_after_mutation_frequency( core::Real const & setting ){ just_min_after_mutation_frequency_ = setting; }
		core::Real just_min_after_mutation_frequency() const{ return just_min_after_mutation_frequency_; }

		void set_temperature( core::Real const & setting ){ temperature_ = setting; }
		core::Real temperature() const{ return temperature_; }

		void set_sample_res( utility::vector1<Size> const & setting ){ sample_res_ = setting; }
		utility::vector1<Size> sample_res() const{ return sample_res_; }


	private:

		void initialize_movers();

		void show_scores( core::pose::Pose & pose, std::string const tag );

		bool
		switch_focus_among_poses_randomly( pose::Pose & pose ) const;

	private:
		core::scoring::ScoreFunctionOP scorefxn_;

		bool verbose_scores_;
		bool use_phenix_geo_;
		bool skip_deletions_;
		bool erraser_;
		bool allow_internal_moves_;
		core::Size num_random_samples_;
		core::Size cycles_;
		core::Real add_delete_frequency_;
		core::Real minimize_single_res_frequency_;
		core::Real switch_focus_frequency_;
		core::Real just_min_after_mutation_frequency_;
		core::Real temperature_;
		utility::vector1<Size> sample_res_;

		RNA_DeleteMoverOP rna_delete_mover_;
		RNA_AddMoverOP rna_add_mover_;
		RNA_AddOrDeleteMoverOP rna_add_or_delete_mover_;
		RNA_ResampleMoverOP rna_resample_mover_;
		Real rmsd_weight_;
		Real max_missing_weight_;

	};

} //monte_carlo
} //swa
} //protocols

#endif
