// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/StepWiseProteinMultiPosePacker.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinMultiPosePacker_HH
#define INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinMultiPosePacker_HH

#include <protocols/stepwise/sampling/protein/StepWiseProteinMultiPosePacker.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.fwd.hh>
#include <protocols/stepwise/sampling/protein/checker/ProteinAtrRepChecker.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class StepWiseProteinMultiPosePacker: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseProteinMultiPosePacker( StepWiseProteinPackerOP packer,
																		rotamer_sampler::RotamerSizedOP sampler );

		//destructor
		~StepWiseProteinMultiPosePacker();

	public:

    /// @brief Apply the minimizer to one pose
    virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		void
		set_screener( protocols::stepwise::screener::StepWiseScreenerOP pose_screener ){ pose_screener_ = pose_screener; }

		void
		set_prepack_based_on_moving_res_list( utility::vector1< Size > const & moving_res_list );

		utility::vector1< pose::PoseOP > pose_list() { return pose_list_; }

		void
		sample_residues( pose::Pose & pose );

		void set_atr_rep_check( bool const & setting ){ atr_rep_check_ = setting; }

		void set_choose_random( bool const & setting ){ choose_random_ = setting; }
		bool choose_random() const{ return choose_random_; }

		void set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }
		Size num_random_samples() const{ return num_random_samples_; }

	private:

		void
		initialize_atr_rep_checker( pose::Pose const & pose );

	private:

		utility::vector1< pose::PoseOP >  pose_list_;

		StepWiseProteinPackerOP packer_;
		rotamer_sampler::RotamerSizedOP sampler_;

		screener::StepWiseScreenerOP pose_screener_;
		checker::ProteinAtrRepCheckerOP atr_rep_checker_;
		utility::vector1< Size > moving_res_list_;
		bool prepack_;
		bool atr_rep_check_;
		bool choose_random_;
		Size num_random_samples_;
		Size max_ntries_;

	};

} //protein
} //sampling
} //stepwise
} //protocols

#endif
