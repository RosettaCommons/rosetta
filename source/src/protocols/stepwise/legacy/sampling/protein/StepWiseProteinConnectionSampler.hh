// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/sampling/protein/StepWiseProteinConnectionSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinConnectionSampler_HH
#define INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinConnectionSampler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/legacy/sampling/protein/StepWiseProteinConnectionSampler.fwd.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.fwd.hh>
#include <protocols/stepwise/sampling/packer/StepWisePacker.fwd.hh>
#include <protocols/stepwise/sampling/align/StepWiseLegacyClusterer.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/stepwise/screener/SimplePoseSelection.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#if defined(WIN32) || defined(PYROSETTA)
	#include <core/scoring/ScoreFunction.hh>
	#include <protocols/stepwise/screener/StepWiseScreener.hh>
#endif

/*
using namespace core;

Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using
*/

namespace protocols {
namespace stepwise {
namespace legacy {
namespace sampling {
namespace protein {

	class StepWiseProteinConnectionSampler: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseProteinConnectionSampler( stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP & working_parameters );

		//destructor
		~StepWiseProteinConnectionSampler();

	public:

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const { return "StepWiseProteinConnectionSampler"; }

		void
		set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ){ scorefxn_ = scorefxn; }

		utility::vector1< core::pose::PoseOP > const & get_pose_list(){ return pose_list_; }

		void
		set_options( modeler_options::StepWiseModelerOptionsCOP options ){ options_ = options; }

		void set_skip_sampling( bool const & setting ){ skip_sampling_ = setting; }

		void
		set_input_streams( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams ){ input_streams_ = input_streams; }

		void
		set_working_parameters( stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP working_parameters );

		utility::vector1< Size > const & working_obligate_pack_res() const { return working_obligate_pack_res_; }


	private:

		void
		initialize_poses_and_checkers( core::pose::Pose & pose );

		void
		initialize_sampler( core::pose::Pose & pose );

		void
		initialize_screeners( core::pose::Pose & pose );

		void
		enable_sampling_of_loop_takeoff( core::pose::Pose & pose,
																		 rotamer_sampler::RotamerSizedOP sampler );

	private:

		stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP working_parameters_;
		utility::vector1< Size > const moving_res_list_;
		utility::vector1< Size > working_obligate_pack_res_;
		utility::vector1< core::pose::PoseOP > pose_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		modeler_options::StepWiseModelerOptionsCOP options_;

		utility::vector1< core::Size > protein_cutpoints_closed_;
		utility::vector1< core::pose::PoseOP > ccd_poses_;
		core::pose::PoseOP atr_rep_screening_pose_;

		rotamer_sampler::RotamerSizedOP sampler_;
		utility::vector1< screener::StepWiseScreenerOP > screeners_;
		utility::vector1< loop_close::StepWiseProteinCCD_CloserOP > stepwise_ccd_closers_;
		packer::StepWisePackerOP stepwise_packer_;
		screener::SimplePoseSelectionOP pose_selection_;

		// updated and passed between sampler setup, loop closing, etc.
		bool pack_all_side_chains_;
		bool skip_sampling_;

		// setup in swa_protein_main, but not implemented/necessary for stepwise monte carlo.
		utility::vector1< InputStreamWithResidueInfoOP > input_streams_; // this is now awkward.


	};

} //protein
} //sampling
} //legacy
} //stepwise
} //protocols

#endif
