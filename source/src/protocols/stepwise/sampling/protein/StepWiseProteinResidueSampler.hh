// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/StepWiseProteinResidueSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinResidueSampler_HH
#define INCLUDED_protocols_stepwise_sampling_protein_StepWiseProteinResidueSampler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinResidueSampler.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinJobParameters.fwd.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.fwd.hh>
#include <protocols/stepwise/sampling/general/StepWiseClusterer.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/stepwise/screener/SimplePoseSelection.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#if defined(WIN32) || defined(PYROSETTA)
	#include <core/scoring/ScoreFunction.hh>
#endif

/*
using namespace core;

Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using
*/

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class StepWiseProteinResidueSampler: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseProteinResidueSampler( StepWiseProteinJobParametersCOP & job_parameters );

		//destructor
		~StepWiseProteinResidueSampler();

	public:

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const { return "StepWiseProteinResidueSampler"; }

		void
		set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){ scorefxn_ = scorefxn; }

		utility::vector1< core::pose::PoseOP > const & get_pose_list(){ return pose_list_; }

		void
		set_options( StepWiseProteinModelerOptionsCOP options ){ options_ = options; }

		void
		set_frag_files( utility::vector1< std::string > const &  frag_files ){ frag_files_ = frag_files; }

		void
		set_input_streams( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams ){ input_streams_ = input_streams; }

		void
		set_job_parameters( StepWiseProteinJobParametersCOP job_parameters );

		void
		set_do_ccd( bool const setting ){ do_ccd_ = setting; }

		void set_moving_res_list( utility::vector1< Size > const & setting ){ moving_res_list_ = setting; }
		utility::vector1< Size > moving_res_list() const { return moving_res_list_; }

		bool full_optimize() const { return full_optimize_; }

	private:

		void
		initialize_sampler( core::pose::Pose & pose );

		void
		initialize_screeners( core::pose::Pose & pose );

		rotamer_sampler::RotamerSizedOP
		get_basic_sampler( core::pose::Pose & pose );

		rotamer_sampler::RotamerSizedOP
		close_loops_in_samples( core::pose::Pose & pose, rotamer_sampler::RotamerSizedOP sampler );

		StepWiseProteinPackerOP
		get_packer( core::pose::Pose & pose );

		void
		enable_sampling_of_loop_takeoff( core::pose::Pose & pose,
																		 rotamer_sampler::RotamerSizedOP sampler );

	private:

		StepWiseProteinJobParametersCOP job_parameters_;
		utility::vector1< core::pose::PoseOP > pose_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		StepWiseProteinModelerOptionsCOP options_;
		utility::vector1< Size > moving_res_list_;

		core::pose::PoseOP screening_pose_;

		rotamer_sampler::RotamerSizedOP sampler_;
		utility::vector1< screener::StepWiseScreenerOP > screeners_;
		screener::SimplePoseSelectionOP pose_selection_;

		// updated and passed between sampler setup, loop closing, etc.
		bool do_ccd_;
		bool close_loops_;
		bool full_optimize_;

		// setup in swa_protein_main, but not implemented/necessary for stepwise monte carlo.
		utility::vector1< std::string > frag_files_; // consider moving into full_model_info_parameters?
		utility::vector1< InputStreamWithResidueInfoOP > input_streams_; // this is now awkward.


	};

} //protein
} //sampling
} //stepwise
} //protocols

#endif
