// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/input_streams/InputStreamStepWiseSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_input_streams_InputStreamStepWiseSampler_HH
#define INCLUDED_protocols_sampler_input_streams_InputStreamStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/input_streams/InputStreamStepWiseSampler.fwd.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinPoseSetup.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace input_streams {

	class InputStreamStepWiseSampler: public sampler::StepWiseSamplerSized {

	public:

		//constructor
		InputStreamStepWiseSampler( stepwise::modeler::protein::InputStreamWithResidueInfoOP input_stream );

		InputStreamStepWiseSampler( stepwise::modeler::protein::InputStreamWithResidueInfoOP input_stream,
												stepwise::legacy::modeler::protein::StepWiseProteinPoseSetupCOP stepwise_pose_setup );

		//destructor
		~InputStreamStepWiseSampler();

	public:

		/// @brief Reset to the first (or random if random()) rotamer.
		virtual void reset();

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const{ return size_;}

		/// @brief set ID -- how StepWiseSamplerSizedComb controls StepWiseSamplerSized. Need some extra work here with InputStreamStepWiseSampler.
		virtual void set_id( Size const setting );

		/// @brief Move to next rotamer
		virtual void operator++();

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose&, core::Size const );

		/// @brief Name of the class
		virtual std::string get_name() const { return "InputStreamStepWiseSampler"; }

		/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
		virtual StepWiseSamplerType type() const { return INPUT_STREAM; }

	private:

		stepwise::modeler::protein::InputStreamWithResidueInfoOP input_stream_;
		stepwise::legacy::modeler::protein::StepWiseProteinPoseSetupCOP stepwise_pose_setup_; // this is really a legacy of SWA protein code.
		Size size_;

	};

} //input_streams
} //sampler
} //stepwise
} //protocols

#endif
