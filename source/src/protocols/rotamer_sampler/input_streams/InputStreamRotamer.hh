// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/input_streams/InputStreamRotamer.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_input_streams_InputStreamRotamer_HH
#define INCLUDED_protocols_rotamer_sampler_input_streams_InputStreamRotamer_HH

#include <protocols/rotamer_sampler/RotamerSized.hh>
#include <protocols/rotamer_sampler/input_streams/InputStreamRotamer.fwd.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/sampling/protein/legacy/StepWiseProteinPoseSetup.hh>

namespace protocols {
namespace rotamer_sampler {
namespace input_streams {

	class InputStreamRotamer: public rotamer_sampler::RotamerSized {

	public:

		//constructor
		InputStreamRotamer( stepwise::sampling::protein::InputStreamWithResidueInfoOP input_stream );

		InputStreamRotamer( stepwise::sampling::protein::InputStreamWithResidueInfoOP input_stream,
												stepwise::sampling::protein::StepWiseProteinPoseSetupCOP stepwise_pose_setup );

		//destructor
		~InputStreamRotamer();

	public:

		/// @brief Reset to the first (or random if random()) rotamer.
		virtual void reset();

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const{ return size_;}

		/// @brief set ID -- how RotamerSizedComb controls RotamerSized. Need some extra work here with InputStreamRotamer.
		virtual void set_id( Size const setting );

		/// @brief Move to next rotamer
		virtual void operator++();

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose&, core::Size const );

		/// @brief Name of the class
		virtual std::string get_name() const { return "InputStreamRotamer"; }

		/// @brief Type of class (see enum in RotamerTypes.hh)
		virtual RotamerType type() const { return INPUT_STREAM; }

	private:

		stepwise::sampling::protein::InputStreamWithResidueInfoOP input_stream_;
		stepwise::sampling::protein::StepWiseProteinPoseSetupCOP stepwise_pose_setup_; // this is really a legacy of SWA protein code.
		Size size_;

	};

} //input_streams
} //rotamer_sampler
} //protocols

#endif
