// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/input_streams/InputStreamStepWiseSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/input_streams/InputStreamStepWiseSampler.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinPoseSetup.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sampler.input_streams.InputStreamStepWiseSampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace input_streams {

	//Constructor
	InputStreamStepWiseSampler::InputStreamStepWiseSampler( stepwise::modeler::protein::InputStreamWithResidueInfoOP input_stream ):
		input_stream_( input_stream ),
		size_( input_stream->compute_size() ) // slow...
  {
		set_random( false );
		init();
  }

	//Constructor
	InputStreamStepWiseSampler::InputStreamStepWiseSampler( stepwise::modeler::protein::InputStreamWithResidueInfoOP input_stream,
																					stepwise::legacy::modeler::protein::StepWiseProteinPoseSetupCOP stepwise_pose_setup ):
		input_stream_( input_stream ),
		stepwise_pose_setup_( stepwise_pose_setup ),
		size_( input_stream->compute_size() ) // slow...
  {
		set_random( false );
		init();
  }

	//Destructor
	InputStreamStepWiseSampler::~InputStreamStepWiseSampler()
	{}

  //////////////////////////////////////////////////////////////////////////
  void
  InputStreamStepWiseSampler::reset(){
		StepWiseSamplerSized::reset(); // rewinds id_ to 1.
		input_stream_->reset();
		input_stream_->advance_to_next_pose_segment(); // need to get first pose set up.
	}

  //////////////////////////////////////////////////////////////////////////
	void
	InputStreamStepWiseSampler::set_id( Size const setting ){
		if ( setting == id_ ) return;
		if ( setting == id_ + 1 ) {
			++( *this );
			return;
		}
		if ( setting == 1 ) { // reset
			reset();
			return;
		}
		runtime_assert( setting == 0 );
		init();
	}

	///////////////////////////////////////////////////////////////////////////
	void
	InputStreamStepWiseSampler::operator++() {
		runtime_assert ( !random() );
		if ( id_ < size() ){
			runtime_assert( input_stream_->has_another_pose() );
			input_stream_->advance_to_next_pose_segment();
		}
		++id_;
	}

	/// @brief Apply the i-th rotamer to pose
	void
	InputStreamStepWiseSampler::apply( core::pose::Pose & pose, core::Size const id ){
		runtime_assert( id == id_ ); // cannot currently request an arbitrary pose from the stream.
		runtime_assert( !random() );

		input_stream_->apply_current_pose_segment( pose );
		// This is annoying but necessary:
		if ( ( stepwise_pose_setup_ != 0 ) && stepwise_pose_setup_->ready_to_align() ){
			//			TR << "ALIGNING POSE! " << std::endl;
			stepwise_pose_setup_->align_pose( pose );
		}

	}

} //input_streams
} //sampler
} //stepwise
} //protocols
