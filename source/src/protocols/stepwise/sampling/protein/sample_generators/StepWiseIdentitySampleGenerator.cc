// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseIdentitySampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/sample_generators/StepWiseIdentitySampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/StepWisePoseSetup.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>
#include <core/types.hh>
#include <utility/exit.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace sample_generators {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseIdentitySampleGenerator::StepWiseIdentitySampleGenerator( InputStreamWithResidueInfoOP input_stream ):
		input_stream_( input_stream )
  {
		reset();
  }

	StepWiseIdentitySampleGenerator::StepWiseIdentitySampleGenerator(
																																	 InputStreamWithResidueInfoOP input_stream,
																																	 StepWisePoseSetupOP stepwise_pose_setup ):
		input_stream_( input_stream ),
		stepwise_pose_setup_( stepwise_pose_setup )
  {
		reset();
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseIdentitySampleGenerator::reset(){
		input_stream_->reset();
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWiseIdentitySampleGenerator::has_another_sample(){
		return ( input_stream_->has_another_pose()  );
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseIdentitySampleGenerator::get_next_sample( core::pose::Pose & pose )
	{

		if ( !input_stream_->has_another_pose() ) utility_exit_with_message( "Asked IdentitySampleGenerator for another sample, but it did not have one!" );
		input_stream_->copy_next_pose_segment( pose );

		// This is annoying but necessary:
		if ( stepwise_pose_setup_ && stepwise_pose_setup_->ready_to_align() ){
			//			std::cout << "ALIGNING POSE! " << std::endl;
			stepwise_pose_setup_->align_pose( pose );
		}

	}

	///////////////////////////////////////////////////////
	Size
	StepWiseIdentitySampleGenerator::size() const{
		//we don't know a priori how many combinations there will be...
		// yes, I guess we could precompute this by looking through both streams ahead of time.
		return 0;
	}


	///////////////////////////////////////////////////////
	void
	StepWiseIdentitySampleGenerator::set_stepwise_pose_setup( StepWisePoseSetupOP stepwise_pose_setup ){
		stepwise_pose_setup_ = stepwise_pose_setup;
	}

} //sample_generators
} //protein
} //sampling
} //stepwise
} //protocols
