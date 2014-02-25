// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWisePoseCombineSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/sample_generators/StepWisePoseCombineSampleGenerator.hh>
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
	StepWisePoseCombineSampleGenerator::StepWisePoseCombineSampleGenerator(
																																				 utility::vector1< InputStreamWithResidueInfoOP > const & input_streams ):
		input_streams_( input_streams ),
		need_to_initialize_( true )
  {
		if ( input_streams_.size() != 2 ) utility_exit_with_message( "PoseCombineSampleGenerator needs two sources of input!" );
		reset();
  }

	StepWisePoseCombineSampleGenerator::StepWisePoseCombineSampleGenerator(
				 utility::vector1< InputStreamWithResidueInfoOP > const & input_streams,
				 StepWisePoseSetupOP stepwise_pose_setup ):
		input_streams_( input_streams ),
		need_to_initialize_( true ),
		stepwise_pose_setup_( stepwise_pose_setup )
  {
		if ( input_streams_.size() != 2 ) utility_exit_with_message( "PoseCombineSampleGenerator needs two sources of input!" );
		reset();
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWisePoseCombineSampleGenerator::reset(){
		input_streams_[ 1 ]->reset();
		input_streams_[ 2 ]->reset();
		need_to_initialize_ = true;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWisePoseCombineSampleGenerator::has_another_sample(){
	 return ( input_streams_[ 1 ]->has_another_pose() ||
						input_streams_[ 2 ]->has_another_pose() );
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseCombineSampleGenerator::get_next_sample( core::pose::Pose & pose )
	{


		////////////////////////////////////////////////////////////
		if( input_streams_[ 2 ]->has_another_pose() ) {

			// This logic doesn't seem totally robust:
			if ( need_to_initialize_ ){
				input_streams_[ 1 ]->copy_next_pose_segment( pose );
				need_to_initialize_ = false;
			}

			input_streams_[ 2 ]->copy_next_pose_segment( pose );

		} else {

			if ( !input_streams_[ 1 ]->has_another_pose() ) utility_exit_with_message( "Asked PoseCombineSampleGenerator for another sample, but it did not have one!" );

			input_streams_[ 2 ]->reset();

			input_streams_[ 1 ]->copy_next_pose_segment( pose );
			input_streams_[ 2 ]->copy_next_pose_segment( pose );

		}

		// This is annoying but necessary:
		if ( stepwise_pose_setup_ && stepwise_pose_setup_->ready_to_align() ){
			//			std::cout << "ALIGNING POSE! " << std::endl;
			stepwise_pose_setup_->align_pose( pose );
			//			static Size count( 0 );
			//			pose.dump_pdb( "blah_"+ObjexxFCL::string_of( count++ )+".pdb");
		}

	}

	///////////////////////////////////////////////////////
	Size
	StepWisePoseCombineSampleGenerator::size() const{
		//we don't know a priori how many combinations there will be...
		// yes, I guess we could precompute this by looking through both streams ahead of time.
		return 0;
	}


	///////////////////////////////////////////////////////
	void
	StepWisePoseCombineSampleGenerator::set_stepwise_pose_setup( StepWisePoseSetupOP stepwise_pose_setup ){
		stepwise_pose_setup_ = stepwise_pose_setup;
	}

} //sample_generators
} //protein
} //sampling
} //stepwise
} //protocols
