// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseCombineSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/sample_generators/StepWiseCombineSampleGenerator.hh>
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
	StepWiseCombineSampleGenerator::StepWiseCombineSampleGenerator(
																																				 sample_generators::StepWisePoseSampleGeneratorOP first_generator,
																																				 sample_generators::StepWisePoseSampleGeneratorOP second_generator ):
		first_generator_( first_generator ),
		second_generator_( second_generator ),
		need_to_initialize_( true )
  {
		reset();
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseCombineSampleGenerator::reset(){
		first_generator_->reset();
		second_generator_->reset();
		need_to_initialize_ = true;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWiseCombineSampleGenerator::has_another_sample(){
		return ( first_generator_->has_another_sample() ||
						 second_generator_->has_another_sample() );
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseCombineSampleGenerator::get_next_sample( core::pose::Pose & pose )
	{


		////////////////////////////////////////////////////////////
		if( second_generator_->has_another_sample() ) {

			// This logic doesn't seem totally robust:
			if ( need_to_initialize_ ){
				first_generator_->get_next_sample( pose );
				need_to_initialize_ = false;
			}

			second_generator_->get_next_sample( pose );

		} else {

			if ( !first_generator_->has_another_sample() ) utility_exit_with_message( "Asked CombineSampleGenerator for another sample, but it did not have one!" );

			second_generator_->reset();

			first_generator_->get_next_sample( pose );
			second_generator_->get_next_sample( pose );

		}

	}

	///////////////////////////////////////////////////////
	Size
	StepWiseCombineSampleGenerator::size() const{
		//we don't know a priori how many combinations there will be...
		// yes, I guess we could precompute this by looking through both streams ahead of time.
		return 0;
	}


} //sample_generators
} //protein
} //sampling
} //stepwise
} //protocols
