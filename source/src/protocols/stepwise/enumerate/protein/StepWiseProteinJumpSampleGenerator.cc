// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinJumpSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseProteinJumpSampleGenerator.hh>
#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.hh>
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseProteinJumpSampleGenerator::StepWiseProteinJumpSampleGenerator(
																		 Size const which_jump,
																		 utility::vector1< core::kinematics::Jump > const & jumps )
  {
		initialize( which_jump, jumps );
		reset();
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinJumpSampleGenerator::initialize(
																		 Size const which_jump,
																		 utility::vector1< core::kinematics::Jump > const & jumps ){
		which_jump_ = which_jump;
		jumps_ = jumps;
		num_samples_ = jumps_.size();
		count_ = 0;
	}

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinJumpSampleGenerator::reset(){
		count_ = 0;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWiseProteinJumpSampleGenerator::has_another_sample(){
		return (count_ < num_samples_);
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinJumpSampleGenerator::get_next_sample( core::pose::Pose & pose )
	{

		count_++;

		if ( count_ > num_samples_ ) utility_exit_with_message( "Asked StepWiseProteinJumpSampleGenerator for another sample but it does not have one!" );

		if ( which_jump_ > 0 )	pose.set_jump( which_jump_, jumps_[ count_ ] );

	}

	///////////////////////////////////////////////////////
	Size
	StepWiseProteinJumpSampleGenerator::size() const{
		return num_samples_;
	}


} //protein
} //enumerate
} //stepwise
} //protocols
