// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseDoNothingSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/sample_generators/StepWiseDoNothingSampleGenerator.hh>
#include <core/types.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

// Fires once. Doesn't actually do anything to the pose.
namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {
namespace sample_generators {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseDoNothingSampleGenerator::StepWiseDoNothingSampleGenerator():
		ready_to_go_( 1 )
  {
		reset();
  }


  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseDoNothingSampleGenerator::reset(){
		ready_to_go_ = 1;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWiseDoNothingSampleGenerator::has_another_sample(){
		return ready_to_go_;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseDoNothingSampleGenerator::get_next_sample( core::pose::Pose &  )
	{
		ready_to_go_ = 0;
	}

	///////////////////////////////////////////////////////
	Size
	StepWiseDoNothingSampleGenerator::size() const{
		return 1;
	}


} //sample_generators
} //protein
} //enumerate
} //stepwise
} //protocols
