// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinFragmentSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/protein/StepWiseProteinFragmentSampleGenerator.hh>
#include <protocols/swa/protein/StepWiseProteinUtil.hh>
#include <protocols/swa/sample_generators/StepWisePoseSampleGenerator.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <utility/exit.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseProteinFragmentSampleGenerator::StepWiseProteinFragmentSampleGenerator(
																					 std::string const frag_file,
																					 utility::vector1< core::Size > const & slice_res,
																					 utility::vector1< core::Size > const & moving_residues	 ):
		count_( 0 )
  {

		core::fragment::ConstantLengthFragSetOP fragset =  new core::fragment::ConstantLengthFragSet( 0 /*frag_length ... is reset by reader*/, frag_file );

		if( fragset->max_frag_length() != moving_residues.size() ) {
			utility_exit_with_message( "Number of -moving_res must match frag size!" );
		}
		fragment_set_slice( fragset, slice_res );


		core::fragment::FrameList frames;
		insert_pos_ = moving_residues[ 1 ];
		fragset->frames( insert_pos_, frames );

		frame = frames[1];
		num_samples_ = frame->nr_frags();

  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinFragmentSampleGenerator::reset(){
		count_ = 0;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWiseProteinFragmentSampleGenerator::has_another_sample(){
		return (count_ < num_samples_);
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinFragmentSampleGenerator::get_next_sample( core::pose::Pose & pose )
	{
		count_++;
		frame->apply( count_, pose );
	}

	///////////////////////////////////////////////////////
	Size
	StepWiseProteinFragmentSampleGenerator::size() const{
		return num_samples_;
	}


}
}
}
