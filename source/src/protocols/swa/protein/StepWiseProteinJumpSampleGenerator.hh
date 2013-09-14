// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_protein_StepWiseProteinJumpSampleGenerator_HH
#define INCLUDED_protocols_swa_protein_StepWiseProteinJumpSampleGenerator_HH

#include <protocols/swa/sample_generators/StepWisePoseSampleGenerator.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/kinematics/Jump.hh>

namespace protocols {
namespace swa {
namespace protein {

	class StepWiseProteinJumpSampleGenerator: public protocols::swa::StepWisePoseSampleGenerator {
	public:

    StepWiseProteinJumpSampleGenerator(
													 Size const which_jump,
													 utility::vector1< core::kinematics::Jump > const & jumps );

		void reset();

		bool has_another_sample();

		void get_next_sample( core::pose::Pose & pose );

		Size size() const;

		void initialize( core::Size const which_jump,
										 utility::vector1< core::kinematics::Jump > const & jumps );


	private:

		Size which_jump_;
		utility::vector1< core::kinematics::Jump > jumps_;
		Size count_;
		Size num_samples_;

  };

} //protein
} //swa
} // protocols

#endif

