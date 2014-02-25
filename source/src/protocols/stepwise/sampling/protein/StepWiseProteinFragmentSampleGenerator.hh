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


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinFragmentSampleGenerator_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinFragmentSampleGenerator_HH

#include <protocols/stepwise/sampling/protein/sample_generators/StepWisePoseSampleGenerator.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class StepWiseProteinFragmentSampleGenerator: public protocols::stepwise::sampling::protein::sample_generators::StepWisePoseSampleGenerator {
	public:

    StepWiseProteinFragmentSampleGenerator(
																					 std::string const frag_file,
																					 utility::vector1< core::Size > const & slice_res,
																					 utility::vector1< core::Size > const & moving_residues	 );

		void reset();

		bool has_another_sample();

		void get_next_sample( core::pose::Pose & pose );

		Size size() const;

	private:

		core::Size count_;
		core::Size insert_pos_;
		core::fragment::FrameOP frame;
		core::Size num_samples_;

  };

} //protein
} //sampling
} //stepwise
} //protocols

#endif

