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


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinMainChainSampleGenerator_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinMainChainSampleGenerator_HH

#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/TorsionID.fwd.hh>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	class StepWiseProteinMainChainSampleGenerator: public protocols::stepwise::enumerate::protein::sample_generators::StepWisePoseSampleGenerator {
	public:

    StepWiseProteinMainChainSampleGenerator( 	 utility::vector1< core::id::TorsionID > const & which_torsions,
																							 utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists );

		void reset();

		bool has_another_sample();

		void get_next_sample( core::pose::Pose & pose );

		Size size() const;

	private:

    utility::vector1< core::id::TorsionID > which_torsions_;
		utility::vector1< utility::vector1< core::Real > > main_chain_torsion_set_lists_;
		Size count_;
		Size num_samples_;

  };

} //protein
} //enumerate
} //stepwise
} //protocols

#endif

