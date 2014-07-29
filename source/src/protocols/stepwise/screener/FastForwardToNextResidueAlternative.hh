// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/FastForwardToNextResidueAlternative.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_FastForwardToNextResidueAlternative_HH
#define INCLUDED_protocols_stepwise_screener_FastForwardToNextResidueAlternative_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/FastForwardToNextResidueAlternative.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

	class FastForwardToNextResidueAlternative: public StepWiseScreener {

	public:

		//constructor
		FastForwardToNextResidueAlternative( Size const moving_res );

		//destructor
		~FastForwardToNextResidueAlternative();

	public:

		//		bool
		//		check_screen();

		virtual
		std::string
		name() const { return "FastForwardToNextResidueAlternative"; }

		virtual
		StepWiseScreenerType
		type() const { return FAST_FORWARD_TO_NEXT_RESIDUE_ALTERNATIVE; }

		virtual
		void
		get_update( sampler::StepWiseSamplerBaseOP sampler );

		// kind of tricky -- put fast_forward above in get_update.
		//		virtual
		//		void
		//		fast_forward( sampler::StepWiseSamplerBaseOP sampler );

	private:

		Size const moving_res_;

	};

} //screener
} //stepwise
} //protocols

#endif
