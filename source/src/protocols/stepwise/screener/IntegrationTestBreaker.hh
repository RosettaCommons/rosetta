// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/IntegrationTestBreaker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_IntegrationTestBreaker_HH
#define INCLUDED_protocols_stepwise_screener_IntegrationTestBreaker_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/AlignRMSD_Screener.fwd.hh>
#include <protocols/stepwise/screener/IntegrationTestBreaker.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class IntegrationTestBreaker: public StepWiseScreener {

public:

	//constructor
	IntegrationTestBreaker( StepWiseScreenerOP screener_whose_counts_to_check,
		StepWiseScreenerOP final_screener /*total_count -- for turning on align screen*/,
		AlignRMSD_ScreenerOP align_rmsd_screener );

	//destructor
	~IntegrationTestBreaker();

public:

	bool
	check_screen();

	std::string
	name() const { return "IntegrationTestBreaker"; }

	StepWiseScreenerType
	type() const { return INTEGRATION_TEST; }

	void
	fast_forward( sampler::StepWiseSamplerOP sampler );

private:

	StepWiseScreenerOP screener_whose_counts_to_check_;
	StepWiseScreenerOP final_screener_;
	AlignRMSD_ScreenerOP align_rmsd_screener_;

};

} //screener
} //stepwise
} //protocols

#endif
