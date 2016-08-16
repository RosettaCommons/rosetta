// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/Scorer.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_Scorer_HH
#define INCLUDED_protocols_stepwise_screener_Scorer_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/Scorer.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class Scorer: public StepWiseScreener {

public:

	//constructor
	Scorer();

	//constructor
	Scorer( core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scorefxn );

	//destructor
	~Scorer();

public:

	virtual
	bool
	check_screen();

	virtual
	std::string
	name() const { return "Scorer"; }

	virtual
	StepWiseScreenerType
	type() const { return SCORER; }

private:

	core::pose::Pose & pose_;
	core::scoring::ScoreFunctionCOP scorefxn_;

};

} //screener
} //stepwise
} //protocols

#endif
