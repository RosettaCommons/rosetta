// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/legacy/screener/ProteinAtrRepScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ProteinAtrRepScreener_HH
#define INCLUDED_protocols_stepwise_screener_ProteinAtrRepScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/legacy/screener/ProteinAtrRepScreener.fwd.hh>
#include <protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

class ProteinAtrRepScreener: public stepwise::screener::SampleApplier {

public:
	//constructor
	ProteinAtrRepScreener( core::pose::Pose & pose_atr_rep_screen,
		protocols::stepwise::modeler::protein::checker::ProteinAtrRepCheckerOP atr_rep_checker );

	//destructor
	~ProteinAtrRepScreener();

public:

	std::string
	name() const { return "ProteinAtrRepScreener"; }

	stepwise::screener::StepWiseScreenerType
	type() const { return stepwise::screener::PROTEIN_ATR_REP; }

	bool
	check_screen();

private:

	protocols::stepwise::modeler::protein::checker::ProteinAtrRepCheckerOP atr_rep_checker_;

};

} //screener
} //legacy
} //stepwise
} //protocols

#endif
