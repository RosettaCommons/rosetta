// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/ProteinCCD_ClosureScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ProteinCCD_ClosureScreener_HH
#define INCLUDED_protocols_stepwise_screener_ProteinCCD_ClosureScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/ProteinCCD_ClosureScreener.fwd.hh>
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class ProteinCCD_ClosureScreener: public SampleApplier {

public:

	//constructor
	ProteinCCD_ClosureScreener( modeler::protein::loop_close::StepWiseProteinCCD_CloserOP ccd_closer,
		pose::Pose & screening_pose );

	//destructor
	~ProteinCCD_ClosureScreener();

public:

	std::string
	name() const { return "ProteinCCD_ClosureScreener"; }

	StepWiseScreenerType
	type() const { return PROTEIN_CCD_CLOSURE; }

	bool
	check_screen();

	void
	add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

private:

	modeler::protein::loop_close::StepWiseProteinCCD_CloserOP ccd_closer_;

};

} //screener
} //stepwise
} //protocols

#endif
