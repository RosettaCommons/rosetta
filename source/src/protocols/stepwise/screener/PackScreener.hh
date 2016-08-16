// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/PackScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_PackScreener_HH
#define INCLUDED_protocols_stepwise_screener_PackScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/PackScreener.fwd.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

class PackScreener:  public SampleApplier {

public:

	//constructor
	PackScreener( pose::Pose & pose,
		modeler::packer::StepWisePackerOP stepwise_packer );

	//destructor
	~PackScreener();

public:

	bool
	check_screen();

	std::string
	name() const { return "PackScreener"; }

	StepWiseScreenerType
	type() const { return PACK; }

	void
	add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

private:

	modeler::packer::StepWisePackerOP stepwise_packer_;

};

} //screener
} //stepwise
} //protocols

#endif
