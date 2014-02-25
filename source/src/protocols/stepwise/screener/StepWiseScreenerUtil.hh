// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StepWiseScreenerUtil.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_StepWiseScreenerUtil_HH
#define INCLUDED_protocols_stepwise_screener_StepWiseScreenerUtil_HH

#include <protocols/stepwise/screener/StepWiseScreenerUtil.fwd.hh>
#include <protocols/rotamer_sampler/RotamerBase.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {


	void
	fast_forward_to_next_residue_pair( rotamer_sampler::RotamerBaseOP sampler,
																		 Size const res1,
																		 Size const res2 );

} //screener
} //stepwise
} //protocols

#endif
