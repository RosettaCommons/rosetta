// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerBase.fwd.hh
/// @brief Abstract Base Class for StepWiseSampler generator.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_StepWiseSamplerBase_fwd_HH
#define INCLUDED_protocols_sampler_StepWiseSamplerBase_fwd_HH

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerBase;
typedef utility::pointer::shared_ptr< StepWiseSamplerBase > StepWiseSamplerBaseOP;

} //sampler
} //stepwise
} //protocols

#endif
