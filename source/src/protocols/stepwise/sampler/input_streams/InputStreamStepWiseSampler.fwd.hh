// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/input_streams/InputStreamStepWiseSampler.fwd.hh
/// @brief 
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_input_streams_InputStreamStepWiseSampler_FWD_HH
#define INCLUDED_protocols_sampler_input_streams_InputStreamStepWiseSampler_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace input_streams {
	
	class InputStreamStepWiseSampler;
	typedef utility::pointer::owning_ptr< InputStreamStepWiseSampler > InputStreamStepWiseSamplerOP;
	typedef utility::pointer::owning_ptr< InputStreamStepWiseSampler const > InputStreamStepWiseSamplerCOP;
	
} //input_streams
} //sampler
} //stepwise
} //protocols

#endif
