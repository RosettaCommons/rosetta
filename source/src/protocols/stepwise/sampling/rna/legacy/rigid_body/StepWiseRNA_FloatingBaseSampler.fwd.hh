// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/legacy/rigid_body/StepWiseRNA_FloatingBaseSampler.fwd.hh
/// @brief 
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSampler_FWD_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSampler_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {
namespace rigid_body {
	
	class StepWiseRNA_FloatingBaseSampler;
	typedef utility::pointer::owning_ptr< StepWiseRNA_FloatingBaseSampler > StepWiseRNA_FloatingBaseSamplerOP;
	typedef utility::pointer::owning_ptr< StepWiseRNA_FloatingBaseSampler const > StepWiseRNA_FloatingBaseSamplerCOP;
	
} //rigid_body
} //legacy
} //rna
} //sampling
} //stepwise
} //protocols

#endif
