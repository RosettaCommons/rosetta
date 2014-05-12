// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_loop_close_StepWiseProteinCCD_Closer_FWD_HH
#define INCLUDED_protocols_stepwise_sampling_protein_loop_close_StepWiseProteinCCD_Closer_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace loop_close {

	class StepWiseProteinCCD_Closer;
	typedef utility::pointer::owning_ptr< StepWiseProteinCCD_Closer > StepWiseProteinCCD_CloserOP;
	typedef utility::pointer::owning_ptr< StepWiseProteinCCD_Closer const > StepWiseProteinCCD_CloserCOP;

} //loop_close
} //protein
} //sampling
} //stepwise
} //protocols

#endif
