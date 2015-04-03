// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/bulge/BulgeApplyMover.fwd.hh
/// @brief 
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_bulge_BulgeApplyMover_FWD_HH
#define INCLUDED_protocols_stepwise_modeler_rna_bulge_BulgeApplyMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace bulge {
	
	class BulgeApplyMover;
	typedef utility::pointer::shared_ptr< BulgeApplyMover > BulgeApplyMoverOP;
	typedef utility::pointer::shared_ptr< BulgeApplyMover const > BulgeApplyMoverCOP;
	
} //bulge
} //rna
} //modeler
} //stepwise
} //protocols

#endif
