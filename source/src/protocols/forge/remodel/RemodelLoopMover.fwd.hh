// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.fwd.hh
/// @brief  forward declaration for RemodelLoopMover
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelLoopMover_fwd_hh
#define INCLUDED_protocols_forge_remodel_RemodelLoopMover_fwd_hh


// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace remodel {


/// @brief forward declaration for protocols::forge::remodel::RemodelLoopMover
class RemodelLoopMover;


/// @brief access pointer for RemodelLoopMover
typedef utility::pointer::weak_ptr< RemodelLoopMover > RemodelLoopMoverAP;


/// @brief const access pointer for RemodelLoopMover
typedef utility::pointer::weak_ptr< RemodelLoopMover const > RemodelLoopMoverCAP;


/// @brief owning pointer for RemodelLoopMover
typedef utility::pointer::shared_ptr< RemodelLoopMover > RemodelLoopMoverOP;


/// @brief const owning pointer for RemodelLoopMover
typedef utility::pointer::shared_ptr< RemodelLoopMover const > RemodelLoopMoverCOP;


} // remodel
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_remodel_RemodelLoopMover_FWD_HH */
