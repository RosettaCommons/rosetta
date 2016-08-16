// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/components/RemodelMover.fwd.hh
/// @brief  forward declaration for RemodelMover
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelMover_fwd_hh
#define INCLUDED_protocols_forge_remodel_RemodelMover_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace remodel {


/// @brief forward declaration for RemodelMover
class RemodelMover;


/// @brief RemodelMover owning pointer
typedef utility::pointer::shared_ptr< RemodelMover > RemodelMover_OP;


/// @brief RemodelMover const owning pointer
typedef utility::pointer::shared_ptr< RemodelMover const > RemodelMover_COP;


/// @brief RemodelMover access pointer
typedef utility::pointer::weak_ptr< RemodelMover > RemodelMover_AP;


/// @brief RemodelMover const access pointer
typedef utility::pointer::weak_ptr< RemodelMover const > RemodelMover_CAP;


} // namespace remodel
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_remodel_RemodelMover_fwd_hh */
