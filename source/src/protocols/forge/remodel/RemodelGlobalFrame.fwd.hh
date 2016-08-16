// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/components/RemodelGlobalFrame.fwd.hh
/// @brief  forward declaration for RemodelGlobalFrame
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelGlobalFrame_fwd_hh
#define INCLUDED_protocols_forge_remodel_RemodelGlobalFrame_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace remodel {


/// @brief forward declaration for RemodelGlobalFrame
class RemodelGlobalFrame;


/// @brief RemodelGlobalFrame owning pointer
typedef utility::pointer::shared_ptr< RemodelGlobalFrame > RemodelGlobalFrame_OP;


/// @brief RemodelGlobalFrame const owning pointer
typedef utility::pointer::shared_ptr< RemodelGlobalFrame const > RemodelGlobalFrame_COP;


/// @brief RemodelGlobalFrame access pointer
typedef utility::pointer::shared_ptr< RemodelGlobalFrame > RemodelGlobalFrame_AP;


/// @brief RemodelGlobalFrame const access pointer
typedef utility::pointer::shared_ptr< RemodelGlobalFrame const > RemodelGlobalFrame_CAP;


} // namespace remodel
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_RemodelGlobalFrame_FWD_HH */
