// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/components/RemodelLigandHandler.fwd.hh
/// @brief  forward declaration for RemodelLigandHandler
/// @author possu huang (possu@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelLigandHandler_fwd_hh
#define INCLUDED_protocols_forge_remodel_RemodelLigandHandler_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace forge {
namespace remodel {


/// @brief forward declaration for RemodelLigandHandler
class RemodelLigandHandler;


/// @brief RemodelLigandHandler owning pointer
typedef utility::pointer::shared_ptr< RemodelLigandHandler > RemodelLigandHandler_OP;


/// @brief RemodelLigandHandler const owning pointer
typedef utility::pointer::shared_ptr< RemodelLigandHandler const > RemodelLigandHandler_COP;


/// @brief RemodelLigandHandler access pointer
typedef utility::pointer::shared_ptr< RemodelLigandHandler > RemodelLigandHandler_AP;


/// @brief RemodelLigandHandler const access pointer
typedef utility::pointer::shared_ptr< RemodelLigandHandler const > RemodelLigandHandler_CAP;


} // namespace remodel
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_RemodelLigandHandler_FWD_HH */
