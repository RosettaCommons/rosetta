// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file PseudocontactShiftEnergyController.fwd.hh
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz & Oliver Lange
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_topology_broker_PseudocontactShiftEnergyController_fwd_hh
#define INCLUDED_protocols_topology_broker_PseudocontactShiftEnergyController_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace topology_broker {

// Forward
class PseudocontactShiftEnergyController;

// Types
typedef  utility::pointer::shared_ptr< PseudocontactShiftEnergyController >  PseudocontactShiftEnergyControllerOP;
typedef  utility::pointer::shared_ptr< PseudocontactShiftEnergyController const >  PseudocontactShiftEnergyControllerCOP;

typedef  utility::pointer::weak_ptr< PseudocontactShiftEnergyController >  PseudocontactShiftEnergyControllerAP;
typedef  utility::pointer::weak_ptr< PseudocontactShiftEnergyController const >  PseudocontactShiftEnergyControllerCAP;

} // namespace topology_broker
} // namespace protocols

#endif
