// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 ///
 /// @file protocols/topology_broker/PcsEnergyController.fwd.hh
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

#ifndef INCLUDED_protocols_topology_broker_PcsEnergyController_fwd_hh
#define INCLUDED_protocols_topology_broker_PcsEnergyController_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace topology_broker {

class PcsEnergyController;

typedef  utility::pointer::shared_ptr< PcsEnergyController >  PcsEnergyControllerOP;
typedef  utility::pointer::shared_ptr< PcsEnergyController const >  PcsEnergyControllerCOP;

typedef  utility::pointer::weak_ptr< PcsEnergyController >  PcsEnergyControllerAP;
typedef  utility::pointer::weak_ptr< PcsEnergyController const >  PcsEnergyControllerCAP;

} // namespace topology_broker
} // namespace protocols

#endif
