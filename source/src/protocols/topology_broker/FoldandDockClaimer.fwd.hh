// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   FoldandDockClaimer.fwd.hh
/// @brief  FoldandDockClaimer forward declarations header
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_topology_broker_FoldandDockClaimer_fwd_hh
#define INCLUDED_protocols_topology_broker_FoldandDockClaimer_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace topology_broker {

// Forward
class FoldandDockClaimer;

// Types
typedef  utility::pointer::shared_ptr< FoldandDockClaimer >  FoldandDockClaimerOP;
typedef  utility::pointer::shared_ptr< FoldandDockClaimer const >  FoldandDockClaimerCOP;

typedef  utility::pointer::weak_ptr< FoldandDockClaimer >  FoldandDockClaimerAP;
typedef  utility::pointer::weak_ptr< FoldandDockClaimer const >  FoldandDockClaimerCAP;


} // namespace kinematics
} // namespace core

#endif
