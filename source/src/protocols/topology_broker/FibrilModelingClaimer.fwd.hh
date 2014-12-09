// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   FibrilModelingClaimer.fwd.hh
/// @brief  FibrilModelingClaimer forward declarations header
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_topology_broker_FibrilModelingClaimer_fwd_hh
#define INCLUDED_protocols_topology_broker_FibrilModelingClaimer_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace topology_broker {

// Forward
class FibrilModelingClaimer;

// Types
typedef  utility::pointer::shared_ptr< FibrilModelingClaimer >  FibrilModelingClaimerOP;
typedef  utility::pointer::shared_ptr< FibrilModelingClaimer const >  FibrilModelingClaimerCOP;

typedef  utility::pointer::weak_ptr< FibrilModelingClaimer >  FibrilModelingClaimerAP;
typedef  utility::pointer::weak_ptr< FibrilModelingClaimer const >  FibrilModelingClaimerCAP;


} // namespace kinematics
} // namespace core

#endif
