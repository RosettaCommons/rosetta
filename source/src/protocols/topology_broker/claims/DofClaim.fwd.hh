// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_claims_DofClaim_fwd_hh
#define INCLUDED_protocols_topology_broker_claims_DofClaim_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace topology_broker {
namespace claims {

// Forward
class DofClaim;

// Types
typedef  utility::pointer::shared_ptr< DofClaim >  DofClaimOP;
typedef  utility::pointer::shared_ptr< DofClaim const >  DofClaimCOP;

typedef  utility::pointer::weak_ptr< DofClaim >  DofClaimAP;
typedef  utility::pointer::weak_ptr< DofClaim const >  DofClaimCAP;

typedef utility::vector1< DofClaimOP > DofClaims;


} // namespace claims
} // namespace kinematics
} // namespace core

#endif
