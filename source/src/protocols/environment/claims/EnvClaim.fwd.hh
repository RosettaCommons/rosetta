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
/// @author Justin Porter


#ifndef INCLUDED_protocols_environment_claims_EnvClaim_fwd_hh
#define INCLUDED_protocols_environment_claims_EnvClaim_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <list>

namespace protocols {
namespace environment {
namespace claims {

// Forward
class EnvClaim;

// Types
typedef  utility::pointer::shared_ptr< EnvClaim >  EnvClaimOP;
typedef  utility::pointer::shared_ptr< EnvClaim const >  EnvClaimCOP;

typedef  utility::pointer::weak_ptr< EnvClaim >  EnvClaimAP;
typedef  utility::pointer::weak_ptr< EnvClaim const >  EnvClaimCAP;

typedef std::list< EnvClaimOP > EnvClaims;

} // namespace claims
} // namespace environment
} // namespace protocols

#endif
