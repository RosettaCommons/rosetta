// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/claims/JumpClaim.fwd.hh
/// @brief Forward declaration for JumpClaim in environment
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_claims_JumpClaim_fwd_hh
#define INCLUDED_protocols_environment_claims_JumpClaim_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>

namespace protocols {
namespace environment {
namespace claims {


class JumpClaim;

// Types
typedef  utility::pointer::owning_ptr< JumpClaim >  JumpClaimOP;
typedef  utility::pointer::owning_ptr< JumpClaim const >  JumpClaimCOP;

typedef  utility::pointer::access_ptr< JumpClaim >  JumpClaimAP;
typedef  utility::pointer::access_ptr< JumpClaim const >  JumpClaimCAP;

typedef utility::vector1< JumpClaimOP > JumpClaims;

}
}
}

#endif
