// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/claims/VirtResClaim.fwd.hh
/// @brief Forward declaration for VirtResClaim in environment
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_claims_VirtResClaim_fwd_hh
#define INCLUDED_protocols_environment_claims_VirtResClaim_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>

namespace protocols {
namespace environment {
namespace claims {


class VirtResClaim;

// Types
typedef  utility::pointer::owning_ptr< VirtResClaim >  VirtResClaimOP;
typedef  utility::pointer::owning_ptr< VirtResClaim const >  VirtResClaimCOP;

typedef  utility::pointer::access_ptr< VirtResClaim >  VirtResClaimAP;
typedef  utility::pointer::access_ptr< VirtResClaim const >  VirtResClaimCAP;

typedef utility::vector1< VirtResClaimOP > VirtResClaims;

}
}
}

#endif
