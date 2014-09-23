// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/ImplicitFastClashCheck.fwd.hh
/// @brief  does implicit fast clash checking WRT the provided pose
/// @author Will Sheffler (will@sheffler.me)

#ifndef INCLUDED_protocols_scoring_ImplicitFastClashCheck_fwd_hh
#define INCLUDED_protocols_scoring_ImplicitFastClashCheck_fwd_hh


#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace scoring {

class ImplicitFastClashCheck;

typedef utility::pointer::shared_ptr< ImplicitFastClashCheck > ImplicitFastClashCheckOP;
typedef utility::pointer::shared_ptr< ImplicitFastClashCheck const > ImplicitFastClashCheckCOP;

typedef utility::pointer::weak_ptr< ImplicitFastClashCheck > ImplicitFastClashCheckAP;
typedef utility::pointer::weak_ptr< ImplicitFastClashCheck const > ImplicitFastClashCheckCAP;

} // namespace scoring
} // namespace protocols

#endif // INCLUDED_protocols_scoring_methods_ImplicitFastClashCheck_fwd_hh
