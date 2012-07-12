// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/swa/SubToFullInfo.fwd.hh
/// @brief  Mapping from a working pose into a bigger pose, for swa monte carlo stuff.
/// @author Rhiju Das

#ifndef INCLUDED_protocols_swa_monte_carlo_SubToFullInfo_fwd_hh
#define INCLUDED_protocols_swa_monte_carlo_SubToFullInfo_fwd_hh


#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


// C++

namespace protocols {
namespace swa {
namespace monte_carlo {

	class SubToFullInfo;
	typedef utility::pointer::owning_ptr< SubToFullInfo > SubToFullInfoOP;

}
}
}

#endif
