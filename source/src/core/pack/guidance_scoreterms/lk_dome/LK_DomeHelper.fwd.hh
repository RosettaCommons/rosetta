// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/pack/guidance_scoreterms/lk_dome/LK_DomeHelper.fwd.hh
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_scoring_lkball_LK_DomeHelper_fwd_hh
#define INCLUDED_core_scoring_lkball_LK_DomeHelper_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace lk_dome {


class LK_DomeHelper;

typedef utility::pointer::shared_ptr< LK_DomeHelper > LK_DomeHelperOP;
typedef utility::pointer::shared_ptr< LK_DomeHelper const > LK_DomeHelperCOP;


}
} //guidance_scoreterms
} //pack
} //core

#endif
