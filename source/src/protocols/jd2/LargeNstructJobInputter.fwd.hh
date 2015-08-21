// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/LargeNstructJobInputter.fwd.hh
/// @brief  Owning pointers and so forth for a JobInputter for cases where the number of jobs is large enough to fill up memory.
/// @author Vikram K. Mulligan, Baker Laboratory (vmullig@uw.edu)

#ifndef INCLUDED_protocols_jd2_LargeNstructJobInputter_fwd_hh
#define INCLUDED_protocols_jd2_LargeNstructJobInputter_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class LargeNstructJobInputter;
typedef utility::pointer::weak_ptr< LargeNstructJobInputter > LargeNstructJobInputterAP;
typedef utility::pointer::weak_ptr< LargeNstructJobInputter const > LargeNstructJobInputterCAP;
typedef utility::pointer::shared_ptr< LargeNstructJobInputter > LargeNstructJobInputterOP;
typedef utility::pointer::shared_ptr< LargeNstructJobInputter const > LargeNstructJobInputterCOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_LargeNstructJobInputter_FWD_HH
