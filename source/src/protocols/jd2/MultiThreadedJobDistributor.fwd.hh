// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MultiThreadedJobDistributor.fwd.hh
/// @brief  Job distributor that luanches threads to carry out independent jobs
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_protocols_jd2_MultiThreadedJobDistributor_fwd_hh
#define INCLUDED_protocols_jd2_MultiThreadedJobDistributor_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class RunningJob;
class MTJob;
class MTJobGroup;

class MultiThreadedJobDistributor;

typedef utility::pointer::owning_ptr< MTJob > MTJobOP;
typedef utility::pointer::owning_ptr< MTJob const > MTJobCOP;

typedef utility::pointer::owning_ptr< MTJobGroup > MTJobGroupOP;
typedef utility::pointer::owning_ptr< MTJobGroup const > MTJobGroupCOP;

typedef utility::pointer::owning_ptr< RunningJob > RunningJobOP;
typedef utility::pointer::owning_ptr< RunningJob const > RunningJobCOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_MultiThreadedJobDistributor_FWD_HH
