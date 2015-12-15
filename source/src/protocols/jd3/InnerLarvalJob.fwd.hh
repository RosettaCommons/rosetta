// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/InnerLarvalJob.fwd.hh
/// @brief  Class declaration for the InnerLarvalJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_InnerLarvalJob_FWD_HH
#define INCLUDED_protocols_jd3_InnerLarvalJob_FWD_HH

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {

class InnerLarvalJob;

typedef utility::pointer::shared_ptr< InnerLarvalJob > InnerLarvalJobOP;
typedef utility::pointer::shared_ptr< InnerLarvalJob const > InnerLarvalJobCOP;

typedef utility::vector1< InnerLarvalJobOP > InnerLarvalJobs;

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_InnerLarvalJob_FWD_HH
