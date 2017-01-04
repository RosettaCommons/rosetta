// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardInnerLarvalJob.fwd.hh
/// @brief  Class declaration for the StandardInnerLarvalJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_standard_StandardInnerLarvalJob_FWD_HH
#define INCLUDED_protocols_jd3_standard_StandardInnerLarvalJob_FWD_HH

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {
namespace standard {

class StandardInnerLarvalJob;

typedef utility::pointer::shared_ptr< StandardInnerLarvalJob > StandardInnerLarvalJobOP;
typedef utility::pointer::shared_ptr< StandardInnerLarvalJob const > StandardInnerLarvalJobCOP;

typedef utility::vector1< StandardInnerLarvalJobOP > StandardInnerLarvalJobs;

} // namesapce standard
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_standard_StandardInnerLarvalJob_FWD_HH
