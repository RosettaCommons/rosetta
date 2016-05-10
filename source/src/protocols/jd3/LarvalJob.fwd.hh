// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/LarvalJob.fwd.hh
/// @brief  The declaration for class protocols::jd3::LarvalJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_LarvalJob_FWD_HH
#define INCLUDED_protocols_jd3_LarvalJob_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <list>

namespace protocols {
namespace jd3 {

class LarvalJob;

typedef utility::pointer::shared_ptr< LarvalJob > LarvalJobOP;
typedef utility::pointer::shared_ptr< LarvalJob const > LarvalJobCOP;

typedef std::list< LarvalJobOP > LarvalJobs;

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_LarvalJob_FWD_HH
