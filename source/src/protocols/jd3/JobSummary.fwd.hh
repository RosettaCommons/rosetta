// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobSummary.fwd.hh
/// @brief  The forward declaration for class protocols::jd3::JobSummary
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobSummary_FWD_HH
#define INCLUDED_protocols_jd3_JobSummary_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {

class JobSummary;

typedef utility::pointer::shared_ptr< JobSummary > JobSummaryOP;
typedef utility::pointer::shared_ptr< JobSummary const > JobSummaryCOP;

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_JobSummary_FWD_HH
