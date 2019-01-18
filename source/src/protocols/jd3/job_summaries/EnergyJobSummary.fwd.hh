// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/job_summaries/EnergyJobSummary.fwd.hh
/// @brief A JobSummary that simply stores the energy from a job.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_jd3_job_summaries_EnergyJobSummary_fwd_hh
#define INCLUDED_protocols_jd3_job_summaries_EnergyJobSummary_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace jd3 {
namespace job_summaries {

class EnergyJobSummary;

typedef utility::pointer::shared_ptr< EnergyJobSummary > EnergyJobSummaryOP;
typedef utility::pointer::shared_ptr< EnergyJobSummary const > EnergyJobSummaryCOP;

} //protocols
} //jd3
} //job_summaries

#endif //INCLUDED_protocols_jd3_job_summaries_EnergyJobSummary_fwd_hh
