// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/MRSJobSummary.fwd.hh
/// @brief  The declaration for class protocols::multistage_rosetta_scripts::MRSJobSummary
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_MRSJobSummary_FWD_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_MRSJobSummary_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace multistage_rosetta_scripts {

class MRSJobSummary;

typedef utility::pointer::shared_ptr< MRSJobSummary > MRSJobSummaryOP;
typedef utility::pointer::shared_ptr< MRSJobSummary const > MRSJobSummaryCOP;

} // namespace multistage_rosetta_scripts
} // namespace protocols

#endif //INCLUDED
