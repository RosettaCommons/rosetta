// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/mmt_msd/MMTMinPackingJob.fwd.hh
/// @brief  forward declaration for class MMTMinPackingJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_mmt_msd_MMTMinPackingJob_FWD_HH
#define INCLUDED_devel_mmt_msd_MMTMinPackingJob_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace mmt_msd {

class MMTMinPackingJob;

typedef utility::pointer::owning_ptr< MMTMinPackingJob > MMTMinPackingJobOP;
typedef utility::pointer::owning_ptr< MMTMinPackingJob const > MMTMinPackingJobCOP;

}
}

#endif
