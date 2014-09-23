// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/Jobs.fwd.hh
///
/// @brief  Objects representing inputs for the JobDistributor
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_jobdist_Jobs_fwd_hh
#define INCLUDED_protocols_jobdist_Jobs_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jobdist {

class BasicJob; // fwd declaration
typedef utility::pointer::shared_ptr< BasicJob > BasicJobOP;
typedef utility::pointer::shared_ptr< BasicJob const > BasicJobCOP;

} // namespace jobdist
} // namespace protocols

#endif // INCLUDED_protocols_jobdist_Jobs_FWD_HH
