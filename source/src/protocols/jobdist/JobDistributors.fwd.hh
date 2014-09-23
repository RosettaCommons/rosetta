// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/JobDistributors.hh
///
/// @brief
/// @author Ian W. Davis

#ifndef INCLUDED_protocols_jobdist_JobDistributors_fwd_hh
#define INCLUDED_protocols_jobdist_JobDistributors_fwd_hh

#include <utility/pointer/owning_ptr.hh>

#include <protocols/jobdist/Jobs.hh>

namespace protocols {
namespace jobdist {

class BaseJobDistributor;
class PlainSilentFileJobDistributor;

typedef utility::pointer::shared_ptr< BaseJobDistributor > BaseJobDistributorOP;
typedef utility::pointer::shared_ptr< PlainSilentFileJobDistributor > PlainSilentFileJobDistributorOP;

}
}

#endif
