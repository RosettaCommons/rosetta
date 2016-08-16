// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MultiThreadingJob.cc
/// @author James Thompson

#include <protocols/jd2/MultiThreadingJob.hh>
#include <protocols/jd2/MultiThreadingJob.fwd.hh>


#include <core/types.hh>

#include <protocols/jd2/InnerMultiThreadingJob.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

InnerMultiThreadingJobOP MultiThreadingJob::multi_threading_inner_job() {
	return inner_job_;
}

MultiThreadingJob::~MultiThreadingJob() {}

} // namespace jd2
} // namespace protocols
