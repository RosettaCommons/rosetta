// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/ThreadingJob.fwd.hh
/// @brief  header file for ThreadingJob classes, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_ThreadingJob_fwd_hh
#define INCLUDED_protocols_jd2_ThreadingJob_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

namespace protocols {
namespace jd2 {

class ThreadingJob;
typedef utility::pointer::owning_ptr< ThreadingJob > ThreadingJobOP;
typedef utility::pointer::owning_ptr< ThreadingJob const > ThreadingJobCOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_ThreadingJob_FWD_HH
