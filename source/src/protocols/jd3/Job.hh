// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/Job.hh
/// @brief  The definition of class protocols::jd3::Job
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_Job_HH
#define INCLUDED_protocols_jd3_Job_HH

// Unit headers
#include <protocols/jd3/Job.fwd.hh>

// Package headers
#include <protocols/jd3/JobResult.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//C++ headers
#include <string>
#include <list>
#include <map>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace jd3 {

/// @brief %protocols::jd3::Job, which is created by the JobQueen from a LarvalJob, is
/// responsible for executing some protocol from beginning to end in its "run" method.
/// It must write all of its output data to the JobResult object that it creates and
/// returns.
///
/// @details Note that in JD2, the currently running Job was globally accessible.  That is not
/// the case in JD3.  Any data that should be output by a running job must be tucked
/// inside the JobResult to be eventually retrieved by the JobQueen or by a JobDataOutputter
/// at the end of execution.
class Job : public utility::pointer::ReferenceCount
{
public:

	Job();
	virtual ~Job();

	/// @brief This is the main function of the Job object. It will be invoked by the JobDistributor.
	/// The Job will return a JobResult at the conclusion of its execution, and the JobResult will
	/// be serialized and sent to the appropriate JobQueen for processing and output.  The Job itself
	/// will be discarded.  Large constant data that might be shared between multiple jobs can be
	/// held by the Job object (or by classes that the Job object holds) but that data should not be
	/// put into the JobResult -- large data should not be repeatedly serialized and shipped between
	/// nodes.
	virtual
	JobResultOP run() = 0;

}; // Job

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_Job_HH
