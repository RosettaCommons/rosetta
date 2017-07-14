// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/job_distributors/VanillaJobDistributor.cc
/// @brief  VanillaJobDistributor class definition: a single process on a single processor with one thread
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_job_distributors_VanillaJobDistributor_HH
#define INCLUDED_protocols_jd3_job_distributors_VanillaJobDistributor_HH

// Unit headers
#include <protocols/jd3/job_distributors/VanillaJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd3/CompletedJobOutput.fwd.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobQueen.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>


// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.fwd.hh>

// C++ headers
#include <map>

namespace protocols {
namespace jd3 {
namespace job_distributors {

/// @brief The VanillaJobDistributor is a single process running by itself, running a single thread.
class VanillaJobDistributor : public JobDistributor {
public:
	typedef std::list< core::Size > SizeList;
	typedef std::list< JobResultID > JobResultIDList;
	typedef std::map< JobResultID, std::pair< LarvalJobOP, JobResultOP > > JobResultMap;

public:

	VanillaJobDistributor();
	virtual ~VanillaJobDistributor();

	/// @brief The main method for executing a protocol.
	virtual
	void
	go( JobQueenOP queen );

private:

	// subroutines part of go()
	void run_jobs_for_dag_node( core::Size job_node );
	utility::vector1< JobResultCOP > construct_job_result_input_list( LarvalJobCOP larval_job );
	CompletedJobOutput run_mature_job( LarvalJobOP larval_job, JobOP mature_job );
	void potentially_output_some_job_results();
	void potentially_discard_some_job_results();


private:

	JobQueenOP job_queen_;
	JobDigraphOP job_dag_;

	// The big old map that stores all the JobResults that are generated
	// over the course of execution
	JobResultMap job_results_;

};

} // job_distributors
} // jd3
} // protocols

#endif //INCLUDED_protocols_jd3_JobDistributor_HH
