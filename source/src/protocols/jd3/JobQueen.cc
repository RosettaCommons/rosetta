// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/JobQueen.cc
/// @brief  The (mostly empty) definition of the JobQueen class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//unit headers
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/JobTracker.hh>
#include <protocols/jd3/JobDigraph.hh>



namespace protocols {
namespace jd3 {

JobQueen::JobQueen(){
	job_tracker_ = utility::pointer::make_shared< JobTracker >();
}

JobQueen::~JobQueen() = default;

JobDigraphOP
JobQueen::create_and_set_initial_job_dag(){
	JobDigraphOP job_graph = create_initial_job_dag();
	job_graph_ = job_graph;

	return job_graph;
}

void
JobQueen::note_job_completed_and_track(protocols::jd3::LarvalJobCOP job, protocols::jd3::JobStatus status, core::Size nresults){
	JQKey key;
	job_tracker_->track_completed_job(key, *job, status);
	note_job_completed(job, status, nresults);
}

LarvalJobs
JobQueen::determine_job_list_and_track(core::Size job_dag_node_index, core::Size max_n_jobs){
	JQKey key;
	LarvalJobs larval_jobs = determine_job_list(job_dag_node_index, max_n_jobs);
	job_tracker_->track_starting_job_list(key, larval_jobs);
	return larval_jobs;
}

JobDigraph const &
JobQueen::get_job_graph() const{
	return *job_graph_;
}


JobTracker &
JobQueen::get_job_tracker() {
	return *job_tracker_;
}

JobTracker const &
JobQueen::get_job_tracker() const {
	return *job_tracker_;
}

} // namespace jd3
} // namespace protocols



