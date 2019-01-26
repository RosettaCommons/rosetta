// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/jobs/MoverJob.cc
/// @brief A base class for any job that takes and uses a mover.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/jd3/jobs/MoverJob.hh>
#include <protocols/jd3/job_summaries/StandardPoseJobSummary.hh>
#include <protocols/jd3/job_results/PoseJobResult.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/simple_metrics/SimpleMetric.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>



static basic::Tracer TR( "protocols.jd3.jobs.MoverJob" );

namespace protocols {
namespace jd3 {
namespace jobs {

using namespace protocols::jd3;
using namespace protocols::jd3::job_summaries;
using namespace protocols::jd3::job_results;
using namespace core::simple_metrics;

//Constructor
MoverJob::MoverJob()
{}

//Destructor
MoverJob::~MoverJob()
{}

void
MoverJob::pose(core::pose::PoseOP pose){
	pose_ = pose;
}

core::pose::PoseOP
MoverJob::pose() {
	return pose_;
}

core::pose::PoseCOP
MoverJob::pose() const {
	return pose_;
}

moves::MoverOP
MoverJob::mover() {
	return mover_;
}

void
MoverJob::set_mover( moves::MoverOP mover){
	mover_ = mover;
}

void
MoverJob::add_mover( moves::MoverOP mover){

	if ( !mover_ || mover_->name() != "SequenceMover" ) {
		mover_ = moves::SequenceMoverOP( new moves::SequenceMover());
	}
	moves::SequenceMover & seq_mover = dynamic_cast<moves::SequenceMover & >( *mover_ );
	seq_mover.add_mover(mover);
}

void
MoverJob::add_metrics(utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics, std::string const & data_prefix  ) {

	if ( metrics_.count(data_prefix) ) {
		for ( auto const & metric : metrics ) {
			metrics_[data_prefix].push_back(metric);
		}
	} else {
		metrics_[data_prefix] = metrics;
	}
}

void
MoverJob::add_metric(core::simple_metrics::SimpleMetricCOP metric, std::string const & data_prefix){
	metrics_[data_prefix].push_back(metric);
}

jd3::CompletedJobOutput
MoverJob::run() {

	using namespace protocols::moves;

	assert(mover_);
	assert(pose_);

	mover()->apply( *pose_ );

	CompletedJobOutput job_output;

	// prepare the job status
	if ( mover()->get_last_move_status() == MS_SUCCESS ) {
		job_output.status = jd3_job_status_success;
	} else if ( mover()->get_last_move_status() == FAIL_RETRY ) {
		job_output.status = jd3_job_status_failed_retry;
	} else if ( mover()->get_last_move_status() == FAIL ) { // treat fail like fail retry?
		job_output.status = jd3_job_status_failed_retry;
	} else if ( mover()->get_last_move_status() == FAIL_DO_NOT_RETRY ) {
		job_output.status = jd3_job_status_failed_do_not_retry;
	} else if ( mover()->get_last_move_status() == FAIL_BAD_INPUT ) {
		job_output.status = jd3_job_status_inputs_were_bad;
	}

	// Retrieve all Poses created by the Mover and append them individually to
	// the CompletedJobOutput's job_results vector.
	while ( pose_ ) {

		//Run any set metrics in preparation for the JobSummary
		for ( auto const & metric_pair : metrics_ ) {
			run_metrics( *pose_, metric_pair.second, metric_pair.first);
		}

		//Create the JobResult and the JobSummary
		PoseJobResultOP job_result = create_job_result(pose_);
		StandardPoseJobSummaryOP summary = create_job_summary(*pose_);

		job_output.job_results.push_back( std::make_pair( summary, job_result ));

		// keep retrieving Poses from the Mover until the Mover stops returning any.
		// This could perhaps be expensive if the Mover wants to deliver 10K poses and
		// was not anticipating that they would all need to live in memory simultaneously!
		pose_ = mover()->get_additional_output();
	}

	return job_output;
}

PoseJobResultOP
MoverJob::create_job_result(core::pose::PoseOP pose){
	return PoseJobResultOP( new PoseJobResult(pose));
}

StandardPoseJobSummaryOP
MoverJob::create_job_summary( core::pose::Pose const & pose ){
	return StandardPoseJobSummaryOP( new StandardPoseJobSummary( pose ));
}

} //protocols
} //jd3
} //jobs
