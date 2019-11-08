// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/job_distributors/VanillaJobDistributor.cc
/// @brief  VanillaJobDistributor class method definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <protocols/jd3/job_distributors/VanillaJobDistributor.hh>

// Package headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/JobSummary.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/job_distributors/JobExtractor.hh>
#include <protocols/jd3/output/OutputSpecification.hh>
#include <protocols/jd3/output/ResultOutputter.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/crash_report.hh>

// C++ headers
#include <string>
#include <ctime>
#include <fstream>

namespace protocols {
namespace jd3 {
namespace job_distributors {

static basic::Tracer TR( "protocols.jd3.VanillaJobDistributor" );

VanillaJobDistributor::VanillaJobDistributor() = default;

VanillaJobDistributor::~VanillaJobDistributor() = default;

void
VanillaJobDistributor::go( JobQueenOP queen ) {
	using core::Size;
	using namespace utility::graph;

	// Warn about the Vanilla JQ not performing any checkpointing.
	if ( basic::options::option[ basic::options::OptionKeys::jd3::checkpoint ] ||
			basic::options::option[ basic::options::OptionKeys::jd3::checkpoint_period ].user() ||
			basic::options::option[ basic::options::OptionKeys::jd3::checkpoint_directory ].user() ||
			basic::options::option[ basic::options::OptionKeys::jd3::restore_from_checkpoint ] ||
			basic::options::option[ basic::options::OptionKeys::jd3::keep_all_checkpoints ].user() ||
			basic::options::option[ basic::options::OptionKeys::jd3::keep_checkpoint ].user() ) {
		TR.Warning << "Checkpointing options detected on the command line, but VanillaJobDistributor does not perform any checkpointing" << std::endl;
	}

	if ( basic::options::option[ basic::options::OptionKeys::jd3::job_definition_schema ].user() ) {
		std::string xsd = queen->job_definition_xsd();
		std::string outfile_name = basic::options::option[ basic::options::OptionKeys::jd3::job_definition_schema ]();
		std::ofstream ofs( outfile_name.c_str() );
		ofs << xsd;
		return;
	}

	job_queen_ = queen;
	auto job_extractor = std::make_shared< JobExtractor >();
	job_extractor->set_job_queen( queen );
	job_dag_ = job_extractor->get_initial_job_dag_and_queue();

	while ( job_extractor->not_done() ) {
		LarvalJobOP larval_job = job_extractor->pop_job_from_queue();

		// Check if this job has already been completed -- or perhaps already been started
		// with a temporary file having been written to the file system -- so that this
		// job can be skipped.  Not all JobQueens will want to allow already-performed jobs
		// to be skipped, in which case, they can just categorically return false to this query.
		if ( job_queen_->has_job_previously_been_output( larval_job ) ) {
			job_queen_->note_job_completed_and_track( larval_job, jd3_job_previously_executed, 0 );
			job_extractor->note_job_no_longer_running( larval_job->job_index() );
			continue;
		}
		utility::vector1< JobResultCOP > input_job_results =
			construct_job_result_input_list( larval_job );

		// OK -- now turn the larval job into a full-fledged Job
		JobOP mature_job = job_queen_->mature_larval_job( larval_job, input_job_results );
		CompletedJobOutput job_output = run_mature_job( larval_job, mature_job );

		// Inform the JobQueen that the job has run and give her the JobStatus
		job_queen_->note_job_completed_and_track(
			larval_job, job_output.status, job_output.job_results.size() );
		// deliver the job summaries one at a time
		for ( ResultIndex ii( 1 ); ii <= job_output.job_results.size(); ++ii ) {
			job_queen_->completed_job_summary( larval_job, ii, job_output.job_results[ ii ].first );
			job_results_[ std::make_pair( larval_job->job_index(), ii ) ] =
				std::make_pair( larval_job, job_output.job_results[ ii ].second );
		}
		// now inform the job extractor the job is complete; this will perhaps
		// trigger a request for more jobs or more job-digraph nodes from
		// the JobQueen
		job_extractor->note_job_no_longer_running( larval_job->job_index() );

		potentially_output_some_job_results();
		potentially_discard_some_job_results();
	}

	job_queen_->flush();
}

utility::vector1< JobResultCOP >
VanillaJobDistributor::construct_job_result_input_list( LarvalJobCOP larval_job )
{
	// construct the (possibly empty) list of job results needed to mature this larval job
	utility::vector1< JobResultCOP > input_job_results;
	input_job_results.reserve( larval_job->input_job_result_indices().size() );
	for ( auto const &result_id : larval_job->input_job_result_indices() ) {
		auto result_iter = job_results_.find( result_id );
		if ( result_iter == job_results_.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to retrieve job result (" +
				utility::to_string( result_id.first ) + ", " + utility::to_string( result_id.second ) +
				") requested by LarvalJob " + larval_job->job_tag() );
		}
		input_job_results.push_back( result_iter->second.second );
	}
	return input_job_results;
}

CompletedJobOutput
VanillaJobDistributor::run_mature_job(
	LarvalJobOP larval_job,
	JobOP mature_job
)
{
	CompletedJobOutput job_output;
	// Run the job!
	try {
		job_output = mature_job->run();
	} catch ( utility::excn::Exception & e ) {
		// An exception thrown by this job.  Inform the JobQueen that it's a badie.
		TR.Error << "Job " << larval_job->job_index() << " named " << larval_job->job_tag() << " threw an exception:\n" << std::endl;
		e.display();

		job_queen_->note_job_completed_and_track( larval_job, jd3_job_status_failed_w_exception, 0 );
	} catch( std::exception const & e ) {
		// An exception thrown by this job.  Inform the JobQueen that it's a badie.
		TR.Error << "Job " << larval_job->job_index() << " named " << larval_job->job_tag() << " threw an exception:\n" << std::endl;
		std::cerr << "        error message: " << e.what() << std::endl;
		utility::save_crash_report( e.what(), typeid( e ).name() );

		job_queen_->note_job_completed_and_track( larval_job, jd3_job_status_failed_w_exception, 0 );
	} catch ( ... ) {
		// An exception thrown by this job.  Inform the JobQueen that it's a badie.
		TR.Error << "Job " << larval_job->job_index() << " named " << larval_job->job_tag() << " threw an unrecognized exception."
			<< std::endl;
		job_queen_->note_job_completed_and_track( larval_job, jd3_job_status_failed_w_exception, 0 );
	}
	return job_output;
}

void
VanillaJobDistributor::potentially_output_some_job_results()
{
	// ask the job queen if she wants to output any results
	std::list< output::OutputSpecificationOP > jobs_to_output =
		job_queen_->jobs_that_should_be_output();
	for ( auto const &output_spec : jobs_to_output ) {
		JobResultID result_id = output_spec->result_id();
		auto result_iter = job_results_.find( result_id );
		if ( result_iter == job_results_.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to retrieve job result (" +
				utility::to_string( result_id.first )  + ", " + utility::to_string( result_id.second ) +
				") for outputting as requested by the JobQeen. Has this job already been output?" );
		}
		std::pair< LarvalJobOP, JobResultOP > const & job_to_output = result_iter->second;
		output::ResultOutputterOP outputter = job_queen_->result_outputter( *output_spec );
		outputter->write_output( *output_spec, *job_to_output.second );

		// and now delete storage for this job result -- once a JobQueen has requested a job for
		// output, then the JobDistributor will discard it.
		job_results_.erase( result_iter );
	}
}

void
VanillaJobDistributor::potentially_discard_some_job_results()
{
	// ask the job queen if she wants to discard any results
	JobResultIDList jobs_to_discard = job_queen_->job_results_that_should_be_discarded();
	for ( auto const &result_id : jobs_to_discard ) {
		auto result_iter = job_results_.find( result_id );

		if ( result_iter == job_results_.end() ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to retrieve job result (" +
				utility::to_string( result_id.first ) + ", " + utility::to_string( result_id.second ) +
				+ ") for discarding as requested by the JobQeen. Has this job already been output?" );
		}
		job_results_.erase( result_iter );
	}
}

} // namespace job_distributors
} // namespace jd3
} // namespace protocols
