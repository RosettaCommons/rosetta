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

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/jd3.OptionKeys.gen.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>
#include <ctime>
#include <fstream>

namespace protocols {
namespace jd3 {
namespace job_distributors {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd3.VanillaJobDistributor" );

VanillaJobDistributor::VanillaJobDistributor()
{}

VanillaJobDistributor::~VanillaJobDistributor() {}

void
VanillaJobDistributor::go( JobQueenOP queen ) {
	using core::Size;
	using namespace utility::graph;
	typedef std::list< Size > SizeList;

	if ( basic::options::option[ basic::options::OptionKeys::jd3::job_definition_schema ].user() ) {
		std::string xsd = queen->job_definition_xsd();
		std::string outfile_name = basic::options::option[ basic::options::OptionKeys::jd3::job_definition_schema ]();
		std::ofstream ofs( outfile_name.c_str() );
		ofs << xsd;
		return;
	}

	job_queen_ = queen;
	job_dag_ = job_queen_->initial_job_dag();


	while ( true ) {
		std::pair< SizeList, bool > job_total_order_pair = topological_sort( *job_dag_ );
		if ( ! job_total_order_pair.second ) {
			throw utility::excn::EXCN_Msg_Exception( "JobQueen produced a non-DAG job graph." );
		}
		SizeList job_total_order = job_total_order_pair.first;
		for ( SizeList::const_iterator node_iter = job_total_order.begin();
				node_iter != job_total_order.end(); ++node_iter ) {
			Size job_node = *node_iter;
			run_jobs_for_dag_node( job_node );
		}

		// ok -- in here, the VJD asks the JQ to update the JobDigraph by adding
		// new nodes and edges to those new nodes.  If the JQ does not add any new nodes,
		// then the JD will exit this loop
		JobDigraphUpdater updater( job_dag_ );
		job_queen_->update_job_dag( updater );
		if ( job_dag_->num_nodes() == updater.orig_num_nodes() ) {
			// no new nodes in the graph, therefore, we are done
			break;
		}
	}

}

void
VanillaJobDistributor::run_jobs_for_dag_node( core::Size job_node )
{
	while ( true ) {
		LarvalJobs jobs_for_node = job_queen_->determine_job_list( job_node, 1000 );
		if ( jobs_for_node.empty() ) break;

		for ( LarvalJobs::const_iterator job_iter = jobs_for_node.begin(); job_iter != jobs_for_node.end(); ++job_iter ) {
			LarvalJobOP larval_job = *job_iter;
			// Check if this job has already been completed -- or perhaps already been started
			// with a temporary file having been written to the file system -- so that this
			// job can be skipped.  Not all JobQueens will want to allow already-performed jobs
			// to be skipped, in which case, they can just categorically return false to this query.
			if ( job_queen_->has_job_completed( larval_job ) ) {
				job_queen_->note_job_completed( larval_job, jd3_job_previously_executed );
				continue;
			}

			utility::vector1< JobResultCOP > input_job_results =
				construct_job_result_input_list( larval_job );

			// OK -- now turn the larval job into a full-fledged Job
			JobOP mature_job = job_queen_->mature_larval_job( larval_job, input_job_results );

			CompletedJobOutput job_output = run_mature_job( larval_job, mature_job );

			// run_mature_job returns a null-pointing job summary if the job failed.
			if ( job_output.first ) {
				// Inform the JobQueen that the job has run and give her the JobSummary
				job_queen_->note_job_completed( larval_job, job_output.first->status() );
				job_queen_->completed_job_summary( larval_job, job_output.first );

				job_results_[ larval_job->job_index() ] = std::make_pair( larval_job, job_output.second );
			}

			potentially_output_some_job_results();
			potentially_discard_some_job_results();

		}
	}

}

utility::vector1< JobResultCOP >
VanillaJobDistributor::construct_job_result_input_list( LarvalJobCOP larval_job )
{
	// construct the (possibly empty) list of job results needed to mature this larval job
	utility::vector1< JobResultCOP > input_job_results;
	input_job_results.reserve( larval_job->input_job_result_indices().size() );
	for ( utility::vector1< core::Size >::const_iterator result_index_iter = larval_job->input_job_result_indices().begin();
			result_index_iter != larval_job->input_job_result_indices().end(); ++result_index_iter ) {
		if ( job_results_.count( *result_index_iter ) == 0 ) {
			throw utility::excn::EXCN_Msg_Exception( "Failed to retrieve job result " +
				utility::to_string( *result_index_iter ) + " requested by LarvalJob " + larval_job->job_tag() );
		}
		input_job_results.push_back( job_results_[ * result_index_iter ].second );
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
	} catch ( utility::excn::EXCN_Base & e ) {
		// An exception thrown by this job.  Inform the JobQueen that it's a badie.
		TR.Error << "Job " << larval_job->job_index() << " named " << larval_job->job_tag() << " threw an exception:\n" <<
			e.msg() << std::endl;
		job_queen_->note_job_completed( larval_job, jd3_job_status_failed_w_exception );
	} catch ( ... ) {
		// An exception thrown by this job.  Inform the JobQueen that it's a badie.
		job_queen_->note_job_completed( larval_job, jd3_job_status_failed_w_exception );
		TR.Error << "Job " << larval_job->job_index() << " named " << larval_job->job_tag() << " threw an unrecognized exception."
			<< std::endl;
	}
	return job_output;
}

void
VanillaJobDistributor::potentially_output_some_job_results()
{
	// ask the job queen if she wants to output any results
	SizeList jobs_to_output = job_queen_->jobs_that_should_be_output();
	for ( SizeList::const_iterator output_iter = jobs_to_output.begin();
			output_iter != jobs_to_output.end(); ++output_iter ) {
		if ( job_results_.count( *output_iter ) == 0 ) {
			throw utility::excn::EXCN_Msg_Exception( "Failed to retrieve job result " + utility::to_string( *output_iter ) +
				" for outputting as requested by the JobQeen. Has this job already been output?" );
		}
		std::pair< LarvalJobOP, JobResultOP > job_to_output = job_results_[ *output_iter ];
		job_queen_->completed_job_result( job_to_output.first, job_to_output.second );
		// and now delete storage for this job result
		job_results_.erase( job_results_.find( *output_iter ));
	}
}

void
VanillaJobDistributor::potentially_discard_some_job_results()
{
	// ask the job queen if she wants to discard any results
	SizeList jobs_to_discard = job_queen_->job_results_that_should_be_discarded();
	for ( SizeList::const_iterator discard_iter = jobs_to_discard.begin();
			discard_iter != jobs_to_discard.end(); ++discard_iter ) {
		if ( job_results_.count( *discard_iter ) == 0 ) {
			throw utility::excn::EXCN_Msg_Exception( "Failed to retrieve job result " + utility::to_string( *discard_iter )
				+ " for discardting as requested by the JobQeen. Has this job already been output?" );
		}
		job_results_.erase( job_results_.find( *discard_iter ));
	}
}

} // namespace job_distributors
} // namespace jd3
} // namespace protocols
