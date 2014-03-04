// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobDistributor.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Base class
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com
/// @author Mike Tyka
/// @author Oliver Lange

// Unit headers
#include <protocols/jd2/JobDistributor.hh>

// Package headers
#include <protocols/jd2/BatchJobInputter.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Parser.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/evaluation/TimeEvaluator.hh>
#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#include <protocols/moves/PyMolMover.hh>

// Utility headers
#include <utility/thread/threadsafe_creation.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>

#ifdef BOINC
#include <protocols/boinc/boinc.hh>
#endif

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#ifndef __native_client__
#include <csignal>
#endif

#ifdef WIN32
#include <iterator>
#endif

namespace protocols {
namespace jd2 {
static basic::Tracer tr("protocols.jd2.JobDistributor");
} //jd2
} //protocols

//// APL Disalbing this code for the time being.  I need to find out who uses it!
////multithreaded case requires special pointers
//#ifdef MULTITHREADED
//
//#include <boost/thread/tss.hpp>
//
//namespace protocols
//{
//	namespace jd2
//	{
//		boost::thread_specific_pointer< JobDistributor > jd_ptr;
//
//		JobDistributor *
//		JobDistributor::get_instance()
//		{
//			if ( jd_ptr.get() == 0 )
//			{
//				jd_ptr.reset( JobDistributorFactory::create_job_distributor() );
//			}
//			return jd_ptr.get();
//		}
//
//	} //jd2
//} //protocols
//
////non-multithreaded case behaves like a singleton
//#else

namespace protocols {
namespace jd2 {

JobDistributor * JobDistributor::instance_ = 0; //this pointer starts null

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex JobDistributor::singleton_mutex_;

std::mutex & JobDistributor::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
JobDistributor * JobDistributor::get_instance()
{
	boost::function< JobDistributor * () > creator = boost::bind( &JobDistributor::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

JobDistributor *
JobDistributor::create_singleton_instance()
{
	return JobDistributorFactory::create_job_distributor();
}


} //jd2
} //protocols

namespace protocols
{
namespace jd2
{

JobDistributor::JobDistributor() :
		//	job_inputter_( JobDistributorFactory::create_job_inputter() ),
		// non-NULL starting state for this pointer; this makes calls to the
		// JobDistributor safe even when not inside go() (of course you will get a
		// stupid object, but at least it won't segfault).  This object deliberately
		// goes away once it's not used.
		current_job_(JD2_BOGUS_JOB->copy_without_output()), current_job_id_(0), last_completed_job_(0), current_batch_id_(
				0)
{
	init_jd();
}

JobDistributor::JobDistributor(bool empty) :
		//	job_inputter_( JobDistributorFactory::create_job_inputter() ),
		// non-NULL starting state for this pointer; this makes calls to the
		// JobDistributor safe even when not inside go() (of course you will get a
		// stupid object, but at least it won't segfault).  This object deliberately
		// goes away once it's not used.
		current_job_(JD2_BOGUS_JOB->copy_without_output()), current_job_id_(0), last_completed_job_(0), current_batch_id_(
				0)
{
	if (!empty)
	{
		init_jd();
	}
	else
	{
		job_inputter_ = NULL;
		job_outputter_ = new NoOutputJobOutputter;
		parser_ = JobDistributorFactory::create_parser();
	}
}

void JobDistributor::init_jd()
{
	instance_ = this; //important so that calls to get_instance in JobInputters or JobOutputters don't lead to a infinite recursion

	// are there batches?
	populate_batch_list_from_cmd();
	if (batches_.size() > 0)
	{
		tr.Debug << "batches present... " << std::endl;
		current_batch_id_ = 1;
		job_inputter_ = new BatchJobInputter(batches_[1]);
	}
	else
	{ //no batches...
		try{
			job_inputter_ = JobDistributorFactory::create_job_inputter();
		} catch (utility::excn::EXCN_Base & excn) {
			basic::Error()
				<< "ERROR: Exception caught by JobDistributor while trying to initialize the JobInputter of type '"
				<< JobInputter::job_inputter_input_source_to_string(
					job_inputter_->input_source())
				<< "'" << std::endl;
			basic::Error()
				<< excn << std::endl;
			utility_exit();
		}
	}

	// get jobs
	try {
		job_inputter_->fill_jobs(jobs_);
	} catch (utility::excn::EXCN_Base & excn) {
			basic::Error()
				<< "ERROR: Exception caught by JobDistributor while trying to fill the input jobs with JobInputter of type type '"
				<< JobInputter::job_inputter_input_source_to_string(
					job_inputter_->input_source())
				<< "'" << std::endl;
			basic::Error()
				<< excn << std::endl;
			utility_exit();
		}

	// have to initialize these AFTER BatchJobInputter->fill_jobs since a new batch might change options
	job_outputter_ = JobDistributorFactory::create_job_outputter();
	parser_ = JobDistributorFactory::create_parser();
}

/// @details read -run:batches and put it into batches_ vector.
void JobDistributor::populate_batch_list_from_cmd()
{
	if (basic::options::option[basic::options::OptionKeys::run::batches].user())
	{
		//typedef utility::vector1< utility::file::FileName >::const_iterator iterator
		utility::vector1<utility::file::FileName> const& fns(
				basic::options::option[basic::options::OptionKeys::run::batches]);
		std::copy(fns.begin(), fns.end(), std::back_inserter(batches_));
	}
}

/// @details restart job-distribution from beginning -- useful if you need a second pass over decoys...
void JobDistributor::restart()
{
	jobs_.clear();
	current_job_id_ = 0;
	last_completed_job_ = 0;
	current_job_ = JD2_BOGUS_JOB->copy_without_output();
	current_batch_id_ = 0;
	init_jd();
}

///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
JobDistributor::~JobDistributor()
{
}

void JobDistributor::go(protocols::moves::MoverOP mover)
{
	go_main(mover);
}

void JobDistributor::go(protocols::moves::MoverOP mover, JobOutputterOP jo)
{
	job_outputter_ = jo;
	go(mover);
}

void JobDistributor::go_main(protocols::moves::MoverOP mover)
{
	using namespace basic::options;
	time_t const allstarttime = time(NULL);
	core::Size tried_jobs(0); //did we try any jobs?

	protocols::moves::MoverOP mover_copy(mover);
	std::string last_inner_job_tag, last_output_tag;
	core::Size last_batch_id = 0; //this will trigger a mover->fresh_instance if we run with batches
	core::Size retries_this_job(0);
	bool first_job(true);

	bool using_parser(false);
	if (parser_)
	{ //if not NULL, we have a parser
		using_parser = true;
		tr.Info
				<< "Parser is present.  Input mover will be overwritten with whatever the parser creates."
				<< std::endl;
	}PROF_START( basic::JD2);
	core::pose::Pose pose;
#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif
	while (obtain_new_job())
	{
		++tried_jobs; //yes, we tried at least one job

		//timing information
		time_t const jobstarttime = time(NULL);
		core::Size const elapsedtime(jobstarttime - allstarttime);

		if ((option[OptionKeys::run::maxruntime].user())
				&& (option[OptionKeys::run::maxruntime]() > 0)
				&& (option[OptionKeys::run::maxruntime]() < int(elapsedtime)))
		{

			basic::Error() << "Run terminating because runtime of "
					<< elapsedtime << " s exceeded maxruntime of "
					<< option[OptionKeys::run::maxruntime]() << " s "
					<< std::endl;
			break; //let it clean up in case there's useful prof information or something
		}

		// setup profiling
		evaluation::TimeEvaluatorOP run_time(NULL);
		if (!option[OptionKeys::run::no_prof_info_in_silentout])
		{
			job_outputter_->add_evaluation(run_time =
					new evaluation::TimeEvaluator); //just don't use this in integration tests!
		}

		tr.Debug << "Starting job " << job_outputter_->output_name(current_job_)
				<< std::endl; //x seconds?

		//Get a copy of the starting pose - this must be done early because the pose is read in on first use, we need to
		//guaruntee it's been read in before the Parser gets a stab at it
		pose.data().clear();
		try
		{
			// Can we add the PyMOL mover here?
			if (option[OptionKeys::run::show_simulation_in_pymol].user()
					&& option[OptionKeys::run::show_simulation_in_pymol].value()
							> 0.0)
			{
				moves::AddPyMolObserver(
						pose,
						option[OptionKeys::run::keep_pymol_simulation_history](),
						option[OptionKeys::run::show_simulation_in_pymol].value());
			}
			job_inputter_->pose_from_job(pose, current_job_);

#ifdef BOINC_GRAPHICS
			// attach boinc graphics pose observer
			// do it here because pose_from_job may replace the pose conformation
			protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

		} catch (utility::excn::EXCN_Base& excn)
		{
			basic::Error()
					<< "ERROR: Exception caught by JobDistributor while trying to get pose from job "
					<< "'" << job_outputter_->output_name(current_job_) << "'" << std::endl
					<< excn << std::endl;
			basic::Error()
					<< "Treating failure as bad input; canceling similar jobs"
					<< std::endl;
			remove_bad_inputs_from_job_list();
			job_failed(pose, false);
			pose = core::pose::Pose();
			continue;
		}

		//These if statements determine whether a new creation of the mover is appropriate
		bool reinitialize_new_input(false);
		bool new_input(false);
		if (current_job_->input_tag() != last_inner_job_tag)
		{
			//this means we have just changed inputs - the next pdb on -l, etc
			tr.Debug << "new input detected, is: " << current_job_->input_tag()
					<< ", last was: " << last_inner_job_tag << std::endl;
			last_inner_job_tag = current_job_->input_tag();
			new_input = true;

			//do we need to reinitialize because of the new input? - yes if mover says, or cmdline says
			if (mover
					&& (mover->reinitialize_for_new_input()
							|| option[OptionKeys::run::reinitialize_mover_for_new_input]))
			{
				reinitialize_new_input = true;
			} //if we need to reinitialize

		} //if the input pose is about to change

		// Are we on a new output structure? (Or are we repeating a failed job?)
		if( job_outputter_->output_name(current_job_) != last_output_tag ) {
			last_output_tag = job_outputter_->output_name(current_job_);
			retries_this_job = 0;
		}

		if (option[OptionKeys::jd2::delete_old_poses].user())
		{
			//to improve jd2 memory performance, we will delete the last
			//input's starting pose. (Previous to this, jd2 never deleted
			//input poses and would accumulate memory over large input sets
			//- not a memory leak but certainly a nasty spot in the
			//basement.) SML 8/7/09

			//This was applied in r32237 but it had problems with special
			//uses of the job distributor and it was reverted. This should
			//probably be applied by default, once the issues with the
			//special uses are worked out.

			if (!first_job && last_completed_job_ != 0)
			{
				tr.Debug << "deleting pose from job " << last_completed_job_ <<std::endl;
				jobs_[last_completed_job_]->inner_job_nonconst()->set_pose(
						NULL);
			}
		}
		//delete pointer to pose of last input; if that was last pointer
		//to pose (and it should have been) this will free the memory


		if (current_batch_id() != last_batch_id)
		{
			tr.Debug << "new batch detected: get fresh instance from mover"
					<< std::endl;
			new_input = true;
			reinitialize_new_input = true;
			last_batch_id = current_batch_id();
		}

		//for regular movers, reinitialize if desired
		if (!using_parser
				&& (reinitialize_new_input || mover->reinitialize_for_each_job()
						|| option[OptionKeys::run::reinitialize_mover_for_each_job]))
		{
			tr.Trace << "reinitializing the mover" << std::endl;
			PROF_STOP( basic::JD2);
			PROF_START( basic::JD2_INIT_MOVER);
			mover_copy = mover->fresh_instance();
			PROF_STOP( basic::JD2_INIT_MOVER);
			PROF_START( basic::JD2);
		}
		else if (using_parser)
		{ //call the parser
			tr.Trace << "Allowing the Parser to create a new Mover if desired"
					<< std::endl;
			try
			{
				parser_->generate_mover_from_job(current_job_, mover_copy,
						new_input);
			} catch (utility::excn::EXCN_Base& excn)
			{
				basic::Error()
					<< "ERROR: Exception caught by JobDistributor while trying to get pose from job '"
					<< job_outputter_->output_name(current_job_) << "'" << std::endl
					<< excn
					<< std::endl;
				basic::Error()
					<< "Treating failure as bad input; canceling similar jobs"
					<< std::endl;
				remove_bad_inputs_from_job_list();
				job_failed(pose, false);
				continue;
			}
			//the Parser might have modified the starting pose (with constraints) - so we'll refresh our copy
			job_inputter_->pose_from_job(pose, current_job_);

#ifdef BOINC_GRAPHICS
			// attach boinc graphics pose observer
			// do it here because pose_from_job may replace the pose conformation
			protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif
#ifdef GL_GRAPHICS
			//nonboinc viewer
			protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" );
#endif

		}
		else
		{
			tr.Trace << "not reinitializing mover" << std::endl;
			//mover_copy = mover; //This breaks when reinitializing only on new input, because non-reinitializing cycles will
			//revert to the wrong place.  If a mover_copy = something is desireable for all options, we need a second
			//mover_copy.
		}

		//use the mover
		mover_copy->reset_status();
		// clear old string info from previous apply calls
		mover_copy->clear_info();
		if (run_time)
			run_time->reset(); //reset clock of TimeEvaluator

		// notify JobOutputter of starting pose, for comparison purposes and/or as interface for initializing evaluators. (Currently does nothing in the base class.)
		job_outputter_->starting_pose(pose);

		protocols::moves::MoverStatus status;
		PROF_STOP( basic::JD2);
		try
		{
			if (basic::options::option[basic::options::OptionKeys::out::std_IO_exit_error_code]()
					> 0)
				std::cout.exceptions(std::ios_base::badbit);

			tr.Debug << "run mover... " << std::endl;
			mover_copy->set_current_tag(
					job_outputter_->output_name(current_job_));
			mover_copy->apply(pose);
			status = mover_copy->get_last_move_status();
			// Job collects (optional) string info from the mover.
			// This info may be output later by a JobOutputter.
			current_job_->add_strings(mover_copy->info());

		} catch (std::ios_base::failure& ex)
		{
			std::cerr << "std::IO error detected... exiting..." << std::endl; // We can not longer use Tracer's at this point
			std::exit(
					basic::options::option[basic::options::OptionKeys::out::std_IO_exit_error_code]()); // Using pure exit instead of utility_exit_with_status to avoid recursion

		} catch (utility::excn::EXCN_BadInput& excn)
		{
			tr.Error
					<< "\n\n[ERROR] Exception caught by JobDistributor for job "
					<< job_outputter_->output_name(current_job_) << excn
					<< std::endl;
			status = protocols::moves::FAIL_BAD_INPUT;

		} catch (utility::excn::EXCN_Base& excn)
		{
			tr.Error
					<< "\n\n[ERROR] Exception caught by JobDistributor for job "
					<< job_outputter_->output_name(current_job_) << excn
					<< std::endl;
			status = protocols::moves::FAIL_DO_NOT_RETRY;
		}
		std::cout.exceptions(std::ios_base::goodbit); // Disabling std::IO exceptions

		PROF_START( basic::JD2);
		core::Size jobtime(time(NULL) - jobstarttime);

		begin_critical_section();
		PROF_START( basic::JD2_OUTPUT);
		// check cases: SUCCESS, FAIL_RETRY, FAIL_DO_NOT_RETRY, FAIL_BAD_INPUT
		if (status == protocols::moves::MS_SUCCESS)
		{
			last_completed_job_ = current_job_id_;
			job_succeeded(pose, jobtime);
			//			tr.Info << job_outputter_->output_name( current_job_ ) << " reported success in " << jobtime << " seconds" << std::endl;
		}
		else if (status == protocols::moves::FAIL_RETRY)
		{
			using namespace basic::options::OptionKeys::jd2;
			++retries_this_job;
			if (option[ntrials].user()
					&& (retries_this_job >= (core::Size) option[ntrials].value()))
			{
				//this represents too many FAIL_RETRY - we will roll over into FAIL_DO_NOT_RETRY
				tr.Warning << job_outputter_->output_name(current_job_)
						<< " reported failure " << retries_this_job
						<< " times and will no longer retry (permanent failure)"
						<< std::endl;
				job_failed(pose, false /* will not retry */);
			}
			else
			{
				mark_current_job_id_for_repetition();
				tr.Warning << job_outputter_->output_name(current_job_)
						<< " reported failure and will retry" << std::endl;
				job_failed(pose, true /* will retry */);
			}
		}
		else if (status == protocols::moves::FAIL_DO_NOT_RETRY)
		{
			tr.Warning << job_outputter_->output_name(current_job_)
					<< " reported failure and will NOT retry" << std::endl;
			job_failed(pose, false /* will not retry */);
		}
		else if (status == protocols::moves::FAIL_BAD_INPUT)
		{
			tr.Warning << job_outputter_->output_name(current_job_)
					<< " reported that its input was bad and will not retry"
					<< std::endl;
			remove_bad_inputs_from_job_list();
			job_failed(pose, false /*will not retry */);
		}
		end_critical_section();
		current_job_finished();
		PROF_STOP( basic::JD2_OUTPUT);
		basic::prof_show();
		first_job = false; //we've finished one by now, and are no longer on the first job edge case
	}PROF_STOP( basic::JD2);
	note_all_jobs_finished();
	if (batches_.size())
	{
		tr.Info << jobs_.size() << " jobs in last batch... in total ";
	}
	else
	{
		tr.Info << jobs_.size() << " jobs considered, ";
	}
	tr.Info << tried_jobs << " jobs attempted in "
			<< (time(NULL) - allstarttime) << " seconds" << std::endl;
	if (tried_jobs == 0)
		tr.Info << "no jobs were attempted, did you forget to pass -overwrite?"
				<< std::endl;
	job_outputter_->flush(); //This call forces out any unprinted data
	basic::prof_show();
}

JobOP JobDistributor::current_job() const
{
	return current_job_;
}

std::string JobDistributor::current_output_name() const
{
	return job_outputter()->output_name(current_job());
}

JobOutputterOP JobDistributor::job_outputter() const
{
	return job_outputter_;
}

/// @brief The input source for the current JobInputter.
JobInputterInputSource::Enum JobDistributor::job_inputter_input_source() const
{
	return job_inputter_->input_source();
}

bool JobDistributor::obtain_new_job(bool reconsider_current_job)
{
	if (reconsider_current_job)
		--current_job_id_;

	if (batches_.size() == 0
			|| get_current_batch() != BatchJobInputter::BOGUS_BATCH_ID)
	{ //batches can be cancelled during computation
		current_job_id_ = get_new_job_id(); //if no batches are present, or current batch still valid
	}
	else
	{
		current_job_id_ = 0; //batch got cancelled... jump to end of batch....
	}

	if (current_job_id_ == 0)
	{
		if (next_batch())
		{ //query if there is a new batch to run after this one has finished
			current_job_id_ = 0;
			return obtain_new_job(); //set to first job of new batch... --- if batch is already computed fully this migh call next_batch() !
		}
		return false;
	}
	else if (current_job_id_ <= jobs_.size())
	{
		current_job_ = jobs_[current_job_id_];
		return true;
	}
	else
	{
		utility_exit_with_message(
				"JobDistributor: nonexistent job returned in obtain_new_job()");
		return false;
	}
}

void JobDistributor::job_succeeded(core::pose::Pose & pose, core::Real run_time)
{
	job_outputter_->final_pose(current_job_, pose);
	//current_job_->set_completed();
	mark_job_as_completed(current_job_id_, run_time);
	return;
}

/// @details no-op implementation in the base class
void JobDistributor::job_failed(core::pose::Pose & /*pose*/,
		bool /*will_retry*/)
{
}

void JobDistributor::mark_job_as_completed(core::Size job_id,
		core::Real run_time)
{
	jobs_[job_id]->set_completed();
	tr.Info << job_outputter_->output_name(jobs_[job_id])
			<< " reported success in " << run_time << " seconds" << std::endl;
	//	tr.Info << "completed job: " << job_outputter_->output_name( jobs_[ job_id ] ) << std::endl;
}

void JobDistributor::mark_job_as_bad(core::Size job_id)
{
	jobs_[job_id]->set_bad();
}

void JobDistributor::remove_bad_inputs_from_job_list()
{
}

void JobDistributor::current_job_finished()
{
}

void JobDistributor::note_all_jobs_finished()
{
}

//This next line prevents accumulation of state within the Job object - should it be within another function?
void JobDistributor::clear_current_job_output()
{
	jobs_[current_job_id_] = current_job_->copy_without_output(); //is this unsafe?  should be its own function? MT: It is now!
}

void JobDistributor::begin_critical_section()
{
}

void JobDistributor::end_critical_section()
{
}

//////////////////////protected accessor functions////////////////////
core::Size JobDistributor::current_job_id() const
{
	return current_job_id_;
}

Jobs const &
JobDistributor::get_jobs() const
{
	return jobs_;
}

// Jobs &
// JobDistributor::get_jobs()  { return jobs_; }

JobInputterOP JobDistributor::job_inputter() const
{
	return job_inputter_;
}

void JobDistributor::mpi_finalize(bool)
{
	//dummy default implementation
}

ParserOP JobDistributor::parser() const
{
	return parser_;
}

/////////////////////////batch stuff //////////////////////////
std::string JobDistributor::get_current_batch() const
{
	if (current_batch_id_ && batches_.size() > 0
			&& current_batch_id_ <= batches_.size())
	{
		return batches_[current_batch_id_];
	}
	else
	{
		return BatchJobInputter::BOGUS_BATCH_ID;
	}
}

void JobDistributor::set_batch_id(core::Size setting)
{
	if (current_batch_id_ == setting)
		return;
	current_batch_id_ = setting;
	if (current_batch_id_ > batches_.size())
		batch_underflow();
	if (current_batch_id_ > batches_.size())
	{
		tr.Error << "[ERROR] illegeal attempt to set batch_id to " << setting
				<< " where we have only " << batches_.size() << " batches"
				<< std::endl;
		utility_exit_with_message("wrong batch_id");
	}
	load_new_batch();
}

bool JobDistributor::next_batch()
{
	++current_batch_id_;

	if (current_batch_id_ > batches_.size())
		batch_underflow();
	if (current_batch_id_ > batches_.size())
	{ //still no new batches.
		tr.Info << "no more batches to process... " << std::endl;
		return false;
	}

	//skip BOGUS_BATCHES ..
	while (current_batch_id_ <= batches_.size()
			&& get_current_batch() == BatchJobInputter::BOGUS_BATCH_ID)
		++current_batch_id_;

	//if ended on BOGUS_BATCH
	if (get_current_batch() == BatchJobInputter::BOGUS_BATCH_ID)
	{
		tr.Trace << "last batch is CANCELLED: run next_batch()" << std::endl;
		return next_batch();
	}

	runtime_assert( current_batch_id_ <= batches_.size());
	load_new_batch();
	return true;
}

/// @detail add new batch to end of batch list... this might be called asynchronous... ie. while we are still in the middle of
/// a current batch, or while we are in non-batch mode
void JobDistributor::add_batch(std::string const& batch, core::Size id)
{
	while (id > batches_.size())
	{
		batches_.push_back(BatchJobInputter::BOGUS_BATCH_ID);
	}
	if (id > 0)
	{
		batches_[id] = batch;
	}
	else
	{
		batches_.push_back(batch);
	}
}

/// @detail restart JobDistributor with a new batch the BatchJobInputter loads
/// new flags and sets global options after this we reload Factory dependent
/// Objects (e.g., JobInputter and JobOutputter )
void JobDistributor::load_new_batch()
{
	runtime_assert( current_batch_id_ <= batches_.size());
	//paranoid

	jobs_.clear();
	current_job_id_ = 0;
	current_job_ = JD2_BOGUS_JOB->copy_without_output();

	//remaking job_inputter has the advantage, that we will also get one, if
	//this is the first batch!
	job_inputter_ = NULL; //triggers destructor --> restores options
	tr.Info << "start batch " << batches_[current_batch_id_] << std::endl;
	job_inputter_ = new BatchJobInputter(batches_[current_batch_id_]);

	job_inputter_->fill_jobs(jobs_);
	// have to initialize these AFTER BatchJobInputter->fill_jobs an new batch
	// might change options

	// should we copy the old-evaluators --- of course only if new options don't
	// change to a different silent-file or score-file
	// evaluation::PoseEvaluators const& evaluations( job_outputter->evaluators() );

	job_outputter_ = JobDistributorFactory::create_job_outputter();
	parser_ = JobDistributorFactory::create_parser();
}

/// @brief Movers (or derived classes) may ask for the JobOutputter
void JobDistributor::set_job_outputter(const JobOutputterOP &new_job_outputter)
{
	job_outputter_ = new_job_outputter;
}

// // JobDistributorDestroyer functions
// JobDistributorDestroyer::JobDistributorDestroyer(JobDistributor* jd) {
// 	jd_ = jd;
// }

// JobDistributorDestroyer::~JobDistributorDestroyer() {
// 	delete jd_;
// }
// void JobDistributorDestroyer::set_job_distributor(JobDistributor* jd) {
// 	jd_ = jd;
// }

/// @details Default callback function for signal handling
void JobDistributor::jd2_signal_handler(int signal_)
{
	std::cout << "Got some signal... It is:" << signal_ << std::endl;
	if (signal_ == SIGINT)
		std::cout << "Ctrl-c was pressed!" << std::endl;
	if (signal_ == SIGABRT)
		std::cout << "Process was aborted!" << std::endl;
	if (signal_ == SIGTERM)
		std::cout << "Process was terminated!" << std::endl;

#ifndef WIN32
	if (signal_ == SIGKILL)
		std::cout << "Process was SIGKILL!" << std::endl;
	if (signal_ == SIGQUIT)
		std::cout << "Process was SIGQUIT!" << std::endl;
#endif

	if (get_instance()->job_outputter_)
		get_instance()->job_outputter_->flush(); //This call forces out any unprinted but finished data
	get_instance()->handle_interrupt();

	//utility_exit_with_status(1);
	std::exit(1); // Using pure exit instead of utility_exit_with_status to avoid recursion when compile with EXIT_THROWS_EXCEPTION
	//SML 3/22/13 - I am no longer sure this is necessary given that all apps are try/catch wrapped and all utility_exit throw exceptions
}

/// @details Setting up callback function that will be call when our process is about to terminate.
// /@details This will allow us to exit propely (clean up in_progress_files/tmp files if any).
void JobDistributor::setup_system_signal_handler(void(*signal_fn)(int))
{
#ifndef __native_client__
  // Soooo many way to kill... wait - there is no special signal for HS with SVD? - lame...
	signal(SIGINT, signal_fn);
	signal(SIGABRT, signal_fn);
	signal(SIGTERM, signal_fn);

#ifndef WIN32
	signal(SIGKILL, signal_fn);
	signal(SIGQUIT, signal_fn);
#endif
#endif

}

/// @details Set signal handler back to default state.
void JobDistributor::remove_system_signal_handler()
{
#ifndef __native_client__
	signal(SIGINT, SIG_DFL);
	signal(SIGABRT, SIG_DFL);
	signal(SIGTERM, SIG_DFL);
#ifndef WIN32
	signal(SIGKILL, SIG_DFL);
	signal(SIGQUIT, SIG_DFL);
#endif
#endif
}

} //jd2
} //protocols
