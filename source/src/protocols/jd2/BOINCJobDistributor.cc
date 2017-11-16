// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/BOINCJobDistributor.cc
/// @brief  implementation of BOINCJobDistributor
/// @author Mike Tyka


// -- IMPORTANT --
// This has to come before boinc.hh or we get this error on VC++
// '_read' : is not a member of 'std::basic_istream<_Elem,_Traits>'
#include <utility/io/izstream.hh>


#ifdef BOINC
	#include <protocols/boinc/boinc.hh>
#endif // BOINC

#include <core/io/silent/util.hh>
#include <protocols/checkpoint/Checkpoint.hh>
#ifdef BOINC
#ifdef USEMPI
Throw a compiler error because MPI and BOINC cannot be used together!
If you got this message, something is wrong with your build settings.
#endif
#endif

// Unit headers
#include <protocols/jd2/BOINCJobDistributor.hh>
// Package headers
#ifdef BOINC
#include <protocols/jd2/JobOutputter.hh>
#endif

#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


static basic::Tracer TR( "protocols.jd2.BOINCJobDistributor" );


namespace protocols {
namespace jd2 {

/// @details constructor.  Notice it calls the parent class!
BOINCJobDistributor::BOINCJobDistributor() :
	ShuffleFileSystemJobDistributor(),
	total_completed_nstruct_( 0 )
{
#ifdef BOINC
	checkpoint_read();
	protocols::boinc::Boinc::worker_startup();
#endif
}

/// @brief dtor
///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
BOINCJobDistributor::~BOINCJobDistributor() = default;

void
BOINCJobDistributor::go( protocols::moves::MoverOP mover )
{

#ifdef BOINC
	JobsContainer const & jobs( get_jobs() );
	JobOutputterOP outputter = job_outputter();
	// Count completed jobs
	for( core::Size i = 1; i <= jobs.size(); i ++ ){
		if( outputter->job_has_completed( jobs[ i ] ) ) total_completed_nstruct_++;
	}
#endif

	ShuffleFileSystemJobDistributor::go( mover );

	// gzip the output silent files.
	core::io::silent::gzip();

	// ideally these would be called in the dtor but the way we have the singleton pattern set up the dtors don't get
	// called
#ifdef BOINC
	checkpoint_clear();
	protocols::boinc::Boinc::worker_finish_summary( total_completed_nstruct_, total_completed_nstruct_, jobs.size() );
	protocols::boinc::Boinc::worker_shutdown(); // Does not return.
#endif
}

void
BOINCJobDistributor::checkpoint_write()
{
	begin_critical_section();
	static time_t last_chkpt_time = time(nullptr);
	time_t time_now = time(nullptr);
	// Refuse to checkpoint more than once a minute, no matter what BOINC wants.
	// Random number checkpoint files can be large (100k or more uncompressed).
	if ( time_now - last_chkpt_time > 60 ) {
#ifdef BOINC
		// BOINC automatically handles begin/end_critical_section() calls.
		utility::io::ozstream ozs("rng.state.gz");
		numeric::random::rg().saveState(ozs);
		ozs.close();
#endif // BOINC
		protocols::checkpoint::Timer::reset();
		last_chkpt_time = time_now;
	}
	end_critical_section();
}

void
BOINCJobDistributor::checkpoint_read()
{
	begin_critical_section();
#ifdef BOINC
	if( utility::file::file_exists("rng.state.gz") ) {
		utility::io::izstream izs("rng.state.gz");
		numeric::random::rg().restoreState(izs);
		izs.close();
	}
#endif // BOINC
	end_critical_section();
}

void
BOINCJobDistributor::checkpoint_clear()
{
#ifdef BOINC
	if( utility::file::file_exists("rng.state.gz") ) {
		utility::file::file_delete("rng.state.gz");
	}
#endif // BOINC
}

/// @brief dummy for master/slave version
core::Size
BOINCJobDistributor::get_new_job_id()
{
#ifdef BOINC
	if (protocols::boinc::Boinc::worker_is_finished( total_completed_nstruct_ )) return 0; // we're done no matter what nstruct syays
#endif
	return ShuffleFileSystemJobDistributor::get_new_job_id();
}

void
BOINCJobDistributor::mark_current_job_id_for_repetition()
{
	// do nothing - no repetitions allowed. Behave as if job had succeeded
}

void BOINCJobDistributor::job_failed( core::pose::Pose& pose, bool /*will_retry*/ ) {
	using namespace basic::options;

	current_job()->set_status_prefix("FAILURE");
	job_succeeded( pose, 0, "" ); //i.e., job_outputter_->final_pose( current_job(), pose );
}

void
BOINCJobDistributor::job_succeeded(core::pose::Pose & pose, core::Real run_time, std::string const & tag)
{
	FileSystemJobDistributor::job_succeeded( pose, run_time, tag );
#ifdef BOINC
	checkpoint_write();
	total_completed_nstruct_++;
#endif
}

void BOINCJobDistributor::begin_critical_section() {
#ifdef BOINC
	boinc_begin_critical_section();
#endif // BOINC
}

void BOINCJobDistributor::end_critical_section() {
#ifdef BOINC
	boinc_end_critical_section();
#endif // BOINC
}


}//jd2
}//protocols


