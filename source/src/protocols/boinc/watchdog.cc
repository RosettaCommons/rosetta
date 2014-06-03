// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 7630 $
//  $Date: 2006-03-10 09:37:52 -0800 (Fri, 10 Mar 2006) $
//  $Author: mtyka & rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifdef BOINC
#include <utility/io/izstream.hh>
#include <protocols/boinc/boinc.hh>
#endif

// Unit header
#include <protocols/boinc/watchdog.hh>

#include <iostream>


#ifdef BOINC
#include <protocols/boinc/boinc.hh>
#include <core/io/silent/util.hh>

// Utility Headers
#include <utility/io/ozstream.hh>

#include <utility/basic_sys_util.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <ctime>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/boinc.OptionKeys.gen.hh>


#ifndef _WIN32
#include "pthread.h"
#endif

namespace protocols {
namespace boinc {
namespace watchdog {

// protocols can set this pose as the global bailout - if the watchdog kicks in it will write out *this*
// pose and give it a special label to be identified as the Bailout ( W_xxx )

#ifndef _WIN32
	pthread_mutex_t bailout_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

std::string bailout_silent_structure = "watchdog_failure: \n";
std::string bailout_silent_structure_header = "";

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//  The watchdog thread has a few roles in BOINC builds.
//
//  1. If Rosetta has been going for cpu_run_timeout long than the user's preferred CPU run time,
//     that's way too long. Exit gracefully. Default cpu_run_timeout = 3600*4
//
//  2. Update percentage complete every second or so, for user feedback.
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

	int WATCHDOG_BLINK_TIME = 2; // Watchdog sees if Rosetta is still running every second or so.
	int const PCT_COMPLETE_UPDATE_TIME = 5; // resolution of pct complete status updates in seconds

#ifdef _WIN32
HANDLE watchdogThread;

UINT WINAPI main_watchdog_windows( void* lpParam )
{
	main_watchdog(NULL);
	return 0;
}
#endif

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void
watchdog_start(){
	using namespace basic::options;

	// This is default for
	// BOINC runs -- can turn it off with no_watchdog.
	if ( option[ OptionKeys::boinc::watchdog ] ) {
		std::cerr << "Starting watchdog..." << std::endl;
#ifdef _WIN32
		watchdogThread = (HANDLE)_beginthreadex(
            NULL,              // default security attributes
            0,                 // use default stack size
            main_watchdog_windows,  // thread function
            NULL,              // argument to thread function
            0,                 // use default creation flags
            NULL);             // returns the thread identifier
#else
		pthread_t p_watchdog;
		pthread_create ( &p_watchdog, NULL,
										 &main_watchdog, NULL );
#endif

	}
	return;
}

void
watchdog_finish()
{
	using namespace basic::options;
	// This is default for
	// BOINC runs -- can turn it off with no_watchdog.
	if ( option[ OptionKeys::boinc::watchdog ] ) {
		std::cerr << std::endl;
		std::cerr << "BOINC :: Watchdog shutting down..." << std::endl;

		// should already by set
		protocols::boinc::Boinc::stop_running_worker();
#ifdef _WIN32
		// Wait for the watchdog thread to shutdown before leaving the function.
		WaitForSingleObject(watchdogThread, INFINITE);
		CloseHandle(watchdogThread);
#endif
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////
void
get_the_hell_out(std::string )
{
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	protocols::boinc::Boinc::stop_running_worker();

	boinc_fraction_done(1);
	boinc_begin_critical_section();

	std::cerr << utility::timestamp() << " :: BOINC " << std::endl;
	std::cerr.flush();

	std::string const fullname = option[ out::file::silent ]();

	unsigned int decoy_estimate =  std::max(1,protocols::boinc::Boinc::decoy_count()) + 1;
	// if structures were built, just bail quietly and return those structures.
	if( utility::file::file_exists( fullname )  ){
		std::cerr << "InternalDecoyCount: " << protocols::boinc::Boinc::decoy_count()  << std::endl;

		// Count ?
		// print a fake message to satisfy validator
		protocols::boinc::Boinc::worker_finish_summary( decoy_estimate, decoy_estimate , 2  );

		core::io::silent::gzip();
		boinc_finish(0); // process terminates right here, nothing else needed.
	}

	// or maybe it's already zipped ?
	if( utility::file::file_exists( fullname + ".gz" ) ){
		// in that case just bail
		std::cerr << "Output exists: " << fullname + ".gz" << " Size: " << utility::file::file_size( fullname + ".gz" ) << std::endl;
		std::cerr << "InternalDecoyCount: " << protocols::boinc::Boinc::decoy_count() << " (GZ)" << std::endl;

		// Check that it's readablie
		std::string wholename( fullname + ".gz" );
		utility::io::izstream infile( wholename.c_str() );

		std::cerr << "-----" << std::endl;
		std::cerr << infile << std::endl;
		std::cerr << "-----" << std::endl;

		if (infile.good()) {
			infile.close();
			// print a fake message to satisfy validator
			protocols::boinc::Boinc::worker_finish_summary( decoy_estimate, decoy_estimate , 2  );
			boinc_finish(0); // process terminates right here, nothing else needed.
		}else{
			std::cerr << "Stream information inconsistent." << std::endl;
		}
	}

	// ONLY IF output file does not exist, make a dummy one !

	std::cerr << "Writing W_0000001" << std::endl;

	// otherwise write a watchdig structure - however this tends to fail the
	// validator currently.

	// change Rhiju's original logic a little bit. Instead of creating a blank
	// .out file, we will write a watchdog_failure line so that the validator
	// can take operators accordingly. This will avoid throwing a lot of client
	// errors on BOINC as a blank file can not be gzipped. -- chu

	utility::io::ozstream pdb_out_checkstream( fullname,
			std::ios_base::in|std::ios_base::out );
	utility::io::ozstream pdb_out_stream;
#ifdef BOINC_GRAPHICS
	if (!Boinc::trywait_semaphore()) {
#else
	pthread_mutex_lock(&bailout_mutex);
#endif
	if (!pdb_out_checkstream ) {
		pdb_out_checkstream.close();
		pdb_out_checkstream.clear();
		pdb_out_stream.open( fullname );
		if ( !pdb_out_stream ) {
			std::cout << "Open failed for file: " << fullname << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		pdb_out_stream << bailout_silent_structure_header;
		pdb_out_stream << bailout_silent_structure;
		//pdb_out_stream << "REMARK  " <<  moreinfostring << std::endl;
	} else {
		pdb_out_checkstream.close();
		pdb_out_checkstream.clear();
		pdb_out_stream.open_append( fullname );
		if ( !pdb_out_stream ) {
			std::cout << " Append failed for file: " << fullname << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		pdb_out_stream << bailout_silent_structure;
		//pdb_out_stream << "REMARK  " <<  moreinfostring << std::endl;
	}

#ifdef BOINC_GRAPHICS
		Boinc::unlock_semaphore();
#else
	pthread_mutex_unlock(&bailout_mutex);
#endif

#ifdef BOINC_GRAPHICS
	}
#endif

	pdb_out_stream.close();
	pdb_out_stream.clear();

	// gzip silent file(s)
	core::io::silent::gzip();

	boinc_end_critical_section();
	protocols::boinc::Boinc::worker_finish_summary( protocols::boinc::Boinc::decoy_count() + 2 , protocols::boinc::Boinc::decoy_count() + 2 , 2  );
	boinc_finish(0); // process terminates right here, nothing else needed.

	return;
}


//////////////////////////////////////////////////////////////////////////////////
// Percentage complete = cpu_time/ cpu_run_time.
// There is an exception, though. If we're getting close to the user's preferred
// run time, always set percentage complete so that it appears that we
// have about ten minutes left...
void
update_pct_complete()
{
	static int pct_complete_blink_counter( 1 );
	if ( pct_complete_blink_counter++ %  PCT_COMPLETE_UPDATE_TIME == 0)
		protocols::boinc::Boinc::update_pct_complete();
}

//////////////////////////////////////////////////////////////////////////////////
void
watchdog_sleep(int const watchdog_time){
#ifdef _WIN32
	Sleep( 1000 * watchdog_time ); // windows -- Sleep function takes milliseconds.
#else
	sleep(watchdog_time); // mac/linux -- Posix thread function takes seconds.
#endif
}

//////////////////////////////////////////////////////////////////////////////////
void*
main_watchdog( void* )
{
	using namespace basic::options;
	using namespace protocols::boinc;

	double current_cpu_time = 0.0;
	std::stringstream moreinfo;

	int watchdog_time = option[ OptionKeys::boinc::watchdog_time ];
	int cpu_run_timeout = option[ OptionKeys::boinc::cpu_run_timeout ];

	int startup_time = 30;
	bool init = false;
	int count_blinks = 0;

	std::cerr << "Watchdog active." << std::endl;
	Boinc boinc_wu = Boinc::instance();

#ifdef BOINC_GRAPHICS
	/* Open the Semaphore */
	// for data sychronization
	Boinc::get_semaphore();
#endif

	// Monitor main and quit when it sends the signal ("watchdog_finish" turns
	// off worker_running flag)
	while (Boinc::is_worker_running()){
		watchdog_sleep(WATCHDOG_BLINK_TIME);

		//A new role for the watchdog.
		// Every time it blinks, update "percentage complete" so that
		// users know that they're making progress.
		update_pct_complete();

		count_blinks++;
		if (count_blinks < watchdog_time) continue;
		//Don't come in too early, allow for about 30 seconds of start up time.
		if (!init && count_blinks < startup_time) continue;

		count_blinks = 0;
		init = true;

		// Rosetta might be suspended or preempted, in which case don't do anything.
		BOINC_STATUS *rosetta_status = new BOINC_STATUS;
		boinc_get_status( rosetta_status );
		if  (!rosetta_status->suspended) {

			// get current working set size (memory use)
			protocols::boinc::Boinc::set_working_set_size( rosetta_status->working_set_size );

			// get current cpu run time
			boinc_wu_cpu_time(current_cpu_time);

			// get max cpu run time user preference
			// user may have updated this
			Boinc::read_and_set_project_prefs();
			int cpu_run_time = boinc_wu.get_project_pref_max_cpu_run_time();

			// Are we taking too long?
			if (current_cpu_time > (cpu_run_timeout + cpu_run_time)  && cpu_run_time > 0) {
				moreinfo << "CPU time: " << current_cpu_time << " seconds. Exceeded timeout" << cpu_run_timeout << " + " <<  cpu_run_time << " seconds";
				std::cerr << "BOINC:: CPU time: " << current_cpu_time << "s, " << cpu_run_timeout << "s + " <<  cpu_run_time << "s";
				get_the_hell_out(moreinfo.str());
				return 0;
			}
		}
	}

	std::cout << "Watchdog finished." << std::endl;
	return 0;
}

} // namespace watchdog
} // namespace boinc
} // namespace protocols

#endif
