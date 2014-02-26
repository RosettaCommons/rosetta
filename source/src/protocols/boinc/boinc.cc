// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/boinc/boinc.cc
/// @brief  Wrappers to make BOINC work
/// @author David Kim, Mike Tyka
/// @note   Diagnostics are defined at http://boinc.berkeley.edu/trac/wiki/DiagnosticsApi

// Project headers

#ifdef BOINC
// This has to come before boinc.hh or we get this error on VC++
// '_read' : is not a member of 'std::basic_istream<_Elem,_Traits>'
#include <utility/io/izstream.hh>
#endif

// Unit header
#include <protocols/boinc/boinc.hh>

#ifdef BOINC
// Project headers
#include <protocols/boinc/watchdog.hh>
#include <utility/io/ozstream.hh>
#include <utility/basic_sys_util.hh>
#include <utility/string_util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#endif

#include <iostream>


#ifdef BOINC

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

#ifdef BOINC_GRAPHICS
#include <core/pose/Pose.hh>
#include <core/io/serialization/serialize_pose.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/boinc/BoincPoseObserver.hh>
#endif



// C++ headers
#include <iostream>

// option key includes

#include <basic/options/keys/boinc.OptionKeys.gen.hh>



namespace protocols {
namespace boinc {

Boinc& Boinc::instance(void) {
		static Boinc theBoinc;
		return theBoinc;
}

void
Boinc::initialize_worker( void )
{
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::exit;

	if (worker_initialized_) return;

	// DIAGNOSTICS
	// http://boinc.berkeley.edu/trac/wiki/DiagnosticsApi
	boinc_init_diagnostics(
		BOINC_DIAG_DUMPCALLSTACKENABLED |
		BOINC_DIAG_HEAPCHECKENABLED |
		BOINC_DIAG_MEMORYLEAKCHECKENABLED |
		BOINC_DIAG_REDIRECTSTDERR |
		BOINC_DIAG_REDIRECTSTDOUT |
		BOINC_DIAG_TRACETOSTDERR
	);

	cerr << utility::timestamp() << " :: BOINC:: Initializing ... ok." << std::endl;
	if ( boinc_init() ) exit( EXIT_FAILURE ); // BOINC API call
	cerr << utility::timestamp() << " :: BOINC :: boinc_init()" << endl;

	// allow upload of stdout.txt
	char filename[256];
	boinc_resolve_filename("stdout.txt", filename, 256);
	freopen(filename, "a", stdout);

	boinc_fraction_done(0.00);

#ifdef BOINC_GRAPHICS
	// create shared mem segment for graphics, and arrange to update it
	//
	// http://boinc.berkeley.edu/trac/wiki/GraphicsApi
	// create shared memory segment
	cerr << "BOINC:: Setting up shared resources ... ok." << std::endl;

	create_shared_memory();

	cerr << "BOINC:: Setting up semaphores ... ok." << std::endl;
	// create samaphore for data synchronization
	create_semaphore();

	cerr << "BOINC:: Updating status ... ok." << std::endl;

	// updated status info
	update_status_shmem();

	cerr << "BOINC:: Registering timer callback... ok." << std::endl;

	boinc_register_timer_callback(update_status_shmem);
#endif

	worker_initialized_ = true;

	cerr << "BOINC:: Worker initialized successfully. " << std::endl;
}


double Boinc::get_project_pref_max_gfx_fps() { return project_pref_max_gfx_fps_; }
double Boinc::get_project_pref_max_gfx_cpu() { return project_pref_max_gfx_cpu_; }
int Boinc::get_project_pref_max_cpu_run_time() { return project_pref_max_cpu_run_time_; }

void Boinc::set_project_pref_max_gfx_fps( double project_pref_max_gfx_fps ) {
	project_pref_max_gfx_fps_ = project_pref_max_gfx_fps;
}
void Boinc::set_project_pref_max_gfx_cpu( double project_pref_max_gfx_cpu ) {
	project_pref_max_gfx_cpu_ = project_pref_max_gfx_cpu;
}
void Boinc::set_project_pref_max_cpu_run_time( int project_pref_max_cpu_run_time ) {
	project_pref_max_cpu_run_time_ = project_pref_max_cpu_run_time;
}

void Boinc::set_working_set_size( double working_set_size ) {
	working_set_size_ = working_set_size;
	if (working_set_size_max_val_ < working_set_size_) {
		working_set_size_max_val_ = working_set_size_;
	}
}

void Boinc::read_and_set_project_prefs() {

	// parse the boinc init data file. may not need to do this.
	//boinc_parse_init_data_file(); // BOINC API call

	/*
Get information from the core client; this information may be useful for graphics.
See http://boinc.berkeley.edu/trac/wiki/BasicApi
	*/

	APP_INIT_DATA app_init_data;
	boinc_get_init_data(app_init_data); // BOINC API call
	if (!app_init_data.project_preferences) return;

	double max_gfx_fpstmp = 0.0;
	double max_gfx_cputmp = 0.0;
	int cpu_run_timetmp = 0;
	parse_double(app_init_data.project_preferences, "<max_fps>", max_gfx_fpstmp);
	if (max_gfx_fpstmp > 0) project_pref_max_gfx_fps_ = max_gfx_fpstmp;
	parse_double(app_init_data.project_preferences, "<max_cpu>", max_gfx_cputmp);
	if (max_gfx_cputmp > 0) project_pref_max_gfx_cpu_ = max_gfx_cputmp;
	parse_int(app_init_data.project_preferences, "<cpu_run_time>", cpu_run_timetmp);
	if (cpu_run_timetmp > 0 && cpu_run_timetmp != project_pref_max_cpu_run_time_) {
		project_pref_max_cpu_run_time_ = cpu_run_timetmp;
		std::cerr << "# cpu_run_time_pref: " << project_pref_max_cpu_run_time_ << std::endl;std::cerr.flush();
	}
}

double Boinc::get_remaining_cpu_time() {
	read_and_set_project_prefs(); // get user max cpu run time preference
	boinc_wu_cpu_time(cpu_time_); // BOINC API call for wu cpu run time
	return project_pref_max_cpu_run_time_ - cpu_time_;
}

void Boinc::update_pct_complete() {
	read_and_set_project_prefs(); // get user max cpu run time preference
	boinc_wu_cpu_time(cpu_time_); // BOINC API call for wu cpu run time
	// Estimate the percentage done depending on the number of decoys already produced:
	// If we're still working on the first decoy, assume that it may well take more then the user's preferred time, in
	// fact it may take up to the watchdog timeout in theory (3x runtime!)
	// If we've done many decoys then safely make a more accurate estimate of the time remaining.
	using namespace basic::options;
	int cpu_run_timeout = option[ OptionKeys::boinc::cpu_run_timeout ]();
	double estimated_wu_length = (double)project_pref_max_cpu_run_time_  + (double(cpu_run_timeout) / (double(decoy_count_) + 1.0) );

	if( ( estimated_wu_length - cpu_time_ ) < 10*60 ) estimated_wu_length = cpu_time_ + 10*60;
	double new_frac = cpu_time_ / ( estimated_wu_length );

	//Only update if increasing
	if (new_frac > fraction_done_) {
		fraction_done_ = new_frac;
		boinc_fraction_done(fraction_done_);
/*
		std::cout << utility::timestamp() << " :: BOINC" <<
					" :: cpu_time_pref: " << project_pref_max_cpu_run_time_ <<
					" :: cpu_time: " << cpu_time_ <<
					" :: fraction_done: " << ObjexxFCL::format::F( 7, 2, fraction_done_ ) << std::endl;
*/
	}
}

//call to BOINC API call for wu cpu run time
double Boinc::get_boinc_wu_cpu_time(){
		double tmp_cpu_time;
		boinc_wu_cpu_time(tmp_cpu_time);
		return(tmp_cpu_time);
}


bool Boinc::worker_is_finished( const int & total_nstruct ){
	// finished if too many restarts without progress
	boinc_wu_cpu_time(cpu_time_);
	if ( ( !boinc_is_standalone() ) &&                                     // don't do this in standalone mode
			 ( no_progress_init_cnt_ >= BOINC_MAX_NO_PROGRESS_INIT_CNT ) &&    // nor if the limit of restarts is reached
			 ( cpu_time_	> (20*60) )                                          // nor in the first 20 minutes (costhe validator will kill it!)
		) {
		std::cerr << "Too many restarts with no progress. Keep application in memory while preempted." << std::endl;
		return true;
	}

	// did we just complete a decoy ?
	if( decoy_count_ < total_nstruct ){
		// then reset the no_progress_init_cnt_
		// counter.
		no_progress_init_cnt_ = 0;
	}
	decoy_count_ = std::max(0, total_nstruct);
#ifdef BOINC_GRAPHICS
	if (shmem_) {
		shmem_->model_count = total_nstruct + 1;
		shmem_->total_mc_trial_count = 0; // reset total monte carlo trial count
		if (shmem_->low_energy) {
			shmem_->model_low_energy = shmem_->low_energy;
			if (shmem_->native_pose_exists) {
				// get native
				core::pose::PoseOP nativeposeOP = new core::pose::Pose;
				core::io::serialization::BUFFER bn((char*)(&shmem_->native_pose_buf ),protocols::boinc::POSE_BUFSIZE);
				core::io::serialization::read_binary(*nativeposeOP,bn);
				if (shmem_->low_energy_pose_exists) {
					core::pose::PoseOP lowenergyposeOP = new core::pose::Pose;
					core::io::serialization::BUFFER bl((char*)(&shmem_->low_energy_pose_buf ),protocols::boinc::POSE_BUFSIZE);
					core::io::serialization::read_binary(*lowenergyposeOP,bl);
					shmem_->model_low_energy_rmsd = core::scoring::native_CA_rmsd( *nativeposeOP, *lowenergyposeOP);
				}
			}
		}
	}
#endif
	// not finished if a structure is not made yet
	if (total_nstruct <= 0) return false;
	double time_left = get_remaining_cpu_time();
	using namespace basic::options;
	if (
			// finished if max_cpu_run_time_ reached
			time_left < 1 ||
			// finished if there is not enough time to make another structure
			time_left < cpu_time_/double(total_nstruct) ||
			// finished if too many models
			total_nstruct >= option[ OptionKeys::boinc::max_nstruct ]
		) {
		// FINISHED!
		// update status for graphics
		fraction_done_ = 1.0;
		boinc_fraction_done(fraction_done_);
#ifdef BOINC_GRAPHICS
		update_status_shmem();
#endif
		return true;
	}
	return false;
}

void Boinc::worker_startup() {
	std::cerr << "BOINC:: Worker startup. " << std::endl;
	using namespace basic::options;
	if (!worker_initialized_)
		utility_exit_with_message( "Must initialize Boinc before calling startup." );

	if ( option[ OptionKeys::boinc::frame_rate ].user() )
		project_pref_max_gfx_fps_ = option[ OptionKeys::boinc::frame_rate ];
	if ( option[ OptionKeys::boinc::cpu_frac].user() )
		project_pref_max_gfx_cpu_ = option[ OptionKeys::boinc::cpu_frac ];
	if ( option[ OptionKeys::boinc::cpu_run_time].user() )
		project_pref_max_cpu_run_time_ = option[ OptionKeys::boinc::cpu_run_time ];

/**
The following code is meant to keep track of progress upon restarts.
Clients that do not keep tasks in memory may stop jobs before they
checkpoint which can cause an endless cycle of restarts wasting
cpu time. If there are too many restarts without a checkpoint, the
job will end (see Boinc::is_finished()).
**/
	// restore checkpoint count
	utility::io::izstream chckptcnt_istream( BOINC_CHECKPOINT_COUNT_FILE );
	if ( chckptcnt_istream ) {
		chckptcnt_istream >> checkpoint_count_;
		chckptcnt_istream.close();
		chckptcnt_istream.clear();
	}

	// check number of restarts
	// read last init count
	int prev_initcnt = 0;
	int prev_chckptcnt = 0;
	utility::io::izstream initcnt_istream( BOINC_INIT_COUNT_FILE );
	if ( initcnt_istream ) {
		initcnt_istream >> prev_initcnt >> prev_chckptcnt;
		initcnt_istream.close();
		initcnt_istream.clear();
	}
	// iterate no_progress_init_cnt if checkpoint count has not increased
	no_progress_init_cnt_ = ( checkpoint_count_ > prev_chckptcnt ) ? prev_initcnt : prev_initcnt + 1;

	// update counts
	utility::io::ozstream initcnt_ostream( BOINC_INIT_COUNT_FILE );
	if ( initcnt_ostream ) {
		initcnt_ostream << no_progress_init_cnt_ << " " << checkpoint_count_;
		initcnt_ostream.close();
		initcnt_ostream.clear();
	}

	worker_running_ = true;
	watchdog::watchdog_start();

}

void Boinc::worker_shutdown() {

	// already did this in worker_is_finished, but oh well
	fraction_done_ = 1.0;
	boinc_fraction_done(fraction_done_);
#ifdef BOINC_GRAPHICS
	update_status_shmem();
#endif

	worker_running_ = false;

	// report memory usage
	// do not edit this because it is parsed from stderr by the
	// RALPH@home sub batch page to report memory usage
	std::cerr << "BOINC :: WS_max " << working_set_size_max_val_ << std::endl;

	watchdog::watchdog_finish();
	std::cerr << "BOINC :: BOINC support services shutting down cleanly ..." << std::endl;
	boinc_finish(0); // BOINC API call. Does not return.
}

// the validator requires this format printed to stderr upon job completion
void Boinc::worker_finish_summary( const int & num_decoys, const int & attempted_decoys, const int & starting_structs ) {
	using namespace ObjexxFCL::format;
	// job distributor iterates decoy count before ending so decrement it for the
	// actual count
	int decoy_cnt = (num_decoys > 0) ? num_decoys-1 : num_decoys;
	int attempt_cnt = (attempted_decoys > 0) ? attempted_decoys-1 : attempted_decoys;
	double cputime = 0.0;
	boinc_wu_cpu_time(cputime); // get wu cpu run time
	// prevent invalid cpu times (validator doesnt accept things that took less then 1200 seconds ):
	if(cputime < 1200 ) cputime = 1201;

	std::cerr << "======================================================" << std::endl;
	std::cerr << "DONE ::" << I( 6, std::min( decoy_cnt, std::max(starting_structs,1) ) ) <<
		" starting structures " << I( 8, cputime ) << " cpu seconds" << std::endl;
	std::cerr << "This process generated " <<
		I( 6, decoy_cnt ) << " decoys from " << I( 7, attempt_cnt ) <<
		" attempts" << std::endl;
	std::cerr << "======================================================" << std::endl;
}

int Boinc::iterate_checkpoint_count() {
	checkpoint_count_++;
	return checkpoint_count_;
}

// must be called after a checkpoint to keep track of progress
void Boinc::checkpoint_completed() {

	int checkpointcount = iterate_checkpoint_count();
	// save checkpoint count to keep track of progress
	utility::io::ozstream chckptcnt_ostream( BOINC_CHECKPOINT_COUNT_FILE );
	if ( chckptcnt_ostream ) {
		chckptcnt_ostream << checkpointcount;
		chckptcnt_ostream.close();
		chckptcnt_ostream.clear();
	}
	boinc_checkpoint_completed();

}


// for watchdog
bool Boinc::is_worker_running() {
	return worker_running_;
}

void Boinc::stop_running_worker() {
	worker_running_ = false;
}


#ifdef BOINC_GRAPHICS

//////////////////////////////////////////////////////////////////////////////////////
// FOR GRAPHICS SHARED MEMORY
//////////////////////////////////////////////////////////////////////////////////////

// called by worker
void Boinc::set_wu_desc() {
  using namespace basic::options;
	if (!shmem_) return;
  static bool init = false;
  if (!init){
		// read description file if one exists and keep row format
		std::string description_file;
		if ( option[ OptionKeys::boinc::description_file ].user() ) {
			std::string description_file = option[ OptionKeys::boinc::description_file ];
			utility::io::izstream desc_stream( description_file );
			if ( desc_stream ) {
				std::vector<std::string> wu_desc_rows;
				std::string tmpstr;
				while (desc_stream && !desc_stream.eof()) {
					desc_stream.getline( tmpstr );
					utility::trim( tmpstr );
					if ( tmpstr.length() > 0)
						wu_desc_rows.push_back( tmpstr );
				}
				if (!wait_semaphore()) {
					core::io::serialization::BUFFER b((char*)(&shmem_->wu_desc_buf),POSE_BUFSIZE);
					core::io::serialization::write_binary(wu_desc_rows,b);
					shmem_->wu_desc_exists = 1;
					unlock_semaphore();
				}
			}
			desc_stream.close();
			desc_stream.clear();
		}
		init = true;
  }
}

// CURRENT POSE
// called by worker
const int Boinc::attach_graphics_current_pose_observer( core::pose::Pose & pose ) {
	if (!shmem_) return 0;
	static BoincCurrentPoseObserverOP bpo = new BoincCurrentPoseObserver();
	bpo->attach_to( pose );
	return 1;
}

// NATIVE POSE
// called by worker
const int Boinc::set_graphics_native_pose( core::pose::Pose & pose ) {
	if (!shmem_) return 0;
	if (!wait_semaphore()) {
		boinc_begin_critical_section();
		core::io::serialization::BUFFER b((char*)(&shmem_->native_pose_buf),POSE_BUFSIZE);
		core::io::serialization::write_binary(pose,b);
		shmem_->native_pose_exists = 1;
		boinc_end_critical_section();
		unlock_semaphore();
	}
	return 1;
}

void  Boinc::update_mc_trial_info( const int & trial_cnt, const std::string & mover_type ) {
	if (shmem_) {
		// stage info
		core::io::serialization::BUFFER b((char*)(&shmem_->mover_type_text),TEXT_BUFSIZE);
		core::io::serialization::write_binary(mover_type,b);
		// monte carlo mover step count
		shmem_->mover_mc_trial_count = trial_cnt;
		shmem_->total_mc_trial_count++;
	}
}


// MC LOW ENERGY POSE
// called by worker mc object
void Boinc::update_graphics_current( core::pose::Pose & pose ) {
/*
	static int count = 0;
  if (count > 0) {
		count--;
		return;
  }
  count = SKIP_FOR_EFFICIENCY;
*/
	if (!shmem_) return;
	if (!trywait_semaphore()) {
		boinc_begin_critical_section();
		if (pose.total_residue() > 0) {
			core::io::serialization::BUFFER b((char*)(&shmem_->current_pose_buf ),POSE_BUFSIZE);
			core::io::serialization::write_binary(pose,b);
			shmem_->current_pose_exists = 1;
		}
		boinc_end_critical_section();
		unlock_semaphore();
  }
}


// MC LOW ENERGY POSE
// called by worker mc object
void Boinc::update_graphics_low_energy( core::pose::Pose & pose, core::Real low_energy, update_mode_enum update_mode ) {
	static bool persist = false;
	if (update_mode == RESET) {
		persist = false;
	} else if (update_mode == PERSIST) {
		persist = true;
	}
/*
	static int count = 0;
  if (count > 0) {
		count--;
		return;
  }
  count = SKIP_FOR_EFFICIENCY;
*/
	if (!shmem_ || (persist && update_mode == DEFAULT)) return;
  if (!trywait_semaphore()) {
		boinc_begin_critical_section();
		shmem_->low_energy = low_energy;
		core::io::serialization::BUFFER b((char*)(&shmem_->low_energy_pose_buf ),POSE_BUFSIZE);
		core::io::serialization::write_binary(pose,b);
		shmem_->low_energy_update_cnt++;
		shmem_->low_energy_pose_exists = 1;
		boinc_end_critical_section();
		unlock_semaphore();
  }
}

// MC LAST ACCEPTED POSE
// called by worker mc object
void Boinc::update_graphics_last_accepted( core::pose::Pose & pose, core::Real last_accepted_energy ) {
/*
	static int count = 0;
  if (count > 0) {
		count--;
		return;
  }
  count = SKIP_FOR_EFFICIENCY;
*/
	if (!shmem_) return;
  if (!trywait_semaphore()) {
		boinc_begin_critical_section();
		shmem_->last_accepted_energy = last_accepted_energy;
		core::io::serialization::BUFFER b((char*)(&shmem_->last_accepted_pose_buf ),POSE_BUFSIZE);
		core::io::serialization::write_binary(pose,b);
		shmem_->last_accepted_pose_exists = 1;
		boinc_end_critical_section();
		unlock_semaphore();
  }
}


void Boinc::create_shared_memory() {
	if (shmem_ == NULL) {
		shmem_ = (BoincSharedMemory*)boinc_graphics_make_shmem(BOINC_SHMEM_NAME, sizeof(BoincSharedMemory));
		if (!shmem_) {
			std::cerr << "failed to create shared mem segment: " << BOINC_SHMEM_NAME << " Size: " << sizeof(BoincSharedMemory) << std::endl;
		} else {
			std::cout << "Created shared memory segment " << std::endl;
		}
	}
}

void Boinc::attach_shared_memory() {
	if (shmem_ == NULL) {
		int attempts = 15;
		while (attempts>1) {
			shmem_ = (protocols::boinc::BoincSharedMemory*)boinc_graphics_get_shmem(BOINC_SHMEM_NAME);
			if (!shmem_) {
#ifdef _WIN32
				Sleep( 1000 );
#else
				sleep(1);
#endif
				attempts--;
			} else {
				std::cout << "Attached shared memory segment " << std::endl;
				return;
			}
		}
		std::cerr << "failed to attach shared mem segment " << std::endl;
	}
}

BoincSharedMemory* Boinc::get_shmem() { return shmem_; }

void Boinc::update_status_shmem() {
	if (!shmem_) return;
	shmem_->fraction_done = fraction_done_; //boinc_get_fraction_done();
	boinc_wu_cpu_time(cpu_time_);
	shmem_->cpu_time = cpu_time_; //boinc_worker_thread_cpu_time();
	shmem_->update_time = dtime();
	boinc_get_status(&shmem_->status);
}


// SEMAPHORE STUFF
// This stuff should be moved into its own class eventually.
//
// WINDOWS IMPLEMENTATION MAY HAVE SECURITY ISSUES WHEN RUN AS A SERVICE
// I don't know how to set the attributes and what to set it to in CreateSemaphore
//
// For graphics shared memory data synchronization

#ifndef USE_SYSV_SEMAPHORE

void Boinc::get_sema_name( char * name ) {
	// make it slot specific since multiple tasks may run simultaneously
	APP_INIT_DATA aid;
	int retval = boinc_get_init_data(aid);
	if (retval) aid.slot = 0;
	sprintf(name, "boinc_minirosetta_sema_%d", aid.slot);
}

#else

const key_t Boinc::get_sema_key() {
	// make it slot specific since multiple tasks may run simultaneously
	static int sema_key;
  if (sema_key) return sema_key;
  APP_INIT_DATA aid;
  int retval = boinc_get_init_data(aid);
  if (retval) aid.slot = 0;
  sema_key = SEMA_KEY_PREFIX + aid.slot;
	return sema_key;
}

const int Boinc::destroy_semaphore() {
	int id = semget(get_sema_key(), 0, 0);
	if (id < 0) return 0;
	if (semctl(id, 1, IPC_RMID, 0)) return 0;
	sem_id_ = 0;
	return 1;
}

#endif // USE_SYSV_SEMAPHORE

int Boinc::create_semaphore() {

#ifdef _WIN32
  // Windows semaphore

	char name[256];
	get_sema_name(name);

	CloseHandle(sem_des_);
	// to do: use a specific security attribute instead of null (default attribute)?

// WINDOWS SECURITY DESCRIPTOR FROM BOINC

	// Win9X doesn't like any reference to a security descriptor.
	// So if we detect that we are running on the Win9X platform pass
	// a NULL value for it.
	//

	DWORD dwError = 0;
	DWORD dwRes = 0;
	PSID pEveryoneSID = NULL;
	PACL pACL = NULL;
	PSECURITY_DESCRIPTOR pSD = NULL;
	EXPLICIT_ACCESS ea;
	SID_IDENTIFIER_AUTHORITY SIDAuthWorld = SECURITY_WORLD_SID_AUTHORITY;
	SECURITY_ATTRIBUTES sa;
	OSVERSIONINFO osvi;
	char global_sema_name[256];

	osvi.dwOSVersionInfoSize = sizeof(osvi);
	GetVersionEx(&osvi);
	if (osvi.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS) {
		sem_des_ = CreateSemaphore(NULL,1,1, name);
	} else {
		// Create a well-known SID for the Everyone group.
		if(!AllocateAndInitializeSid(&SIDAuthWorld, 1,
			SECURITY_WORLD_RID,
			0, 0, 0, 0, 0, 0, 0,
			&pEveryoneSID)
		) {
			std::cerr << "CreateSemaphore failure: AllocateAndInitializeSid Error " << GetLastError() << std::endl;
			goto Cleanup;
		}

		// Initialize an EXPLICIT_ACCESS structure for an ACE.
		// The ACE will allow Everyone all access to the shared memory object.
		ZeroMemory(&ea, sizeof(EXPLICIT_ACCESS));
		ea.grfAccessPermissions = SEMAPHORE_ALL_ACCESS;
		ea.grfAccessMode = SET_ACCESS;
		ea.grfInheritance= NO_INHERITANCE;
		ea.Trustee.TrusteeForm = TRUSTEE_IS_SID;
		ea.Trustee.TrusteeType = TRUSTEE_IS_WELL_KNOWN_GROUP;
		ea.Trustee.ptstrName  = (LPTSTR) pEveryoneSID;

		// Create a new ACL that contains the new ACEs.
		dwRes = SetEntriesInAcl(1, &ea, NULL, &pACL);
		if (ERROR_SUCCESS != dwRes) {
				std::cerr <<  "CreateSemaphore failure: SetEntriesInAcl " << GetLastError() << std::endl;
				goto Cleanup;
		}

		// Initialize a security descriptor.
		pSD = (PSECURITY_DESCRIPTOR) LocalAlloc(LPTR, SECURITY_DESCRIPTOR_MIN_LENGTH);
		if (NULL == pSD) {
			std::cerr << "CreateSemaphore failure: LocalAlloc Error " << GetLastError() << std::endl;
			goto Cleanup;
		}

		if (!InitializeSecurityDescriptor(pSD, SECURITY_DESCRIPTOR_REVISION)) {
			std::cerr << "CreateSemaphore failure: InitializeSecurityDescriptor  Error " << GetLastError() << std::endl;
			goto Cleanup;
		}

		// Add the ACL to the security descriptor.
		if (!SetSecurityDescriptorDacl(pSD,
				TRUE,     // bDaclPresent flag
				pACL,
		    FALSE) // not a default DACL
		) {
			std::cerr << "CreateSemaphore failure: SetSecurityDescriptorDacl Error " << GetLastError() <<  std::endl;
			goto Cleanup;
		}

		// Initialize a security attributes structure.
		sa.nLength = sizeof (SECURITY_ATTRIBUTES);
		sa.lpSecurityDescriptor = pSD;
		sa.bInheritHandle = FALSE;

		// Use the security attributes to set the security descriptor
		// when you create a shared file mapping.

		// Try using 'Global' so that it can cross terminal server sessions
		// The 'Global' prefix must be included in the shared memory
		// name if the shared memory segment is going to cross
		// terminal server session boundaries.
		//
		bool try_global = true;
		if (try_global) {
			sprintf(global_sema_name, "Global\\%s", name);
			sem_des_ = CreateSemaphore(&sa, 1, 1, global_sema_name);
			dwError = GetLastError();
			if (!sem_des_ && (ERROR_ACCESS_DENIED == dwError)) {
				// Couldn't use the 'Global' tag, so try the original name.
				try_global = false;
			}
		}
		if (!try_global) {
			sem_des_ = CreateSemaphore(&sa, 1, 1, name);
			dwError = GetLastError();
		}
	}

	if (sem_des_) {
		if (GetLastError() == ERROR_ALREADY_EXISTS) {
			CloseHandle(sem_des_);
			sem_des_ = NULL;
		}
	}

	if (!sem_des_) {
		std::cerr << "CreateSemaphore failure! Cannot create semaphore!" << std::endl;
		return 0;
	} else {
		std::cout << "Created semaphore" << std::endl;
		return 1;
	}

Cleanup:

	if (osvi.dwPlatformId != VER_PLATFORM_WIN32_WINDOWS) {
		if (pEveryoneSID)
			FreeSid(pEveryoneSID);
		if (pACL)
			LocalFree(pACL);
		if (pSD)
			LocalFree(pSD);
	}
	return 0;

#else

#ifndef USE_SYSV_SEMAPHORE
// POSIX semaphore

	char name[256];
	get_sema_name(name);

	sem_close(sem_des_);
	sem_unlink(name);
	sem_des_ = sem_open(name, O_CREAT|O_EXCL, 0777, 1);
	if(sem_des_ == (void*)-1) {
		std::cerr << "sem_open create failure! Cannot create semaphore!" << std::endl;
		return 0;
	} else {
		std::cout << "Created semaphore" << std::endl;
		return 1;
	}

#else
// System V semaphore
// POSIX semaphore gave a seg fault when sem_wait was called on linux

  destroy_semaphore();

	union semun
	{
		int val;
		struct semid_ds *buf;
		unsigned short *array;
	} arg;
	arg.val = 1;

	int semid = semget(get_sema_key(), 1, IPC_CREAT|IPC_EXCL|0777);
	if (semid < 0) {
    std::cerr << "semget failure! Cannot create semaphore!" << std::endl;
    return 0;
  }
	if( semctl(semid, 0, SETVAL, arg) < 0) {
    std::cerr << "semctl failure! Cannot set semaphore value!" << std::endl;
		return 0;
	}
	sem_id_ = semid;
	std::cout << "Created semaphore" << std::endl;
	return 1;

#endif // USE_SYSV_SEMAPHORE
#endif // _WIN32

}

int Boinc::get_semaphore() {

#ifdef _WIN32
	// already opened
	if (sem_des_) return 1;

	char name[256];
	get_sema_name(name);

	char global_sema_name[256];

	// The 'Global' prefix must be included in the shared memory
	// name if the shared memory segment is going to cross
	// terminal server session boundries.
	//
	sprintf(global_sema_name, "Global\\%s", name);
	sem_des_ = OpenSemaphore(SEMAPHORE_ALL_ACCESS, NULL, global_sema_name);
	if (!sem_des_) {
		// Couldn't use the 'Global' tag, so just attempt to use the original
		// name.
		//
		sem_des_ = OpenSemaphore(SEMAPHORE_ALL_ACCESS, NULL, name);
	}
	if (!sem_des_) {
		std::cerr << "OpenSemaphore failure" << std::endl;
		return 0;
	}
	std::cerr << "Opened semaphore" << std::endl;
	return 1;
#else
#ifndef USE_SYSV_SEMAPHORE

	// already opened
	if (sem_des_) return 1;

	char name[256];
	get_sema_name(name);
	sem_des_ = sem_open(name, 0, 0777, 0);
	if(sem_des_ == (void*) -1){
		std::cerr << "sem_open failure" << std::endl;
		return 0;
	}
	std::cout << "Opened semaphore" << std::endl;
	return 1;
#else
	if (sem_id_) return 1;

	int semid = semget(get_sema_key(), 0, 0);
	if (semid < 0) {
		std::cerr << "semget failure! Cannot open semaphore!" << std::endl;
		return 1;
	}
	sem_id_ = semid;
	std::cout << "Opened semaphore" << std::endl;
	return 0;
#endif // USE_SYSV_SEMAPHORE
#endif // _WIN32
}

int Boinc::wait_semaphore() {
#ifdef _WIN32
	if (WaitForSingleObject(sem_des_, INFINITE) == WAIT_OBJECT_0) return 0;
	return 1;
#else
#ifndef USE_SYSV_SEMAPHORE
	return sem_wait(sem_des_);
#else
	if (sem_id_ <= 0) return 1;
	struct sembuf ops[] = {{0, -1, 0}};
	if (semop(sem_id_, ops, 1)) return 1;
	return 0;
#endif // USE_SYSV_SEMAPHORE
#endif // _WIN32
}

int Boinc::trywait_semaphore() {
#ifdef _WIN32
	if (WaitForSingleObject(sem_des_, 0) == WAIT_OBJECT_0) return 0;
	return 1;
#else
#ifndef USE_SYSV_SEMAPHORE
	return sem_trywait(sem_des_);
#else
	if (sem_id_ <= 0) return 1;
	struct sembuf ops[] = {{0, -1, IPC_NOWAIT}};
	if (semop(sem_id_, ops, 1)) return 1;
	return 0;
#endif // USE_SYSV_SEMAPHORE
#endif // _WIN32
}

int Boinc::unlock_semaphore() {
#ifdef _WIN32
	// returns non-zero on success
	if (ReleaseSemaphore(sem_des_, 1, NULL)) return 0;
	return 1;
#else
#ifndef USE_SYSV_SEMAPHORE
	return sem_post(sem_des_);  // returns zero on success
#else
	if (sem_id_ <= 0) return 1;
	struct sembuf ops[] = {{0, 1, 0}};
	if (semop(sem_id_, ops, 1)) return 1;
	return 0;
#endif // USE_SYSV_SEMAPHORE
#endif // _WIN32
}

#ifdef _WIN32
HANDLE Boinc::sem_des_ = NULL;
#else
#ifndef USE_SYSV_SEMAPHORE
sem_t* Boinc::sem_des_ = NULL;
#else
int Boinc::sem_id_;
#endif // USE_SYSV_SEMAPHORE
#endif // _WIN32

// SEMAPHORE STUFF */

BoincSharedMemory* Boinc::shmem_ = NULL;

#endif  // BOINC_GRAPHICS

// DATA

// worker status
bool Boinc::worker_initialized_ = false;
double Boinc::fraction_done_ = 0.0;
double Boinc::cpu_time_ = 0.0;
double Boinc::working_set_size_ = 0.0;
double Boinc::working_set_size_max_val_ = 0.0;

// for checking progress
int Boinc::checkpoint_count_ = 0;
int Boinc::decoy_count_ = 0;
int Boinc::no_progress_init_cnt_ = 0;

// for watchdog
bool Boinc::worker_running_ = false;

// default project prefs
// don't change these defaults or boinc users may get upset
// these are user changable preferences
double Boinc::project_pref_max_gfx_fps_ = BOINC_MAX_GFX_FPS;
double Boinc::project_pref_max_gfx_cpu_ = BOINC_MAX_GFX_CPU;
int Boinc::project_pref_max_cpu_run_time_ = BOINC_DEFAULT_MAX_CPU_RUN_TIME;

} // namespace boinc
} // namespace protocols

#endif
