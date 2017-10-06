// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/boinc/boinc.hh
/// @brief  Wrappers to make BOINC work
/// @author David Kim (dekim@u.washington.edu)


#ifndef INCLUDED_protocols_boinc_boinc_hh
#define INCLUDED_protocols_boinc_boinc_hh

// External headers: BOINC
#ifdef BOINC
#ifdef WIN32
#include <boinc_win.h> // WIN32 INCLUDE
#include <str_util.h> // WIN32 INCLUDE
#endif
#include <boinc_api.h>
#include <diagnostics.h>
#include <error_numbers.h>
#ifdef BOINC_GRAPHICS
#include <graphics2.h>
#include <util.h>
#include <shmem.h>

// for pose serialization
#define SKIP_FOR_EFFICIENCY 200

#ifndef _WIN32
#ifdef __APPLE__
#include <semaphore.h>
#else
#define USE_SYSV_SEMAPHORE
#define SEMA_KEY_PREFIX 54321
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#endif
#endif


#include <protocols/boinc/boinc_shmem.hh>
#include <core/pose/Pose.fwd.hh>

#endif // BOINC_GRAPHICS
#endif // BOINC


// C++ headers

#define BOINC_MAX_NO_PROGRESS_INIT_CNT 25 // maximum allowed restarts w/ no progress
#define BOINC_MAX_NSTRUCT 99  // maximum nstruct that the assimilator should have to handle. More than that is likely too much to sensibly upload anyway!
#define BOINC_DEFAULT_MAX_CPU_RUN_TIME 10800
#define BOINC_MAX_GFX_FPS 5.0
#define BOINC_MAX_GFX_CPU 10.0
#define BOINC_CHECKPOINT_COUNT_FILE "boinc_checkpoint_count.txt"
#define BOINC_INIT_COUNT_FILE "boinc_init_count.txt"

#ifdef BOINC_GRAPHICS
#define BOINC_SHMEM_NAME  "rosetta"
#endif

namespace protocols {
namespace boinc {

	enum background_type { BLACK_BG=1, BLUE_GRADIENT_BG };

class Boinc {

public:

	/// @brief Access to singleton Boinc class.
	static Boinc& instance(void);

	/// @brief Initialize BOINC. Sets up diagnostics, calls boinc_init(), etc.
	/// Should be used at the start of main in worker app.
	void initialize_worker( void );

	///////////////////////////////////////////////////////////////////////////
	////// User project prefs

	/// @brief Get maximum graphics frames per second user preference
	static double get_project_pref_max_gfx_fps( void );
	/// @brief Get maximum percent cpu for graphics user preference
	static double get_project_pref_max_gfx_cpu( void );
	/// @brief Get maximum target cpu run time user preference
	static int get_project_pref_max_cpu_run_time( void );
	/// @brief Get graphics background type user preference
	static background_type get_project_pref_bg(void);
	/// @brief Set maximum graphics frames per second. Standalone mode only.
	static void set_project_pref_max_gfx_fps( double project_pref_max_gfx_fps );
	/// @brief Set maximum percent cpu for graphics. Standalone mode only.
	static void set_project_pref_max_gfx_cpu( double project_pref_max_gfx_cpu );
	/// @brief Set maximum target cpu run time. Standalone mode only.
	static void set_project_pref_max_cpu_run_time( int project_pref_max_cpu_run_time );
	/// @brief Set working set size to keep track of memory usage
	static void set_working_set_size( double working_set_size );
	/// @brief Rereads project prefs and stores info
	static void read_and_set_project_prefs(void);

	///////////////////////////////////////////////////////////////////////////
	/////// Run status - cpu time and percent complete

	/// @brief Gets the remaining cpu time which depends on the user
	/// run time preference
	static double get_remaining_cpu_time(void);
	/// @brief Determines the fraction complete based on the max_cpu_run_time
	/// preference and the total wu_cpu_run_time. Calls BOINC API
	/// boinc_fraction_done().
	static void update_pct_complete(void);
	/// @brief gets the current runtime for the boinc wu
	static double get_boinc_wu_cpu_time();
	/// @brief Is the task finished?
	static bool worker_is_finished(const int & total_nstruct);

	///////////////////////////////////////////////////////////////////////////
	//////// Startup and shutdown

	/// @brief Start BOINC.
	static void worker_startup( void );
	/// @brief Stop BOINC. Does not return. Uses BOINC API boinc_finish(0).
	static void worker_shutdown( void );
	/// @brief Brief summary. Formatted for the R@h validator so don't change.
	static void worker_finish_summary( const int & num_decoys, const int & attempted_decoys, const int & starting_structs );

	///////////////////////////////////////////////////////////////////////////
	//////// Run progress

	/// @brief Keeps track of number of checkpoints made.
	/// Necessary for determining whether application is making progress on the client
	static int iterate_checkpoint_count( void );
	/// @brief Must be called after a checkpoint is made
	static void checkpoint_completed(void);

	// for watchdog
	static bool is_worker_running( void );
	static void stop_running_worker( void );

	///////////////////////////////////////////////////////////////////////////
	//////// For graphics

#ifdef BOINC_GRAPHICS

	/////////////////////////////////////////////////////////////////////////////
	// WORKER APP GRAPHICS FUNCTIONS

	// called by worker
	static void set_randomly_cycle_appearance( bool const setting );

	static void set_wu_desc( void );

	// The current pose shared memory data will get updated with the observer.
	static int attach_graphics_current_pose_observer( core::pose::Pose & pose );

	// The current pose "ghost" shared memory data will get updated with the observer.
	static int attach_graphics_current_pose_ghost_observer( core::pose::Pose & pose );

	// Sets the native pose
	static int set_graphics_native_pose( core::pose::Pose & pose );

	// Update monte carlo mover trial info
	static void  update_mc_trial_info( const int & trial_cnt, const std::string & mover_type );

	// Update the current Pose buffer directly.
	static void update_graphics_current( core::pose::Pose & pose );

	// Update the current Pose ghost buffer directly.
	static void update_graphics_current_ghost( core::pose::Pose & pose );

	// Low energy pose and energy - set by monte carlo object
	// If PERSIST is set, later calls with DEFAULT will be ignored until RESET is used.
	enum update_mode_enum { DEFAULT, PERSIST, RESET };
	static void update_graphics_low_energy( core::pose::Pose & pose, core::Real low_energy, update_mode_enum update_mode = DEFAULT );

	// Last accepted pose and energy - set by monte carlo object
	static void update_graphics_last_accepted( core::pose::Pose & pose, core::Real last_accepted_energy );

	/////////////////////////////////////////////////////////////////////////////


	// shared memory
	static void create_shared_memory( void );

	/// @brief Signals that the shared memory object has been fully initialized by the worker,
	/// so that the graphics app may read from it.
	/// @details The shared memory object must be fully created and initialized before calling this.
	/// @author Vikram K. Mulligan, Baker lab.
	static void set_shared_memory_fully_initialized();

	/// @brief Waits for the signal that the shared memory object has been fully initialized by the worker,
	/// so that the graphics app may read from it.
	/// @author Vikram K. Mulligan, Baker lab.
	static void wait_for_shared_memory_initialization();

	static void attach_shared_memory( void );
	static void update_status_shmem( void );
	static BoincSharedMemory* get_shmem( void );
#endif // BOINC_GRAPHICS

	// data synchronization
#ifdef USE_SYSV_SEMAPHORE
	// integer key based on the boinc client slot run directory
	static key_t get_sema_key( void );

	static int destroy_semaphore( void );
#endif

	// name based on the boinc client slot run directory
	static void get_sema_name( char * name);

	static int create_semaphore( void );
	static int get_semaphore( void );
	static int wait_semaphore( void );
	static int trywait_semaphore( void );
	static int unlock_semaphore( void );

	static int decoy_count()  { return decoy_count_; };

private:

	Boinc(void) = default;

	static bool worker_initialized_;
	static double fraction_done_;
	static double cpu_time_;
	static double working_set_size_;
	static double working_set_size_max_val_;
	static int checkpoint_count_;
	static int decoy_count_;
	static int no_progress_init_cnt_;
	static bool worker_running_;
	static double project_pref_max_gfx_fps_;
	static double project_pref_max_gfx_cpu_;
	static int project_pref_max_cpu_run_time_;
	static background_type project_pref_bg_;

#ifdef BOINC_GRAPHICS

#ifdef WIN32
	static HANDLE sem_des_;
#else
#ifndef USE_SYSV_SEMAPHORE
	static sem_t* sem_des_;
#else
	static int sem_id_;
#endif
#endif

	static BoincSharedMemory* shmem_;

#endif // BOINC_GRAPHICS

};


} // namespace boinc
} // namespace protocols


#endif // INCLUDED_protocols_boinc_boinc_HH
