// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/MPI_WorkUnitManager.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_MPI_WorkUnitManager_hh
#define INCLUDED_protocols_wum_MPI_WorkUnitManager_hh


#ifdef USEMPI
#include <mpi.h> //keep this first
#else
#define MPI_ANY_SOURCE 0
#endif

#include <protocols/wum/MPI_WorkUnitManager.fwd.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <string>
#include <vector>

#include <utility/vector1.hh>


namespace protocols {
namespace wum {


/// @brief Helper function - returns rank of the current node.
int mpi_rank();

/// @brief Helper function - returns total number of nodes
int mpi_npes();

/// @brief Get a time in seconds. If MPI is enabled you'll get an accuracy of microsecs
core::Real get_time();


class MPI_WorkUnitManager : public WorkUnitManager {
public:
	MPI_WorkUnitManager(
		char machine_letter
	);

	virtual ~MPI_WorkUnitManager(){}

	// @brief Main start function. This gets overloaded of course
	virtual void go() = 0;

protected: // Emperor specific fucntions

	virtual void init() = 0;

	virtual void process_inbound_wus()=0;

	virtual void process_outbound_wus()=0;

protected:

	/// @brief Print a line with general run statistics, now
	virtual void print_stats();

	/// @brief Print a line with general run statistics, only if enough time has pased since the last statistics printout
	virtual void print_stats_auto();

	/// @brief Reset all the stats counters
	void reset_timing_stats();

protected:  /// MPI Communication function

	/// @brief Check for messages and process them accordingly
	void process_incoming_msgs( bool wait_until_message  = false);

	/// @brief Send a workunit to an arbitrary recipient
	void send_MPI_workunit( const WorkUnitBaseOP &wu, int dest_rank  ) const;

	/// @brief Receive a workunit and add it to the inbound queue. By default accept any workunit
	///        or accept a specific source rank. Note that this function is *blocking* and will onyl return once a workunit has been received.
	void receive_MPI_workunit( core::Size node_rank = MPI_ANY_SOURCE );

	void send_next_WU_on_request();

	/// @brief How many masters are there in total ?
	//core::Size n_masters() const;

	/// @brief This is for masters only - returns a serial master number, i.e. 0 is the first master, N-1 is the last master.
	/// Returns the running serial number (0 based) of a master. Returns 0 if emperor
	/// The first master returns 0, the second returns 1, etc..
	//core::Size my_master_rank() const;

	// @brief Return a single char depending on type: E = emperor, M = master and S = Slave
	char get_machine_letter();

	enum MPI_TIMING {
		TIMING_WAIT=0,            // Waiting is defined as waiting for a responce
		TIMING_TRANSFER_SEND,     // Raw sending time, time spent in MPI_send essentially
		TIMING_TRANSFER_RECV,     // Raw receiving time, time spent in MPI_Recv
		TIMING_CPU,               // CPU time
		TIMING_IO_WRITE,          // Disk writing
		TIMING_IO_READ,           // Disk reaading
		TIMING_IDLE,              // Explicit idling (whenever sleep(x) is invoked). This is differnet from WAITing as in it indicated a lack of stuff to do rather then waiting for some other player.
		TIMING_end };             // just an end marker for loops.

	/// @brief This initiates a new timer block. Note that there is no end_timer() function - you just keep calling start_timer, which automatically ends the previous block (and records times etc)
	core::Real start_timer( MPI_TIMING timing_mode ) const;

	/// @brief Display the timing statistics
	void print_timing_stats();

	/// @brief Return the total life time of this class in seconds.
	long wall_time() const;

private:  // ------- variables of the timing system --------------

	/// @brief Minimum time in seconds between statistics updates (to screen)
	core::Size  print_statw_interval_;

	/// @brief when did the last timing block start ? (unixtime)
	mutable core::Real timing_last_start_time_;

	/// @brief Timing type that is currently active
	mutable MPI_TIMING timing_last_type_;

	/// @brief Running sums of total time spent in for any given timing type.
	mutable core::Real timing_total_[TIMING_end];

	/// @brief Start unixtime of this class (set once in contructor)
	long start_time_wall_clock_;

	/// @brief Last unixtime when stats were printed
	mutable long last_stats_;

	/// @brief Total amount of traffic received, in bytes
	mutable core::Size traffic_total_received_;

	/// @brief Total amount of traffic sent, in bytes
	mutable core::Size traffic_total_sent_;

	/// @brief Total time spent in MPI_Send
	mutable core::Real send_wu_time_;

	/// @brief Total number of sends
	mutable core::Real send_wu_time_n_;

	/// @brief Total time spent in MPI_Recv
	mutable core::Real recv_wu_time_;

	/// @brief Total time of receives
	mutable core::Real recv_wu_time_n_;

	char machine_letter_;
};


}
}

#endif

