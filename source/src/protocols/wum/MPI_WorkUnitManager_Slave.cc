// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/wum/MPI_WorkUnitManager_Slave.cc
/// @brief
/// @author Mike Tyka

#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <core/io/silent/SilentStruct.hh>

#include <basic/Tracer.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

using namespace ObjexxFCL::format;

namespace protocols {
namespace wum {

static basic::Tracer TR( "MPI_WUM_Slave" );


MPI_WorkUnitManager_Slave::MPI_WorkUnitManager_Slave( core::Size my_master ):
	MPI_WorkUnitManager( 'S' ),   // this is the one-letter identifier used in printing statistics
	my_master_( my_master )
{
	terminate_ = false;

}


void
MPI_WorkUnitManager_Slave::go()
{
	TR << "Init Slave: " << mpi_rank() << std::endl;
	init();

	do {
		TRDEBUG << "Requesting new jobs.." << std::endl;
		request_new_jobs();

		TRDEBUG << "Slave: Waiting for and processing incoming messages... " << std::endl;

		process_incoming_msgs( true );

		TRDEBUG << "Running jobs" << std::endl;
		process_inbound_wus();

		TRDEBUG << "Slave: Processing outbound queue..." << std::endl;
		process_outbound_wus();

		if ( terminate_ ) break;

		print_stats_auto();
	} while ( true );

	MPI_WorkUnitManager::print_stats();

	TR << "Process terminating!" << std::endl;
}

/// @brief Process workunits on the slave node. I.e. Execute the WU jobs.
void
MPI_WorkUnitManager_Slave::process_inbound_wus()
{
	start_timer( TIMING_CPU );
	// for every workunit on the stack send execute it and transfer it to the outbound queue
	while ( inbound().size() > 0 ) {

		TRDEBUG << "inbound: " << inbound().next()->get_wu_type() << std::endl;

		// Terminate if gets termination WU from Master
		if ( inbound().next()->get_wu_type() == "terminate" ) {
			TRDEBUG << "Received termination sigal from Master." << std::endl;
			//Important to send back to Master that I'm terminating...
			outbound().add( inbound().next() );

			terminate_ = true;
			inbound().pop_next(); // Clean them before you leave!
			continue;
		}

		// run the WU
		runtime_assert( inbound().next() );
		TRDEBUG << "Slave: " << mpi_rank() << " running WU...";
		TRDEBUG << inbound().next()->get_wu_type() << std::endl;

		// special rule for wait workunit. Basically count execution of a wait workunit as
		// IDLING and that of every other one as CPU time
		if ( inbound().next()->get_wu_type() == "waitwu" ) start_timer( TIMING_IDLE );

		inbound().next()->set_run_start();        // record starting time of run for stat analysis
		inbound().next()->run();                  // execute the workunit - i.e. do the work
		inbound().next()->set_run_stop();         // stop the timer just set up above

		start_timer( TIMING_CPU );                // Now we're definitely back in CPU mode

		// Transfer results to the outbound queue
		TRDEBUG << "Slave: " << mpi_rank() << " transferring WU...";
		TRDEBUG << inbound().next()->get_wu_type() << std::endl;

		// Don't add waitwu on outbound!
		// waitwu will attempt to send useless info to master and will prevent having any "real" inbound from master
		//if( inbound().next()->get_wu_type() != "waitwu" ) outbound().add( inbound().next() );
		outbound().add( inbound().next() );
		inbound().pop_next();
	}
}


void
MPI_WorkUnitManager_Slave::process_outbound_wus()
{
	// For each item in the outbound queue, send the WorkUnit back to who ever sent it
	while ( outbound().size() > 0 ) {
		send_MPI_workunit( outbound().next(), (int)outbound().next()->last_received_from() );
		// ERROR behaviour here ?
		outbound().pop_next();  // This clears the WU just sent from the stack
	}
}


void
MPI_WorkUnitManager_Slave::request_new_jobs()
{
#ifdef USEMPI
	unsigned int data;

	// The first blocking send will only return when the master has recevied the data. Since the data transfer is negligible,
	// We'll count that as waiting time, not sending time
	start_timer( TIMING_WAIT );
	TR.Trace << "Requesting new job.." << std::endl;
	core::Size dest_rank = get_my_master();
	runtime_assert( (int)dest_rank != mpi_rank() ); // no self-sending
	MPI_Send( &data, 1, MPI_UNSIGNED, dest_rank,              WUM_MPI_REQUEST_WU, MPI_COMM_WORLD );
	start_timer( TIMING_CPU );
#endif
}


} // namespace wum
} // namespace protocols


