// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopHashMap.cc
/// @brief
/// @author Mike Tyka

#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <utility/assert.hh> //MPI_ONLY macro

#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>


// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
/// ObjexxFCL headers
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>



#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

using namespace ObjexxFCL::format;

namespace protocols {
namespace wum {

static thread_local basic::Tracer TR( "MPI_WUM" );

int mpi_rank(){
	int mpi_rank_=0;
	#ifdef USEMPI
		MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	#else
		utility_exit_with_message( "ERROR: The MPI_WorkUnitManager will not work unless you have compiled using extras=mpi" );
	#endif
	return mpi_rank_;
}

int mpi_npes(){
	int mpi_npes_=0;
	#ifdef USEMPI
		MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );
	#else
		utility_exit_with_message( "ERROR: The MPI_WorkUnitManager will not work unless you have compiled using extras=mpi" );
	#endif
	return mpi_npes_;
}



core::Real get_time(){
	#ifdef USEMPI
	 return MPI_Wtime();
	#else
	 return (core::Real) time(NULL);
	#endif
}








MPI_WorkUnitManager::MPI_WorkUnitManager( char machine_letter  ):
WorkUnitManager(),
last_stats_(0),
traffic_total_received_(0),
traffic_total_sent_(0),
send_wu_time_(0),
send_wu_time_n_(0),
recv_wu_time_(0),
recv_wu_time_n_(0),
machine_letter_( machine_letter )
{
	TR << "Starting MPI_WorkUnitManager.." << std::endl;
 	using namespace basic::options;
  using namespace basic::options::OptionKeys;

	start_time_wall_clock_ = time(NULL);
	timing_last_start_time_ = 0;
	timing_last_type_ = TIMING_CPU;
	reset_timing_stats();
	start_timer( TIMING_CPU );

	TR << "This is node " << mpi_rank() << " Nprocs: " << mpi_npes() << std::endl;

	outbound().set_memory_limit(  option[ OptionKeys::wum::memory_limit ]() * 1000 );  // memory_limit option is in Kilobytes!
	inbound().set_memory_limit(  option[ OptionKeys::wum::memory_limit ]() * 1000 );
}




char MPI_WorkUnitManager::get_machine_letter(){
	return machine_letter_;
}





void
MPI_WorkUnitManager::process_incoming_msgs( bool MPI_ONLY( wait_until_message ) )
{
#ifdef USEMPI
	while(true){
		// Check if there's anything on the line
		start_timer( TIMING_IDLE );
		MPI_Status status;
		int result;
		TR.Trace << "Probing for incoming messages .." << std::endl;

		core::Size before_size =  outbound().size();
		while( outbound().size() > 0 ){
			TR.Trace << "Fulfilling work requests since we have outbound work" << std::endl;
			MPI_Iprobe( MPI_ANY_SOURCE, WUM_MPI_REQUEST_WU, MPI_COMM_WORLD, &result, &status);
			if( !result){ // If there are no work requests
				break;      // break out and continue to accept *any* messages
			}
			// if we're here that means we got a work request - deal with that
			// sanity check - this should absolutely be true here
			if(  status.MPI_TAG != WUM_MPI_REQUEST_WU ){
				TR.Error << "ERROR: status.MPI_TAG != WUM_MPI_REQUEST_WU" << std::endl;
				break;
			}
		  TR.Trace << "Someone requested work!" << std::endl;
			send_next_WU_on_request();
		}

		if(  outbound().size() != before_size ){
			TR.Trace << "Present " << (before_size - outbound().size()) << " units" << std::endl;
		}

	  TR.Trace << "Now listening for any inbound work"  << std::endl;
		MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &result, &status);
		if( !result){ // if there's nothing on the line...
			if( !wait_until_message ){  // depending on this flag
				start_timer( TIMING_CPU );
				return;                 // return to caller
			} else {                    // or
				// try again
				continue;
			}
		}
		start_timer( TIMING_CPU );

		// interpret what's there
		switch( status.MPI_TAG ){
			case WUM_MPI_REQUEST_WU:
				TR.Debug << "Sending WU on request to " << status.MPI_SOURCE << std::endl;
				send_next_WU_on_request();
				break;
			case WUM_MPI_SEND_WU:
				TR.Debug << "Receiving WU from " << status.MPI_SOURCE << std::endl;
				receive_MPI_workunit();
				return; // now,  surely ther eis work
			default:
				TR.Error << "Unknown MPI_Message waiting from" << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;

//				// receive the message and discard result.
//				int data;
//				TR << "Cleaning out unknown data.." << std::endl;
//				MPI_Recv( &data, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
//				TR << "Cleaned out unknown data.." << std::endl;

				// can't break here because this could lock the master in an infinite loop
				// should we somehow clear the queue ?
				//return;
		}
		print_stats_auto();

	}
#endif
}





void MPI_WorkUnitManager::send_MPI_workunit( const WorkUnitBaseOP& MPI_ONLY(wu), int MPI_ONLY( dest_rank ) ) const {
#ifdef USEMPI
	const core::Real warning_threshold = 0.050000; // 50ms. usually sends take on the order of 600us so that limit is quite loose. if the limit is breached something is seriously going arwy.

	runtime_assert( dest_rank < mpi_npes() );
	runtime_assert( (int)dest_rank != mpi_rank() ); // no self-sending
	// serialize data
	wu->serialize();
	// now send data
	int size_of_raw_data;
	unsigned char * raw_data_ptr=NULL;
	size_of_raw_data = wu->raw_data_dump( &raw_data_ptr );


	// The first blocking send will only return when the master has recevied the data. Since the data transfer is negligible,
	// We'll count that as waiting time, not sending time
	start_timer( TIMING_WAIT );
	TRDEBUG << "Sending workunit to " << dest_rank << std::endl;
	// announce that you're about to send data and the size of it
	double  start_wait = get_time();
	MPI_Ssend( &size_of_raw_data,    1,                MPI_UNSIGNED, dest_rank, WUM_MPI_SEND_WU,    MPI_COMM_WORLD );
	TR.Debug << "Sent header announcing incoming WU of size " << F(5,1,(size_of_raw_data/1024.0)) << "kB" << " to node " << dest_rank << std::endl;
	start_timer( TIMING_TRANSFER_SEND );
	double  start_send = get_time();
	TR.Debug << "Sending WU" << std::endl;
	MPI_Ssend( (char*) raw_data_ptr, size_of_raw_data, MPI_CHAR,     dest_rank, WUM_MPI_DATA_BLOCK, MPI_COMM_WORLD );
	TR.Debug << "WU sent" << std::endl;
	double  end_send = get_time();
	send_wu_time_ += end_send - start_send;
	send_wu_time_n_++;
	//Unused parameter commented out to silence warning
	//double last_spent = start_timer( TIMING_CPU );
	TRDEBUG << "MPI_SEND: " << dest_rank << "  " << F(7,1,(start_send - start_wait)*1000.0)  << "ms  " << F(7,1,(end_send - start_send)*1000.0)
	<< "ms  " << F(5,1,(size_of_raw_data/1024.0)) << "kB" << std::endl;
	if( (end_send - start_wait) >  warning_threshold){
		TR << "WARNING LONG MPI_SEND: " << dest_rank << "  " << F(7,1,(start_send - start_wait)*1000.0)  << "ms  " << F(7,1,(end_send - start_send)*1000.0)
		   << "ms  " << F(5,1,(size_of_raw_data/1024.0)) << "kB" << std::endl;
	}

	delete [] raw_data_ptr;
	TR.Trace << "  Delete temp data.. " << std::endl;
	start_timer( TIMING_CPU );
	wu->clear_serial_data();

	traffic_total_sent_ += size_of_raw_data + sizeof( int );
#endif
}


void MPI_WorkUnitManager::receive_MPI_workunit( core::Size MPI_ONLY(node_rank) ){
#ifdef USEMPI
	MPI_Status status;
	int size_of_raw_data;
	unsigned char * raw_data_ptr;

	start_timer( TIMING_TRANSFER_RECV );

	double  start_recv = get_time();
	TR.Debug << "Receiving MPI header (ie how big the WU will be)" << std::endl;
	MPI_Recv( &size_of_raw_data, 1, MPI_UNSIGNED, node_rank, WUM_MPI_SEND_WU, MPI_COMM_WORLD, &status);
	TR.Debug << "Received MPI header, receiving WU of " << F(5,1,(size_of_raw_data/1024.0)) << "kB from node " << status.MPI_SOURCE << std::endl;
	raw_data_ptr = new unsigned char [size_of_raw_data];
	// now receive a datablock fromt he very same source
	TRDEBUG << "Confirmed datablock is coming" << std::endl;
	MPI_Recv( (char*) raw_data_ptr, size_of_raw_data, MPI_CHAR,  status.MPI_SOURCE, WUM_MPI_DATA_BLOCK, MPI_COMM_WORLD, &status);
	double  end_recv = get_time();
	TRDEBUG << "MPI_RECV: " << status.MPI_SOURCE  << F(5,0,(end_recv - start_recv)/1000000.0) << "us  " << F(5,1,(size_of_raw_data/1024.0)) << "kB" << std::endl;

	recv_wu_time_ += end_recv - start_recv;
	recv_wu_time_n_++;
	start_timer( TIMING_CPU );

	if( raw_data_ptr[size_of_raw_data-1] != 0){
		TR.Error << "  ERROR: cannot load data - terminal zero not found!" << std::endl;
		return;
	}

	raw_data_ptr[size_of_raw_data-1] = 0;

	traffic_total_received_ += size_of_raw_data + sizeof( int );

	TRDEBUG << "  RECEVIED WU: Data: " << std::endl;

	WorkUnitBaseOP wu( new WorkUnitBase );
  runtime_assert( wu != 0 );
	wu->raw_data_load( raw_data_ptr, size_of_raw_data );
	delete [] raw_data_ptr;

  // Here at this point we have a WorkUnitBaseOP to a workUnitBase.
  // Now we need to interpret the id field and upcast or somehow otherwise
  // create the right type of work unit such that the polymorphic code
  // for the interpretation of the serial data can take place.

	WorkUnitBaseOP qualified_wu = work_unit_list().get_work_unit( *wu )->clone();

  runtime_assert( qualified_wu != 0 );
	// cope over data (the header and the serial data)
	(*qualified_wu) = (*wu);
	(*qualified_wu).last_received_from_ = status.MPI_SOURCE;
	TRDEBUG << "  Received: " << std::endl;
	//if( TRDEBUG.visible() ) qualified_wu->print( TR );

	qualified_wu->deserialize( );
	qualified_wu->clear_serial_data();
	inbound().add( qualified_wu ); // add to stack of WUs to be processed
	TRDEBUG << "DONE Receiving" << std::endl;

	start_timer( TIMING_CPU );
#endif
}

void MPI_WorkUnitManager::send_next_WU_on_request( ){
#ifdef USEMPI
	MPI_Status status;
	int data;

	// Receive the actual data (it's been merely MPI_Probed up till now, just preceeding this function)
	MPI_Recv( &data, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, WUM_MPI_REQUEST_WU, MPI_COMM_WORLD, &status);

	// Find next work unit which does not blacklist the node that's requesting a workunit
	// (indicated by MPI_SOURCE)
	WorkUnitQueue::iterator suitable_work_unit = outbound().begin();
	for( ; suitable_work_unit != outbound().end(); ++suitable_work_unit )
	{
			// break out of loop once a WU is found that matches the above criterion
			if( !( (*suitable_work_unit)->in_blacklist( status.MPI_SOURCE ) ) ) break;
			// blurb some debug output if in debug mode
			TRDEBUG << "WU " << (*suitable_work_unit)->id() << " was not sent to " << status.MPI_SOURCE << " because it was blacklisted." << std::endl;
  }

	if( suitable_work_unit == outbound().end() ) {
				// No more work
				TRDEBUG << "No suitable work for node "  << status.MPI_SOURCE << " ( blacklisted=" << outbound().size() << ")" << std::endl;

				// craete a idling command workunit - since we have no work for the slave node that's requesting work
				WorkUnit_WaitOP wait_wu( new WorkUnit_Wait() );
				wait_wu->set_wu_type("waitwu");
				outbound().push_back( wait_wu );

				// set suitable_work_unit to the work unit just inserted
				suitable_work_unit = outbound().end();
				// after setting suitable_work_unit to outbound().end() we need to
				// decrement the iterator by one to have it point to the last element
				--suitable_work_unit;
	}

	// at this point there *must* be a work unit in the queue. if not we fucked up bad.
	runtime_assert( outbound().size() != 0 );

	TRDEBUG << "Sending next WU on request... "<< std::endl;
	start_timer( TIMING_CPU );

	// if we do, then suitable_work_unit one to the node that requested another job
	send_MPI_workunit( *suitable_work_unit, status.MPI_SOURCE );

	// remove the workunit that was just sent
	outbound().erase( suitable_work_unit );

	// if error free (ERROR CHECKING!)
	TRDEBUG << "END Send-on-request" << std::endl;

#endif
}



core::Real MPI_WorkUnitManager::start_timer( MPI_TIMING timing_mode ) const
{
	core::Real current_time = get_time();
	core::Real elapsed = 0;
	runtime_assert( timing_mode < TIMING_end );
	// analyse old timer
	if( timing_last_start_time_ != 0 ){
		elapsed = current_time - timing_last_start_time_;
		timing_total_[timing_last_type_] += elapsed;
	}

	// set new timer
	timing_last_start_time_ = current_time;
	timing_last_type_ = timing_mode;

	return elapsed;
}

void MPI_WorkUnitManager::print_stats_auto(){
	if( time(NULL) - last_stats_ > 60 ){
		MPI_WorkUnitManager::print_stats();
		last_stats_ = time(NULL);
	}
}

void MPI_WorkUnitManager::reset_timing_stats(){
	for( core::Size i=0;i<TIMING_end;i++) timing_total_[i] = 0;
}


long MPI_WorkUnitManager::wall_time() const{
	return time(NULL) - start_time_wall_clock_;
}

void MPI_WorkUnitManager::print_stats( )
{
	core::Real total_secs = 0;
	total_secs =
								timing_total_[TIMING_CPU] +
	              timing_total_[TIMING_TRANSFER_SEND] +
	              timing_total_[TIMING_TRANSFER_RECV] +
								timing_total_[TIMING_IO_READ] +
								timing_total_[TIMING_IO_WRITE] +
								timing_total_[TIMING_WAIT] +
								timing_total_[TIMING_IDLE];
	TR << "STATW" << get_machine_letter()  << " " <<
        I( (int)7, (int)wall_time() )  << "s " <<
        I( (int)7, total_secs )  << "s " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_CPU]           / total_secs)) << "% " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_TRANSFER_SEND] / total_secs)) << "% " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_TRANSFER_RECV] / total_secs)) << "% " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_IO_READ]       / total_secs)) << "% " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_IO_WRITE]      / total_secs)) << "% " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_WAIT]          / total_secs)) << "% " <<
				F( 4, 1, 100.0f*( timing_total_[TIMING_IDLE]          / total_secs)) << "% " <<
				F( 6, 2, (float(traffic_total_sent_)/1024.0/1024.0) ) << "Mb  " <<
				F( 6, 2, (float(traffic_total_received_)/1024.0/1024.0) ) << "Mb  " <<
				F( 6, 2, (float(send_wu_time_/1000.0/(send_wu_time_n_+1)) ) ) << "ms " <<
				F( 6, 2, (float(recv_wu_time_/1000.0/(recv_wu_time_n_+1)) ) ) << "ms " <<
				"";

	reset_timing_stats();

	core::Size in_total = inbound().size();
	core::Size in_total_structs=0;
	core::Size in_total_structs_memory=0;
	core::Size in_total_WU_memory=0;
	inbound().mem_stats( in_total_structs, in_total_structs_memory, in_total_WU_memory );

	core::Size out_total = outbound().size();
	core::Size out_total_structs=0;
	core::Size out_total_structs_memory=0;
	core::Size out_total_WU_memory=0;
	outbound().mem_stats( out_total_structs, out_total_structs_memory, out_total_WU_memory );

	TR <<
		"IWUs: "    << in_total  << "  " <<
		"IStruc: " << in_total_structs << "  " <<
		"IMem: "   << int((in_total_WU_memory + in_total_structs_memory)/1000.0) << "  " <<
		"OWUs: "    << out_total  << "  " <<
		"OStruc: " << out_total_structs << "  " <<
		"OMem: "   << int((out_total_WU_memory+out_total_structs_memory)/1000.0) << "  " <<
		"OMem: "   << outbound().mem_foot_print() <<
		" kB " <<
		std::endl;

}



} // namespace wum
} // namespace protocols
