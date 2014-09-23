// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

//testing memory

// Unit headers
#include <protocols/jd2/archive/MPIArchiveJobDistributor.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>
// AUTO-REMOVED #include <protocols/jd2/BatchJobInputter.hh> //for BOGUS_BATCH_ID
// Package headers
// AUTO-REMOVED #include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/jd2/MpiFileBuffer.hh>
#include <utility/io/ozstream.hh> //to toggle MPI rerouting

// Utility headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/assert.hh>
#include <basic/prof.hh>
#include <ObjexxFCL/string.functions.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/archive.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// C++ headers
#include <string>
// AUTO-REMOVED #include <ctime>
// AUTO-REMOVED #include <math.h>
#include <basic/prof.hh>
//Auto Headers
#include <utility/vector1.hh>

static thread_local basic::Tracer tr( "protocols.jd2.MPIArchiveJobDistributor" );
using basic::mem_tr;

namespace protocols {
namespace jd2 {
namespace archive {

///our setup of dedicated processes...
int const in_master_rank_( 1 ); //keep const for now
int const in_file_buf_rank_( 0 );
int const in_archive_rank_( 2 );
int const in_min_client_rank_( 3 );

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;


///@details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
MPIArchiveJobDistributor::MPIArchiveJobDistributor() :
  MPIFileBufJobDistributor( in_master_rank_, in_file_buf_rank_, in_min_client_rank_, true /*start empty*/ ),
	nr_notify_( option[ OptionKeys::archive::completion_notify_frequency] ),
	archive_rank_( in_archive_rank_ )
{

	//if we are testing we want to send JOB_COMPLETION more often
	if ( option[ OptionKeys::run::test_cycles ] || option[ OptionKeys::run::dry_run ] ) {
		nr_notify_ = std::min( nr_notify_, Size(10) );
	}
}

void
MPIArchiveJobDistributor::set_archive( ArchiveBaseOP archive ) {
	if ( rank() == archive_rank() ) {
		theArchive_ = archive;
	}
}
///@brief dummy for master/slave version -- start the appropriate process depending on rank()
void
MPIArchiveJobDistributor::go( protocols::moves::MoverOP mover )
{
	//copied MPIFileJobDistributor, because in this case the archive - process sends stop to FileBuf.
	utility::io::ozstream::enable_MPI_reroute( min_client_rank(), file_buf_rank() );
	mem_tr << "MPIArchiveJobDistributor::go" << std::endl;
	/// JD
	if ( rank() == master_rank() ) {
		tr.Warning << "Master JD starts" << std::endl;
    master_go( mover );
	} else if ( rank() == file_buf_rank() ) {
		/// FileBuffer
		protocols::jd2::WriteOut_MpiFileBuffer buffer( file_buf_rank() );
		tr.Warning << "FileBuffer starts " << std::endl;
		buffer.run();
	} else if ( rank() == archive_rank() ) {
		/// Archive
		tr.Warning << "Archive starts... " << std::endl;
		archive::ArchiveManager archive( archive_rank(), master_rank(), file_buf_rank() );
		runtime_assert( theArchive_ != 0 );
		archive.go( theArchive_ );
		tr.Warning << "send STOP to FileBuffer " << std::endl;
		protocols::jd2::WriteOut_MpiFileBuffer buffer( file_buf_rank() );
		buffer.stop();
	} else if( rank() >= min_client_rank() ){
		/// Slave/Runner/Worker
		go_main( mover );
  }

	// ideally these would be called in the dtor but the way we have the singleton pattern set up the dtors don't get
	// called
#ifdef USEMPI
 	MPI_Barrier( MPI_COMM_WORLD );
 	MPI_Finalize();
#endif
	if ( rank() == master_rank() ) {
		std::cerr << "MPI FINALIZED closing down... " << std::endl;
		std::cout << "MPI FINALIZED closing down... " << std::endl;
	}
}

///@detail receive a new batch from ArchiveManager -- interpret batch_nr == 0 as STOP
bool
MPIArchiveJobDistributor::receive_batch( Size MPI_ONLY( source_rank ) ) {
basic::prof_show();
#ifdef USEMPI
	MPI_Status status;
	int buf[ 2 ];
	//receive size of string
	MPI_Recv( buf, 2, MPI_INT, source_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &status );
	Size size( buf[ 0 ]);
	Size id( buf[ 1 ] );
	//receive string
	std::string new_batch;
	char *cbuf = new char[ size+1 ];
	MPI_Recv( cbuf, size, MPI_CHAR, source_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &status );

	//STOP?
	if ( id == 0 ) { //use this as STOP signal!
		tr.Debug << "received STOP signal from Archive " << std::endl;
		delete[] cbuf;
		return false;
	}

	///assign C++ string to cbuf
	new_batch.assign( cbuf, size );
	delete[] cbuf;

	tr.Info << "received new batch " << new_batch << " with id " << id << std::endl;
	add_batch( new_batch, id );
#endif
	return true;
}

///@detail sync batches with worker nodes.. this is called if they get a job for a batch they don't know yet...
/// this method will send ALL batches they don't have yet.
void
MPIArchiveJobDistributor::sync_batches( Size MPI_ONLY( slave_rank ) ) {
	PROF_START( basic::ARCHIVE_SYNC_BATCHES );
#ifdef USEMPI
	tr.Trace << "Node " << rank() << " sync batches with " << slave_rank << std::endl;
	int buf[ 4 ];
	buf[ 1 ] = ADD_BATCH;
	MPI_Status status;

	///send last known batch from SLAVE --> MASTER
	Size slave_batch_size( nr_batches() );
	Size nr_to_have;
	if ( rank() != master_rank() ) { //SLAVE -- SEND
		buf[ 0 ] = slave_batch_size;
		MPI_Send( &buf, 1, MPI_INT, master_rank(), MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	} else {  //MASTER -- RECEIVE
		MPI_Recv( &buf, 1, MPI_INT, slave_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &status );
		slave_batch_size = buf[ 0 ];
	}

	tr.Trace << "Node " << rank() << " slave_batch_size " << slave_batch_size << std::endl;

	//MASTER --> SLAVE how many batches will be sent
	nr_to_have = nr_batches();
	if ( rank() != master_rank() ) { //SLAVE
		MPI_Recv( &buf, 1, MPI_INT, master_rank(), MPI_JOB_DIST_TAG, MPI_COMM_WORLD, &status );
		nr_to_have = buf[ 0 ];
	} else {  //MASTER
		buf[ 0 ] = nr_to_have;
		MPI_Send( &buf, 1, MPI_INT, slave_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	}
	tr.Trace << "Node " << rank() << " master_batch_size " << nr_to_have << std::endl;

	//MASTER --> SLAVE now send the individual batches
	for ( Size send_id = slave_batch_size + 1; send_id <= nr_to_have; ++send_id ) {
		if ( rank() != master_rank() ) { //SLAVE
			receive_batch( master_rank() );
			tr.Trace << "nr_batches() " << nr_batches() << " send_id " << send_id << std::endl;
			runtime_assert( nr_batches() == send_id );
		} else {  //MASTER
			//send size of string
			buf[ 0 ] = batch( send_id ).size();
			buf[ 1 ] = send_id;
			MPI_Send(buf, 2, MPI_INT, slave_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
			//send string
			MPI_Send(const_cast<char*> ( batch( send_id ).data()), batch( send_id ).size(), MPI_CHAR, slave_rank, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
		}
	}
	#endif
	PROF_STOP( basic::ARCHIVE_SYNC_BATCHES );
}

///@detail send message to ArchiveManager .. eg. QueueEmpty
/// always send current_batch_id with the message ... used to determine if QueueEmpty is outdated
void
MPIArchiveJobDistributor::master_to_archive( Size MPI_ONLY(tag) ) {
#ifdef USEMPI
	runtime_assert( rank() == master_rank() );
	runtime_assert( rank() != archive_rank() );
	Size const mpi_size( 6 );
	int mpi_buf[ mpi_size ];
	mpi_buf[ 0 ] = tag;
	mpi_buf[ 1 ] = current_batch_id();
	MPI_Send( &mpi_buf, mpi_size, MPI_INT, archive_rank(), MPI_ARCHIVE_TAG, MPI_COMM_WORLD );
#endif
}

///@detail called if JD is at and of BatchQueue...
/// for a worker node that might mean he needs to sync batches with Master
/// for a master node it means he sends QUEUE-EMPTY to ArchiveManager
void
MPIArchiveJobDistributor::batch_underflow() {
	if ( !( rank() == master_rank() ) ) {
		slave_to_master( BATCH_SYNC );
		sync_batches( rank() );
	} else if ( rank() == master_rank() ) {
		PROF_START( basic::MPI_JD2_WAITS_FOR_ARCHIVE );
		tr.Debug << "no more batches... ask ArchiveManager if there is some more to do... wait..." << std::endl;
		_notify_archive();
		basic::show_time( tr,  "no more batches: send QUEUE_EMPTY to archive" );
		master_to_archive( QUEUE_EMPTY );
		tr.Info << "wait for answer on QUEUE-EMPTY msg... send with " << current_batch_id() << " batch_id " << std::endl;
		eat_signal( ADD_BATCH, archive_rank() );
		receive_batch( archive_rank() ); //how about some time-out
		tr.Debug << "...received " << std::endl;
		basic::show_time( tr,  "refilled queue: received new batches after QUEUE_EMPTY" );
		PROF_STOP( basic::MPI_JD2_WAITS_FOR_ARCHIVE );
	}
}

///@detail process messages... BATCH_SYNC, ADD_BATCH, or delegate to Parent class
/// also send pending CompletionMessages out...
////this is a good place to do it, since CompletionMessages are non-blocking and we are otherwise in blocking communication with WorkerNodes
bool
MPIArchiveJobDistributor::process_message(
   core::Size msg_tag,
	 core::Size slave_rank,
	 core::Size slave_job_id,
	 core::Size slave_batch_id,
	 core::Real run_time
) {
	runtime_assert( rank() == master_rank() );

	//	basic::show_time( tr,  "jd2 main msg-loop: process message..." );

	// send out any pending notifications to archive if present -- this is non-blocking
	_notify_archive(); //we should get here often enough... (basically every finished job
	//-- unless of course we haven't started any jobs yet)

	// now go thru messages
	switch ( msg_tag ) {
	case BATCH_SYNC: //a slave has received a job with an unknown batch and asks for an update of the Batchlist
		sync_batches( slave_rank );
		break;
	case ADD_BATCH: //the ArchiveManager adds a new BATCH
		runtime_assert( slave_rank == archive_rank() );
		receive_batch( archive_rank() );
		break;
	case CANCEL_BATCH: //the ArchiveManager cancels a certain BATCH
		runtime_assert( slave_rank == archive_rank() );
		{
			//bool was_good( get_current_batch() != jd2::BatchJobInputter::BOGUS_BATCH_ID );
			tr.Trace << "currently running batch " << get_current_batch() << std::endl;
			receive_batch( archive_rank() ); //right now we hack it... ArchiveSends BOGUS_BATCH_ID 5 to cancel batch 5
			// now in the batch_list the name of the respective batch willl be changed to BOGUS_BATCH_ID ... next time we
			// get to it we will skip, if we are running it now, we will terminate.
			tr.Trace << "now the current batch is " << get_current_batch() << std::endl;
			// put this functionality into Baseclass: JobDistributor::obtain_new_job -- it was probably causing a dead-lock in the communication...
			//			if ( was_good  && get_current_batch() == jd2::BatchJobInputter::BOGUS_BATCH_ID ) next_batch();
		}
		break;
	default:
		return Parent::process_message( msg_tag, slave_rank, slave_job_id, slave_batch_id, run_time );
	}

	return true;
}

///@detail queue up a CompletionMessage
void
MPIArchiveJobDistributor::notify_archive( CompletionMessage const& msg ) {
	//TODO: check if there are older messages regarding this batch... if so ... remove
	tr.Debug << "add to notification queue " << msg.batch_id << std::endl;
	if ( pending_notifications_.size()
		&& pending_notifications_.back().batch_id == msg.batch_id
		&& pending_notifications_.back().msg_tag == msg.msg_tag
	) {
		pending_notifications_.back() = msg;
	} else {
		pending_notifications_.push_back( msg );
	}
}

/// stuff needed for non-blocking communication
#ifdef USEMPI
MPI_Request notify_request;
int notify_buf[ 6 ];
bool notify_first( true ); ///buffer hasn't been used yet?
#endif

///the private implementation of notify_archive.
/// send JOB_COMPLETION message to Archive if a message is in the message queue.
void MPIArchiveJobDistributor::_notify_archive() {
	PROF_START( basic::MPI_NOTIFY_ARCHIVE );
	basic::Tracer notification_tracer( "protocols.jd2.notifications" );

	//nothing in queue?
	if ( pending_notifications_.size() == 0 ) return;

	//okay, queue is filled send something
#ifdef USEMPI
	int flag( 1 );
	if ( !notify_first ) { //if not first time, make sure last message has been received already...
		notification_tracer.Debug << "test MPI-Send completion of last JOB_COMPLETION ( batch_" << notify_buf[ 1 ] << " ) message...";
		basic::show_time( tr,  "try to send JOB_COMPLETION" );
		MPI_Status status;
		MPI_Test( &notify_request, &flag, &status ); //has last communication succeeded ? --- buffer is free again.
		int flag2;
		MPI_Test_cancelled( &status, &flag2 );
		notification_tracer.Debug << ( flag ? "completed " : "pending " ) << ( !flag2 ? "/ test succeeded " : "/ test cancelled" ) << std::endl;
	}
	if ( flag ) {
		//okay ready to send next message
		CompletionMessage const& msg(	pending_notifications_.front() );
		notification_tracer.Debug << "send out JOB_COMPLETION " << msg.batch_id << std::endl;
		basic::show_time( tr,  "send JOB_COMPLETION" );
		//		int notify_buf[ 6 ];
		notify_buf[ 0 ] = msg.msg_tag;//JOB_COMPLETION, QUEUE_EMPTY;
		notify_buf[ 1 ] = msg.batch_id;
		notify_buf[ 2 ] = msg.final ? 1 : 0;
		notify_buf[ 3 ] = msg.bad;
		notify_buf[ 4 ] = msg.good;
		notify_buf[ 5 ] = msg.njobs;
		MPI_Isend( &notify_buf, 6, MPI_INT,  archive_rank(), MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &notify_request ); //don't block JobDistributor
		pending_notifications_.pop_front();
				notify_first = false;
	}
	basic::show_time( tr,  "finished _notify_archive" );
#endif
	PROF_STOP( basic::MPI_NOTIFY_ARCHIVE );
}

///@detail work out if CompletionMessage should be send... looks at completed/bad decoys
/// send "final" message if all jobs done... sends "update" message if nr_new_completed_ > nr_notify
void
MPIArchiveJobDistributor::notify_archive( core::Size batch_id ) {
	/// send "final" message ?
	tr.Trace << "notify_archive for batch: " << batch_id << " now " << nr_new_completed_[ batch_id ] << " decoys " << std::endl;
	if ( nr_completed_[ batch_id ] + nr_new_completed_[ batch_id ] + nr_bad_[ batch_id ] == nr_jobs_[ batch_id ] ) {
		nr_completed_[ batch_id ] += nr_new_completed_[ batch_id ];
		nr_new_completed_[ batch_id ] = 0;
		//	still send to close files in MPI-FILE-BUF	if ( batch( batch_id ) != BatchJobInputter::BOGUS_BATCH_ID ) { //don't send message if this Batch has ben CANCELLED
		notify_archive(  CompletionMessage( batch_id, true, nr_bad_[ batch_id ], nr_completed_[ batch_id ], nr_jobs_[ batch_id ] ) );
			//		}
			//// send "update" message ?
	} else if ( nr_new_completed_[ batch_id ] >= nr_notify_ ) {
		nr_completed_[ batch_id ] += nr_new_completed_[ batch_id ];
		nr_new_completed_[ batch_id ] = 0;
		//if ( batch( batch_id ) != BatchJobInputter::BOGUS_BATCH_ID ) { //don't send message if this Batch has ben CANCELLED
		notify_archive( CompletionMessage( batch_id, false, nr_bad_[ batch_id ], nr_completed_[ batch_id ], nr_jobs_[ batch_id ] ) );
			//		}
	}
	tr.Trace << "nr_batches " << nr_batches() << " current_job_id() " << current_job_id() << " get_jobs().size() " << get_jobs().size()
					 << " nr_processors " << number_of_processors() << std::endl;
 	//are we quickly running out of jobs? -- checking for equality to reduce number of messages -- is this safe? do we ever skip jobs?
	if ( nr_batches() == batch_id && ( (int) current_job_id() == ( (int) get_jobs().size() - (int) number_of_processors() ) ) ) {
		//tr.Info << "jobs are low... send QUEUE_EMPTY with " << batch_id << " batch_id " << std::endl;
		//		pending_notifications_.push_front( CompletionMessage( batch_id, QUEUE_EMPTY ) );
		//		_notify_archive();
	}
}

///@detail overloaded to update our job-statistics ( needed for CompletionMessages )
void MPIArchiveJobDistributor::mark_job_as_completed( core::Size job_id, core::Size batch_id, core::Real run_time ) {
	tr.Trace << "mark_job_as_completed " << job_id << " batch: " << batch_id << " " << run_time << " seconds" << std::endl;
	Parent::mark_job_as_completed( job_id, batch_id, run_time );
	if ( rank() == master_rank() ) {
		runtime_assert( batch_id <= nr_jobs_.size() );
		nr_new_completed_[ batch_id ] += 1;
		notify_archive( batch_id );
	}
}

void MPIArchiveJobDistributor::mark_job_as_bad( core::Size job_id, core::Size batch_id ) {
	Parent::mark_job_as_bad( job_id, batch_id );
	if ( rank() == master_rank() ) {
		runtime_assert( batch_id <= nr_jobs_.size() );
		nr_bad_[ batch_id ] += nstruct_[ batch_id ];
		notify_archive( batch_id );
	}
}

///@detail load new batch from BatchQueue .. overloaded to setup the statistics for CompletionMessages
void MPIArchiveJobDistributor::load_new_batch() {
	//	if ( current_batch_id() )	notify_archive( current_batch_id() );
	Parent::load_new_batch();
	if ( rank() == master_rank() ) { //in principle I'd rather do this in add_batch() but we need option of new batch for nstruct...
		while( nr_jobs_.size() < current_batch_id() ) {
			nr_jobs_.push_back( get_jobs().size() );
			nr_new_completed_.push_back( 0 );
			nr_completed_.push_back( 0 );
			nr_bad_.push_back( 0 );
			nstruct_.push_back( option[ out::nstruct ] ); ///Assumming that all JobInputters create nstruct jobs per input_tag...
		}
		runtime_assert( nr_jobs_.size() == current_batch_id() );
		mem_tr << "MPIArchiveJobDistributor::load_new_batch()'ed" << std::endl;
	}
}



}//archive
}//jd2
}//protocols
