// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MPIFileBufJobDistributor.cc
/// @brief  implementation of MPIFileBufJobDistributor
/// @author Oliver Lange olange@u.washington.edu


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

// Unit headers
#include <protocols/jd2/archive/ArchiveManager.hh>
#include <protocols/jd2/archive/MPIArchiveJobDistributor.hh>
#include <protocols/jd2/archive/ArchiveBase.hh>
#include <protocols/jd2/BatchJobInputter.hh> //for BOGUS_BATCH_ID

// Package headers
#include <protocols/jd2/MpiFileBuffer.hh>


//for factory
//#include <protocols/abinitio/IterativeAbrelax.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <ObjexxFCL/string.functions.hh>
#include <utility/file/file_sys_util.hh>

// Option headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

//#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/option.cc.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/prof.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <sstream>
#include <iterator>

//Debug headers
#include <fstream> //testing

#if (defined WIN32) //&& (!defined WIN_PYROSETTA)
#include <windows.h>
#endif

//Auto Headers
#include <protocols/checkpoint/CheckPointer.hh>
//#include <protocols/jobdist/Jobs.hh>
//#include <protocols/noesy_assign/CrossPeak.hh>
#include <utility/vector1.hh>
#include <boost/bind.hpp>

static basic::Tracer tr( "protocols.jd2.Archive" );
static basic::MemTracer mem_tr;

OPT_1GRP_KEY( File, iterative, input_pool )
OPT_1GRP_KEY( String, iterative, input_pool_struct_type )

bool protocols::jd2::archive::ArchiveManager::options_registered_( false );

using namespace basic::options;
using namespace basic::options::OptionKeys;
//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::jd2::archive::ArchiveManager::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		NEW_OPT( iterative::input_pool, "read these structures into pool", "" );
		NEW_OPT( iterative::input_pool_struct_type, "specifies the input-silent-struct type", "protein" );
	}
}

namespace protocols {
namespace jd2 {
namespace archive {

#ifdef WIN32
void sleep(int seconds){
	//#if (defined WIN32) && (!defined WIN_PYROSETTA)
	Sleep( seconds * 1000 );
	//#endif
}
#endif

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;


std::string
Batch::batch() const {
	return "batch_" + ObjexxFCL::lead_zero_string_of( id(), 6 );
}

std::string Batch::dir() const {
	return batch() + "/";
}

std::string Batch::silent_out() const {
	return batch() + "/decoys.out";
}

std::string Batch::alternative_decoys_out() const {
	return batch() + "/decoys_stage2.out";
}

std::string Batch::silent_in() const {
	// if ( has_silent_in() )
	return batch() + "/decoys.in";
	// else
	//  return "";
}

std::string Batch::score_file() const {
	return batch() + "/score.fsc";
}

std::string Batch::flag_file() const {
	return batch() + "/flags";
}

std::string Batch::broker_file() const {
	return batch() + "/setup.tpb";
}

std::string Batch::extra_broker_files() const {
	return "";
}

void Batch::show( std::ostream& out, bool single_line ) const {
	std::string eol( single_line ? " " : "\n" );
	out << "ID " << id() << eol
		<< "INPUT " << ( has_silent_in() ? "yes" : "no" ) << eol
		<< "NSTRUCT " << nstruct() << eol
		<< "INTERMEDIATES " << ( intermediate_structs() ? "yes" : "no" ) << eol
		<< "RETURNED " << decoys_returned() << eol
		<< "FINISHED " << ( has_finished() ? "yes" : "no" ) << eol
		<< "CANCELLED " << ( is_cancelled() ? "yes" : "no" ) << eol
		<< "ALLOW_READING_CANCELLED_DECOYS " << ( allow_reading_cancelled_decoys() ? "yes" : "no" ) << eol;
}

std::ostream& operator<< (std::ostream& out, Batch const& batch ) {
	batch.show( out, true );
	return out;
}

void Batch::write_info_file() const {
	utility::io::ozstream out( dir() + "BATCH_INFO" );
	tr.Debug << "write batch info " << dir() << "BATCH_INFO" << std::endl;
	show( out, false /*not single_line*/ );
}

void Batch::read_info_file() {
	core::Size this_id = id(); //to detec errors
	std::string this_batch = batch(); //for error report
	utility::io::izstream in( dir() + "BATCH_INFO" );
	if ( !in.good() ) throw CREATE_EXCEPTION(EXCN_Archive, "cannot find " + dir() + "BATCH_INFO" );
	in >> *this;
	if ( this_id != id() ) {
		throw CREATE_EXCEPTION(EXCN_Archive, "Inconsistency detected when reading BATCH_INFO for "+ this_batch+" ID in BATCH_INFO is " + batch() );
	}
}

//instead of goto statements:   I think goto would be clearer... but there are coding guidlines to adhere...
void report_tag_error( Batch& batch, std::string const& expected_tag, std::string const& tag ) {
	throw CREATE_EXCEPTION(EXCN_Archive, "Error reading batch information for batch: "+batch.batch()+" expected_tag: "+expected_tag+ " found " + tag);
}

void report_value_error( Batch& batch, std::string const& tag ) {
	throw CREATE_EXCEPTION(EXCN_Archive, "Error reading batch information for batch: "+batch.batch()+" wrong value for tag: "+tag );
}

std::istream& operator >> (std::istream& in, Batch &batch ) {
	std::string tag;
	std::string expected_tag;

	in >> tag;
	expected_tag = "ID";
	if ( tag == expected_tag ) {
		in >> batch.batch_id_;
		if ( !in.good() ) report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );

	in >> tag;
	expected_tag = "INPUT";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.has_silent_in_ = true;
		else if ( yesno == "no" ) batch.has_silent_in_ = false;
		else report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );

	in >> tag;
	expected_tag = "NSTRUCT";
	if ( tag == expected_tag ) {
		in >> batch.nstruct_;
		if ( !in.good() ) report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );

	in >> tag;
	expected_tag = "INTERMEDIATES";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.intermediate_structs_ = true;
		else if ( yesno == "no" ) batch.intermediate_structs_ = false;
	} else report_tag_error( batch, expected_tag, tag );

	in >> tag;
	expected_tag = "RETURNED";
	if ( tag == expected_tag ) {
		in >> batch.decoys_returned_to_archive_;
		if ( !in.good() ) report_value_error( batch, tag );
	} else report_tag_error( batch, expected_tag, tag );

	in >> tag;
	expected_tag = "FINISHED";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.has_finished_ = true;
		else if ( yesno == "no" ) batch.has_finished_ = false;
	} else report_tag_error( batch, expected_tag, tag );

	in >> tag;
	expected_tag = "CANCELLED";
	if ( tag == expected_tag ) {
		std::string yesno;
		in >> yesno;
		if ( !in.good() ) report_value_error( batch, tag );
		if ( yesno == "yes" ) batch.is_cancelled_ = true;
		else if ( yesno == "no" ) batch.is_cancelled_ = false;
	} else report_tag_error( batch, expected_tag, tag );

	return in;
}

//#ifndef WIN32

void BaseArchiveManager::set_archive( AbstractArchiveBaseOP anArchive ) {
	theArchive_ = anArchive;
	theArchive_->set_manager( this );
	theArchive_->initialize();
}

/// @details constructor.  Notice it calls the parent class!  It also builds some internal variables for determining
///which processor it is in MPI land.
ArchiveManager::ArchiveManager( core::Size archive_rank, core::Size jd_master_rank, core::Size file_buf_rank ) :
	archive_rank_( archive_rank ),
	jd_master_rank_( jd_master_rank ),
	file_buf_rank_( file_buf_rank ),
	save_archive_time_interval_( 60 )
{
	runtime_assert( options_registered_ );
}

core::Size ArchiveManager::unfinished_batches() const {
	Size unfinished_batches( 0 );
	for ( auto const & it : batches() ) {
		if ( !it.has_finished() && !it.is_cancelled() && it.valid() ) ++unfinished_batches;
	}
	return unfinished_batches;
}


void
ArchiveManager::go( ArchiveBaseOP archive )
{
	set_archive( archive );
	tr.Debug << "starting ArchiveManager ..." << archive_rank_ << " " << jd_master_rank_ << " " << file_buf_rank_ << std::endl;
	mem_tr << "initialized IterativeAbrelax" << std::endl;
	try {
		if ( !restore_archive() ) {
			if ( option[ OptionKeys::iterative::input_pool ].user() ) {
				std::string const& decoys( option[ OptionKeys::iterative::input_pool ]() );
				tr.Info << "reading decoys from " <<  decoys << " into archive " << std::endl;
				core::io::silent::SilentFileOptions opts;
				core::io::silent::SilentFileData sfd( decoys, false, false,  option[ OptionKeys::iterative::input_pool_struct_type ](), opts );
				sfd.read_file( decoys );
				the_archive().init_from_decoy_set( sfd );
			}
		}
		save_archive();
		read_existing_batches();
	} catch ( utility::excn::Exception& excn ) {
		send_stop_to_jobdistributor();
		throw;
	}
	// if ( batches_.size() == 0 ) theArchive_->generate_batch();
	sleep( 5 ); //give JobDistributor time to start up...
#ifdef USEMPI
	MPI_Status status;
	//MPI_Request request;
#endif
	bool stop( false );
	bool print_status( true );
	while ( !stop || unfinished_batches() ) {

		if ( print_status && tr.Debug.visible() ) {
			tr.Debug << "probing for message in ArchiveManager" << std::endl;
			tr.Debug << "\nSTATUS: " << (stop ? "STOP send: " : "" ) << "  ------ unfinished_batches: " << unfinished_batches() << std::endl;
			tr.Debug << "POOL_STATUS: " << std::endl;
			the_archive().save_status( tr.Debug );
			tr.Debug << "END_STATUS\n\n"<< std::endl;
			basic::show_time( tr,  "manager main msg-loop: probe for message..." );
			print_status = false;
		}
		//is there a message ?
		// AMW: cppcheck flags this as being reducible in scope
		// but MPI needs it to be here
		int flag( -1 );
#ifdef USEMPI
		//no idea why... but 4 request seems to be the magical number... to receive the correct answer ... WEIRD
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
#endif
		//  if ( !flag ) { //nothing ...
		//    //tell JobDistributor, that we are ready to receive message
		//    int buf[ 4 ];
		//    buf[ 0 ] = NOTIFICATION_QUERY;
		//    MPI_Send( &buf, 1, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
		//    sleep( 1 );
		//    //check if there is something this time...
		//    MPI_Iprobe( jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &flag, &status );
		//   }

		try {
			//if there is a message -- go get it.
			int buf[ 6 ]={0,0,0,0,0,0};
			if ( flag ) {
#ifdef USEMPI
				int merrno = MPI_Recv( &buf, 6, MPI_INT, jd_master_rank_, MPI_ARCHIVE_TAG, MPI_COMM_WORLD, &status );
				if ( merrno != MPI_SUCCESS ) tr.Error << "MPI_Recv error " << std::endl;
#endif
				//    basic::show_time( tr,  "manager main msg-loop: received message..." );
			} else { //nothing received
				//    basic::show_time( tr,  "manager main msg-loop: no message: idle..." );
				idle();
				continue;
			}
			print_status = true;
			// here if we got a message
			Size const msg_tag( buf[ 0 ]);
			tr.Debug << "received message in ArchiveManager " << msg_tag << std::endl;

			switch( msg_tag ) {
			case JOB_COMPLETION : {
				Size const batch_id( buf[ 1 ] );
				bool const final( buf[ 2 ] == 1 );
				Size const bad( buf[ 3 ] );
				Size const good( buf[ 4 ] );
				Size const total( buf[ 5 ] ); //total nr of jobs
				basic::show_time( tr,  "ArchiveManager receveid job-completion..." );
				tr.Debug << "ArchiveManager received JOB_COMPLETION " << batch_id << " " << bad << " " << good << " " << total << std::endl;
				jobs_completed_[ batch_id ] = CompletionMessage( batch_id, final, bad, good, total );
				break; //switch
			}
			case QUEUE_EMPTY : {
				Size const batch_id( buf[ 1 ] );

				//we ignore QUEUE_EMPTY if we know that a new batch has been submitted after issuing of this signal (i.e., the batch-number
				// coming with the message would be smaller than the currently highest batch number... however, there might be invalid batches...
				// find last valid and unfinished batch in list:
				Size max_working_batch_id( batches_.size() );
				if ( batches_.size() ) {
					while ( max_working_batch_id > 0
							&& ( !batches_[ max_working_batch_id ].valid() || batches_[ max_working_batch_id ].has_finished() ) )
							--max_working_batch_id;
					if ( batch_id <= max_working_batch_id ) {
						tr.Info << "ArchiveManager ignored outdated QUEUE_EMPTY with batch_id " << batch_id << " -- already submitted " << batches_.size() << std::endl;
						break; //switch
					}
				}
				//any job-completions we should work thru before generating a new batch?
				PROF_START( basic::ARCHIVE_CRITICAL_JOBSCOMPLETE );
				while ( jobs_completed_.size() ) {
					jobs_completed(); //get thru these before making job decisions
				}
				PROF_STOP( basic::ARCHIVE_CRITICAL_JOBSCOMPLETE );
				//  the_archive().idle(); why was this in the job-completed loop ?

				PROF_START( basic::ARCHIVE_GEN_BATCH );
				//this is a valid QUEUE_EMPTY request: do something about it
				tr.Info << "ArchiveManager received QUEUE_EMPTY" << std::endl;
				tr.Debug << "JD batch_id: " << batch_id << " max_working_batch_id: " << max_working_batch_id << std::endl;
				basic::show_time( tr,  "manager main msg-loop: queue empty..." );
				if ( !the_archive().finished() ) {
					//if !finished Archive should always generate a batch...
					//but let's make sure by monitoring, since it would be bad if we hang in the communication...
					Size ct( batches_.size() );//monitor number of batches
					if ( !stop ) the_archive().generate_batch();
					if ( ct == batches_.size() ) { //if generate_batch didn't create anything --- we still owe Jobdistributor a signal
						send_stop_to_jobdistributor(); //send stop
						stop = true;
					}
				} else {
					tr.Debug << "archive is finished ... spinning down" << std::endl;
					send_stop_to_jobdistributor();
					stop = true;
				}
				PROF_STOP( basic::ARCHIVE_GEN_BATCH );
				basic::prof_show();
				break; //switch
			}
			default :
				utility_exit_with_message( "unknown msg in ArchiveManager " + ObjexxFCL::string_of( msg_tag ) );
			} //switch
		} catch ( utility::excn::Exception &excn ) {
			basic::show_time( tr,  "Exception in main msg-loop !" );
			tr.Error << excn.msg() << std::endl;
			tr.Error << "spinning down" << std::endl;
			save_archive();
			//this usually doesn't work the jobs always run to completion ... let's hard exit for now.
			utility_exit_with_message( "error detected in ArchiveManager -- spinning down" );
			send_stop_to_jobdistributor();
			stop = true;
		}
	} //while loop
	save_archive();
	tr.Info << "ArchiveManager finished !!!" << std::endl;
}

void
ArchiveManager::idle() {

	{ //save archive
		static time_t last_save( time(nullptr) );
		time_t now( time( nullptr ) );
		Size const elapsedtime( now - last_save );
		if ( elapsedtime > save_archive_time_interval_ ) {
			save_archive();
			last_save = now;
		}
	}

	// tr.Debug << "idle..." << std::endl;
	if ( jobs_completed_.size() ) {
		PROF_START( basic::ARCHIVE_JOBSCOMPLETE );
		jobs_completed();
		PROF_STOP( basic::ARCHIVE_JOBSCOMPLETE );
		return;
	};

	// if ( !the_archive().finished() && the_archive().ready_for_batch() ) {
	//  the_archive().generate_batch();
	// } else {
	time_t before( time(nullptr) );
	the_archive().idle();
	time_t after( time( nullptr ) );
	if ( after-before > 1 ) tr.Debug << "spend " << after-before << " seconds in archives idle method... " << std::endl;
	//sleep some more if idle didn't use much time
	if ( after-before < 5 ) sleep( (5 - ( after - before )) );
	// }
}

void BaseArchiveManager::read_returning_decoys( Batch& batch, bool final ) {
	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileData sfd( opts );
	utility::vector1< std::string > tags_in_file;

	//this keeps order as in file... important since we skip already known tags by just keeping their number
	sfd.read_tags_fast( batch.silent_out(), tags_in_file );

	unlock_file( batch, final );

	tr.Debug << "found " << tags_in_file.size() << " decoys in " << batch.silent_out() << std::endl;

	auto iter = tags_in_file.begin();
	std::string unread_tag = "none";
	for ( Size ct = 1;
			iter != tags_in_file.end() && ct <= batch.decoys_returned();
			++iter, ++ct ) {
		unread_tag = *iter;
	}; //just skipping...
	utility::vector1< std::string > tags_to_read;

	std::copy( iter, tags_in_file.end(), std::back_inserter( tags_to_read ) );
	if ( tags_to_read.size() ) {
		std::cerr << "last_skipped tag: " << unread_tag << " first tag to read: " << tags_to_read[1] << " last tag to read: " << tags_to_read.back() << std::endl;
		try {
			sfd.read_file( batch.silent_out(), tags_to_read );
		} catch ( utility::excn::Exception& excn ) { //or should we be more specific ?
			if ( final ) throw; //rethrow if it is the final version of the file...
			tr.Error << "ignored ERROR: " << excn.msg() << std::endl;
			tr.Error << "this is not the final version of " << batch.silent_out() << "\n... maybe some data is still held in a cache of the filesystem..."
				<< " let's see if it works better the next time we have to read" << std::endl;
			//or sleep( 5 ) and retry as above ?
			return;
		}

		{ //now update our batch information so that this is already known in read_structures
			runtime_assert( batch.decoys_returned()+sfd.size() == tags_in_file.size() );
			batch.set_decoys_returned( tags_in_file.size() );
			if ( final ) {
				batch.mark_as_finished();
			}
		}

		PROF_STOP( basic::ARCHIVE_READ_DECOYS );
		tr.Debug << "add " << tags_to_read.size() << " structures to archive " << std::endl;

		PROF_START( basic::ARCHIVE_EVAL_DECOYS );

		//if alternative decoys present read also these
		SilentFileData alternative_decoys_sfd( opts );
		if ( batch.intermediate_structs() ) {
			alternative_decoys_sfd.read_file( batch.alternative_decoys_out(), tags_to_read );
		}

		//read structures and add to archive
		the_archive().read_structures( sfd, alternative_decoys_sfd, batch );

		PROF_STOP( basic::ARCHIVE_EVAL_DECOYS );
	} else {
		tr.Info << "no more decoys to read from file " << batch.silent_out() << std::endl;
		PROF_STOP( basic::ARCHIVE_READ_DECOYS );
	}
}

void ArchiveManager::unlock_file( Batch const& batch, bool final ) {
	if ( !final ) {
		tr.Debug << "...and release file" << std::endl;
		WriteOut_MpiFileBuffer file_buf( file_buf_rank_ );
		file_buf.release_file( ".//"+batch.silent_out() );
	}
}

void
ArchiveManager::jobs_completed() {// core::Size batch_id, bool final, core::Size bad ) {
	runtime_assert( jobs_completed_.begin() != jobs_completed_.end() );
	CompletionMessage msg = jobs_completed_.begin()->second;
	jobs_completed_.erase( jobs_completed_.begin() );
	Size batch_id( msg.batch_id );
	bool final( msg.final );
	Size bad( msg.bad );
	Size good_decoys( msg.good );
	Batch& batch( batches_[ batch_id ] );

	// here if in integration-test mode, jump out if not final
	if ( option[ run::constant_seed ] && !final ) return;

	tr.Debug << "jobs_completed for " << batch.batch() << "..." << "already "
		<< batch.decoys_returned() << " decoys known" << std::endl;
	runtime_assert( batch.id() == batch_id );
	WriteOut_MpiFileBuffer file_buf( file_buf_rank_ );
	if ( bad ) {
		///there were some bad-jobs --- we might be at the end of this run... hard to tell
	}
	PROF_START( basic::ARCHIVE_BLOCK_FILE );
	if ( !final ) {
		tr.Debug << "not final ... block file" << std::endl;
		//careful if file isn't in FileBuf anymore... create it just for blocking ? don't block... ?
		file_buf.block_file( ".//"+batch.silent_out() ); //destructor will release file automatically
	} else {
		tr.Debug << "final ... close file " << std::endl;
		// file_buf.close_file( ".//"+batch.silent_out() ); //that is not very nice, but file-buf isn't very smart with filenames...
		//  file_buf.close_file( ".//"+batch.score_file() ); // not required since now we have garbage-collection
	}
	if ( batch.is_cancelled() && !batch.allow_reading_cancelled_decoys() ) {
		tr.Debug << "returned decoys of cancelled batch.. ignore..." << std::endl;
		return;
	}
	PROF_STOP( basic::ARCHIVE_BLOCK_FILE );
	//sleep( 5 );
	PROF_START( basic::ARCHIVE_READ_DECOYS );


	if ( good_decoys ) {
		tr.Debug << "read file " << batch.silent_out() << std::endl;
		utility::io::izstream testin( batch.silent_out() );
		tr.Debug << "stream is " << ( testin.good() ? "good " : "bad" ) << std::endl;
		if ( !testin.good() ) { //this happens sometimes... usually it needs a little bit of waiting and then it works -- NFS lag ?
			//let's look at this later...
			jobs_completed_[ batch_id ] = msg;
			sleep( 5 );
			return;
		}

		read_returning_decoys( batch, final );

		PROF_START( basic::SAVE_ARCHIVE );
		if ( jobs_completed_.size() == 0 ) save_archive();
		PROF_STOP( basic::SAVE_ARCHIVE );
	} else { // no good decoys found
		tr.Debug << " no good decoys to read " << std::endl;
		throw CREATE_EXCEPTION(EXCN_Archive, "all decoys returned with FAIL_BAD_INPUT" );
	}

	if ( final ) {
		batch.mark_as_finished();
	}
	batch.write_info_file();
}

void
ArchiveManager::queue_batch( Batch const& batch ) {
	tr.Debug << "queue new batch into MPIArchiveJobDistributor " << batch.flag_file() << std::endl;
#ifdef USEMPI
	Size const size( 3 );
	int buf[ size ];
	buf[ 0 ] = ADD_BATCH;
	buf[ 1 ] = batch.id();
	buf[ 2 ] = batch.nstruct();
//#ifdef USEMPI
	MPI_Send( &buf, size, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//need to have MPI_JOB_DIST_TAG... since it goes into main msg-loop of JobDist

	//send size of string
	std::string strbuf( batch.flag_file() );
	buf[ 0 ] = strbuf.size();
	buf[ 1 ] = batch.id();
	MPI_Send( buf, 2, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//send string
	MPI_Send(const_cast<char*> ( strbuf.data() ), strbuf.size(), MPI_CHAR, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
#else
	protocols::jd2::JobDistributor::get_instance()->add_batch( batch.flag_file() );
#endif

}

void BaseArchiveManager::cancel_batches_previous_to( core::Size batch_id, bool allow_reading_of_decoys ) {
	for ( auto & batche : batches_ ) {
		if ( batche.id() == batch_id ) break;
		cancel_batch( batche, allow_reading_of_decoys );
	}
}

void
BaseArchiveManager::cancel_batch( Batch& batch, bool allow_reading_of_decoys ) {
	batch.mark_as_cancelled( allow_reading_of_decoys );
}

void
ArchiveManager::cancel_batch( Batch& batch, bool allow_reading_of_decoys ) {
	Parent::cancel_batch( batch, allow_reading_of_decoys );

	if ( option[ OptionKeys::run::constant_seed ]() ) {
		tr.Warning << "asked to cancel batch, but ignore in constant_seed mode to enable integration test" << std::endl;
		return;
	}
	tr.Debug << "cancel batch  " << batch.flag_file() << std::endl;
#ifdef USEMPI
	Size const size( 3 );
	int buf[ size ];
	buf[ 0 ] = CANCEL_BATCH;
	buf[ 1 ] = batch.id();
	buf[ 2 ] = batch.nstruct();
//#ifdef USEMPI
	MPI_Send( &buf, size, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//need to have MPI_JOB_DIST_TAG... since it goes into main msg-loop of JobDist

	//send size of string
	std::string strbuf( BatchJobInputter::BOGUS_BATCH_ID );
	buf[ 0 ] = strbuf.size();
	buf[ 1 ] = batch.id();
	MPI_Send( buf, 2, MPI_INT, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
	//send string
	MPI_Send(const_cast<char*> ( strbuf.data() ), strbuf.size(), MPI_CHAR, jd_master_rank_, MPI_JOB_DIST_TAG, MPI_COMM_WORLD );
#else
	protocols::jd2::JobDistributor::get_instance()->add_batch( BatchJobInputter::BOGUS_BATCH_ID );
#endif
	batch.write_info_file();
}

void ArchiveManager::send_stop_to_jobdistributor() {

	// we can also do this the quick way:
	//   a) there is still some problem and sometimes jobs-just don't finish... haven't figured out why.
	//   b) we don't really care, if we are here there was either an error,  or the Archive is converged.
	//       at the current low-acceptance rate the remaining jobs are unlikely to yield anything useful...
	//   c) quick exits saves costly time on the cluster.
	if ( option[ OptionKeys::jd2::mpi_nowait_for_remaining_jobs ]() ) {
		save_archive();
		utility_exit_with_message("quick exit from job-distributor due to flag jd2::mpi_nowait_for_remaining_jobs --- this is not an error " );

		//we do this by sending empty batch.
		tr.Debug << "send STOP signal to JobDistributor " << std::endl;
	}
	Batch stop_batch( 0 );
	queue_batch( stop_batch );
}

void
ArchiveManager::read_existing_batches() {
	using utility::file::file_exists;
	using namespace basic::options::OptionKeys;

	//possible options:
	//reads all structures from batches as if they were newly coming in
	bool b_reread_all_structures( option[ OptionKeys::archive::reread_all_structures ]() /* default false */ );

	//don't know how to probe directory... take counter...
	core::Size id( 1 );
	Batch aBatch( id );
	batches_.clear();
	while ( file_exists( aBatch.flag_file() ) ) {
		Batch& new_batch( start_new_batch() );
		runtime_assert( new_batch.id() == id );
		tr.Info << "found existing batch " << new_batch.batch() << std::endl;
		try {
			finalize_batch( new_batch, true /*reread */ );
			tr.Debug << new_batch << std::endl;
		} catch (EXCN_Archive& excn ) {
			//last started batch must have problems... ignore it
			tr.Warning << new_batch.batch()+" is errorneous: " + excn.msg() << std::endl;
			tr.Warning << "ignoring this batch..." << std::endl;
			//fill batch list if exception state has left us without it
			if ( batches_.size() < id ) batches_.push_back( Batch( id ) );
			batches_.back().mark_as_invalid();
		}
		if ( b_reread_all_structures ) {
			if ( batches_[ id ].decoys_returned() ) {
				jobs_completed_[ id ] =
					CompletionMessage( id, batches_[id ].has_finished(), 0, batches_[ id ].decoys_returned(), batches_[ id ].nstruct() );
			}
			batches_[ id ].set_decoys_returned( 0 );
		}
		aBatch = Batch( ++id );
	}
}

Batch&
BaseArchiveManager::start_new_batch() {
	using utility::file::file_exists;

	core::Size batch_id( last_batch_id() + 1 );
	tr.Debug << "start new batch " << batch_id << std::endl;
	batches_.push_back( Batch( batch_id ) );
	Batch &new_batch( batches_.back() );

	new_batch.set_id( batch_id );
	//make directory:
	utility::file::create_directory( new_batch.dir() );
	new_batch.user_options().add_built_in_options();
	add_all_rosetta_options( new_batch.user_options() );

	//copy the system broker setup --- OBSOLET since Sept 20th 2010. now broker:setup is FileVector option.
	//  if ( !file_exists( new_batch.broker_file() ) && option[ OptionKeys::broker::setup ].user() ) {
	//   utility::io::ozstream batch_broker( new_batch.broker_file() );
	//   utility::io::izstream system_broker( option[ OptionKeys::broker::setup ]() );
	//   std::string line;
	//   while ( getline( system_broker, line ) ) batch_broker << line << std::endl;
	//  }
	new_batch.nstruct() = basic::options::option[ basic::options::OptionKeys::out::nstruct ];
	return batches_.back();
}

void report_batch_inconsistency( Batch& new_batch, std::string const &tag ) {
	throw CREATE_EXCEPTION(EXCN_Archive, "inconsistency detected when re-reading "+new_batch.batch()+" for " + tag);
}

void
BaseArchiveManager::finalize_batch( Batch& new_batch, bool reread ) {
	using utility::file::file_exists;
	using namespace basic::options::OptionKeys;
	tr.Debug << "finalize_batch " << new_batch << std::endl;
	if ( !reread ) new_batch.set_decoys_returned( 0 );

	if ( !utility::file::file_exists( new_batch.broker_file() ) ) {
		utility::io::ozstream broker( new_batch.broker_file() );
		broker << "# NO CLAIMERS PRESENT" << std::endl;
		broker.close();
	}

	if ( file_exists( new_batch.flag_file() ) ) {
		tr.Debug << "checking aBatch.flag_file()... " << std::endl;
		utility::options::OptionCollection batch_opts;
		batch_opts.add_built_in_options();
		add_all_rosetta_options( batch_opts );
		try {
			tr.Debug << "load options from file" << std::endl;
			batch_opts.load_options_from_file_exception( new_batch.flag_file() );
		} catch ( utility::excn::Exception &excn ) {
			tr.Error << "problems with flags in " << new_batch.flag_file() << " aborting... " << std::endl;
			// excn.show( tr.Error );
			batches_.pop_back();
			throw CREATE_EXCEPTION(EXCN_Archive, new_batch.flag_file() + " contains errors: " + excn.msg() );
		}
		if ( !reread )  {
			//access all archive controlled options... so they are not in the "user_flags" anymore
			if ( batch_opts[ in::file::silent ].user() ) {
				tr.Warning << "option -in:file:silent will be overwritten by ArchiveMaster"
					<< " -- control directly via class Batch" << std::endl;
			}
			if ( batch_opts[ out::nstruct ].user() ) {
				tr.Warning << "option -nstruct will be overwritten by ArchiveMaster "
					<< "-- control directly via class Batch" << std::endl;
			}
			if ( batch_opts[ run::intermediate_structures ].user() ) {
				tr.Warning << "option -run::intermediate_structures will be overwritten by ArchiveMaster "
					<< "-- control directly via class Batch" << std::endl;
			}
			if ( batch_opts[ out::file::silent ].user() ) {
				tr.Warning << "option -out:file:silent will be overwritten by ArchiveMaster "
					<< "-- control directly via class Batch" << std::endl;
			}
			if ( batch_opts[ broker::setup ].user() ) {
				tr.Warning << "option -broker:setup will be overwritten by ArchiveMaster "
					<< "-- control directly via class Batch" << std::endl;
			}
			if ( batch_opts[ out::file::scorefile ].user() ) {
				tr.Warning << "option -out:file:scorefile will be overwritten by ArchiveMaster "
					<< "-- control directly via class Batch" << std::endl;
			}
		}

		core::Size nstruct( batch_opts[ out::nstruct ]() );
		//bool intermeds( batch_opts[ run::intermediate_structures ]() );
		std::string silent_out( batch_opts[ out::file::silent ]() );
		utility::vector1< std::string > broker( batch_opts[ broker::setup ]() );
		std::ostringstream broker_files;
		std::copy( broker.begin(), broker.end(), std::ostream_iterator<std::string>( broker_files, " "));
		//std::string score_file( batch_opts[ out::file::scorefile ]() );

		// now the other options are "inaccessed options" and can be dumped to a stream
		std::stringstream user_flags;
		batch_opts.show_inaccessed_user_options( user_flags );
		tr.Debug << "user_options: \n" << user_flags.str() << std::endl;

		// and can be added to the batch-options
		new_batch.user_options().load_options_from_stream( user_flags, "USER_FLAGS" );
		if ( reread ) {
			bool has_silent( batch_opts[ in::file::silent ].user() );
			new_batch.read_info_file();
			//new_batch.set_intermediate_structs( intermeds ); //this is not read from BATCH_INFO

			// for all other values we just double-check consistency
			if ( new_batch.nstruct() != nstruct ) report_batch_inconsistency( new_batch, "NSTRUCT" );
			if ( new_batch.has_silent_in() !=  has_silent ) report_batch_inconsistency( new_batch, "INPUT" );
			if ( silent_out != new_batch.silent_out() ) report_batch_inconsistency( new_batch, "OUTPUT" );
			if ( broker_files.str() != new_batch.all_broker_files() ) report_batch_inconsistency( new_batch, "BROKER_FILE" );
			//TODO: determine how many decoys have been returned to archive...
		}
	}

	//now write the final flag-file
	utility::io::ozstream flag_out( new_batch.flag_file() );
	new_batch.user_options().show_user( flag_out );
	flag_out << "\n\n#Archive controlled flags" << std::endl;
	flag_out << "-out:file:silent " << new_batch.silent_out() << std::endl;
	if ( new_batch.has_silent_in() ) flag_out << "-in:file:silent " << new_batch.silent_in() << std::endl;

	flag_out << "-out:nstruct " << new_batch.nstruct() << std::endl;
	flag_out << "-out:file:scorefile " << new_batch.score_file() << std::endl;
	flag_out << "-broker:setup " << new_batch.all_broker_files() << std::endl;

	if ( new_batch.intermediate_structs() ) flag_out << "-run:intermediate_structures" << std::endl;

	if ( !reread ) {
		new_batch.write_info_file();
	}

	if ( !new_batch.has_finished() && !new_batch.is_cancelled() && the_archive().still_interested( new_batch ) ) {
		tr.Debug << "queue " << new_batch.batch() << " " << new_batch.flag_file() << std::endl;
		queue_batch( new_batch );
	} else {
		new_batch.mark_as_finished();
	}

	tr.Debug << "\n" << std::endl;

}


void
ArchiveManager::save_archive() {
	the_archive().save_to_file();
}


bool
ArchiveManager::restore_archive() {
	return the_archive().restore_from_file();
}

//#endif //ndef WIN32

}//archive
}//jd2
}//protoco
