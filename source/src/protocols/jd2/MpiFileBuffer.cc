// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MpiFileBuffer.hh
/// @brief  header file for MPISilentFileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @details this outputter will send silentstructs via MPI to dedicated node that will collect all structures
/// @author Oliver Lange olange@u.washington.edu

#include <protocols/jd2/MpiFileBuffer.hh>
#include <iostream>
#include <utility/io/mpistream.hh>
#include <utility/exit.hh>
#ifdef USEMPI
#include <utility/file/file_sys_util.hh>
#endif

#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <ctime>

#include <protocols/jd2/SingleFileBuffer.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

using namespace core;
using namespace utility::io::mpi_stream;

static thread_local basic::Tracer tr( "protocols.jd2.MpiFileBuffer" );
using basic::mem_tr;

MpiFileBuffer::MpiFileBuffer( Size file_buf_rank )
: buffer_rank_( file_buf_rank ),
	last_channel_( 0 ),
	bSlaveCanOpenFile_( true ),
	bKeepFilesAlive_( true ),
	seconds_to_keep_files_alive_( 100 ) {

	if ( seconds_to_keep_files_alive_<=0 ) bKeepFilesAlive_ = false;
	last_garbage_collection_ = time(NULL);
#ifdef USEMPI
	int my_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank );/* get current process id */
	my_rank_ = my_rank;
#endif

}

#ifdef USEMPI
Size const MPI_WIND_DOWN( 1001 );
Size const MPI_BLOCK_FILE( 1002 );
Size const MPI_RELEASE_FILE( 1003 );
#endif
Size const MPI_CLOSE_FILE( 1004 );

#ifdef USEMPI
void MpiFileBuffer::receive_str( Size slave, Size size, std::string& line ) {
	char *cbuf = new char[ size+1 ];
	MPI_Status stat;
  MPI_Recv( cbuf, size, MPI_CHAR, slave, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );
  line.assign( cbuf, size );
	delete[] cbuf;
}
#else
void MpiFileBuffer::receive_str( Size , Size , std::string&  ) {}
#endif

MpiFileBuffer::~MpiFileBuffer() {
	while ( blocked_files_.size() ) {
		release_file( blocked_files_.front() );
	}
	garbage_collection();
}

void MpiFileBuffer::show_status( std::ostream& os ) const {
	os << "MpiFileBuffer Status Report...." << std::endl;
	os << "known filenames: " << open_files_.size() << std::endl;
	os << "open buffers:    " << open_buffers_.size() << std::endl;
	os << "garbage_list:    " << garbage_collector_.size() << std::endl;
	os << "blocked_files:   " << blocked_files_.size() << std::endl;
}

void MpiFileBuffer::garbage_collection() {
	time_t const now( time(NULL) );
	//tr.Debug << "last garbage collection " << now-last_garbage_collection_ << " seconds ago" << std::endl;
	if ( now-last_garbage_collection_ < 30 ) return;
	if ( tr.Debug.visible() ) show_status( tr.Debug );
	tr.Debug << "garbage collection active..." << std::endl;
	for ( GarbageList::iterator it=garbage_collector_.begin(); it!=garbage_collector_.end(); ) {
		GarbageList::iterator to_erase( it );
		++it; //move forward now, because to_erase will be erased in "close_file"
		tr.Debug << "marked " << to_erase->first << " " << now-to_erase->second << " seconds ago." << std::endl;
		if ( now-to_erase->second > seconds_to_keep_files_alive_ ) {
			int channel( to_erase->first );
			SingleFileBufferOP buf = open_buffers_[ channel ];
			if ( !buf ) {
				garbage_collector_.erase( to_erase );
				// this is now understood. some files can be closed directly via MPI_CLOSE_FILE (e.g., from ArchiveManager )
				// to avoid this case I could call clear_channel_from_garbage_collector( channel_id ) in close_file()
				// but then one has to be careful to not delete the iterator in this for-loop
				continue;
			}
			if ( !buf->has_open_slaves() ) {
				tr.Debug << "channel "<< channel
					<< " has no more open slaves... and has not been touched again --- close via garbage collector" << std::endl;
				//    garbage_collector_.erase( to_erase ); now removed from garbage_collector in close_file()
				close_file( channel );
				mem_tr << "closed_channel" << std::endl;
			} else {
				runtime_assert( false ); //shouldn't happen anymore
				tr.Debug << "channel " << to_erase->first << " has open slaves again ... not closed, remove from closing list" << std::endl;
				//    garbage_collector_.erase( to_erase );
			}
		}
	}
	last_garbage_collection_ = now;
}

void MpiFileBuffer::run() {
#ifdef USEMPI
	tr.Debug << "MpiFileBuffer " << (( my_rank_ != buffer_rank_ ) ? "not" : "" ) <<  " started on node " << my_rank_ << std::endl;
  if ( my_rank_ != buffer_rank_ ) return;
	MPI_Status stat;
	bStop_ = false; //this might be changed from some where via msg.
  while( !bStop_ ) {
		int buf[ 4 ];
		MPI_Recv( buf, 4, MPI_INT, MPI_ANY_SOURCE, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );
		Size const msg_type( buf[ 2 ] );
		Size const size( buf[ 1 ] );
		Size const slave( buf[ 0 ] );
		Size const channel_id( buf[ 3 ] );
		if ( msg_type == MPI_STREAM_OPEN || msg_type == MPI_STREAM_OPEN_APPEND ) {
			std::string filename;
			receive_str( slave, size, filename );
			Size file_status;
			open_channel( slave, filename, msg_type == MPI_STREAM_OPEN_APPEND, file_status );
		} else if ( msg_type == MPI_STREAM_SEND ) {
			std::string line;
			receive_str( slave, size, line );
			store_to_channel( slave, channel_id, line );
		} else if ( msg_type == MPI_STREAM_CLOSE ) {
			close_channel( slave, channel_id );
		} else if ( msg_type == MPI_STREAM_FLUSH ) {
			tr.Debug << "MPI_STREAM_FLUSH received" << std::endl;
			flush_channel( slave, channel_id );
		} else if ( msg_type == MPI_WIND_DOWN ) {
			bStop_ = true;
		} else if ( msg_type == MPI_BLOCK_FILE ) {
			std::string filename;
			receive_str( slave, size, filename );
			block_file( slave, filename );
		} else if ( msg_type == MPI_STREAM_FILE_EXIST ) {
			std::string filename;
			receive_str( slave, size, filename );
			bool exist( is_open_channel( filename ) );
			if ( !exist ) exist = utility::file::file_exists( filename );
			int iexist( exist ? 1 : 0 );
			MPI_Send( &iexist, 1, MPI_INT, slave, MPI_STREAM_TAG, MPI_COMM_WORLD );
		} else if ( msg_type == MPI_CLOSE_FILE ) {
			std::string filename;
			receive_str( slave, size, filename );
			tr.Debug << "received MPI_CLOSE_FILE " <<filename<<" from Node" << slave << std::endl;
			bool success = close_file( filename );
			//send confirmation...
			int closed = success ? 1 : 0;
			MPI_Send( &closed, 1, MPI_INT, slave, MPI_STREAM_TAG, MPI_COMM_WORLD );
		} else {
			utility_exit_with_message( "unknown msg-id received in MpiFileBuffer.cc");
		}
		garbage_collection();
  }
#endif
}


void MpiFileBuffer::stop() {
#ifdef USEMPI
  if (my_rank_ == buffer_rank_ ) { bStop_=true; return; };
	int buf[ 4 ];
  buf[ 0 ] = my_rank_;
  buf[ 1 ] = 0;
	buf[ 2 ] = MPI_WIND_DOWN;
	buf[ 3 ] = 0;
	MPI_Send( &buf, 4, MPI_INT, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
#endif
}

void MpiFileBuffer::block_file( Size from_node, std::string const& filename ) {
	runtime_assert( my_rank_ == buffer_rank_ );
	Filenames::const_iterator iter = open_files_.find( filename );
	Size channel;

	if ( iter != open_files_.end() ) {
		channel = iter->second;
		SingleFileBufferOP buf = open_buffers_[ channel ];
		runtime_assert( buf != 0 ); //consistent?
		buf->block( from_node ); //closes-file, sends MPI signal back and forth and hangs until release, reopens file
		tr.Debug << "block released... for file " << filename << std::endl;
	} else {
		tr.Warning << "file " << filename << " is not known to MpiFileBuffer " << std::endl;
#ifdef USEMPI
		int status = 0;
		tr.Debug << "send blocking confirmation... " << filename << std::endl;
		MPI_Send( &status, 1, MPI_INT,  from_node, MPI_STREAM_TAG, MPI_COMM_WORLD );
#endif
	}
}

void MpiFileBuffer::block_file( std::string const& MPI_ONLY(filename) ) {
#ifdef USEMPI
	runtime_assert( buffer_rank_ != my_rank_ );
	int buf[ 4 ];
  buf[ 0 ] = my_rank_;
  buf[ 1 ] = filename.size();
	buf[ 2 ] = MPI_BLOCK_FILE;
	buf[ 3 ] = 0;

  MPI_Send(buf, 4, MPI_INT, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
	MPI_Send(const_cast<char*> (filename.data()), filename.size(), MPI_CHAR, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
	tr.Debug << "wait for confirmation of block  " << std::endl;
	MPI_Status stat;
	MPI_Recv( &buf, 1, MPI_INT, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );//file is blocked when we receive this
	if ( buf[ 0 ] == 1 ) {
		tr.Debug << "block confirmed... " << std::endl;
		blocked_files_.push_back( filename );
	} else {
		tr.Debug << "block not accepted ... " << std::endl;
	}
	// now do stuff and then say release_file
#endif
}


bool MpiFileBuffer::remote_close_file( std::string const& filename ) {
	runtime_assert( buffer_rank_ != my_rank_ );
	int buf[ 4 ];
	buf[ 0 ] = my_rank_;
	buf[ 1 ] = filename.size();
	buf[ 2 ] = MPI_CLOSE_FILE;
	buf[ 3 ] = 0;
#ifdef USEMPI

  MPI_Send(buf, 4, MPI_INT, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
	MPI_Send(const_cast<char*> (filename.data()), filename.size(), MPI_CHAR, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );

	MPI_Status stat;
	MPI_Recv( &buf, 1, MPI_INT, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );//file is blocked when we receive this
#endif
	if ( buf[ 0 ] == 1 ) {
		tr.Debug << "remote close confirmed... " << std::endl;
	} else {
		tr.Debug << "close not accepted ... " << std::endl;
	}
	return buf[ 0 ] != 0;
	// now do stuff and then say release_file

}


void MpiFileBuffer::release_file( std::string filename ) {
	//if this is a reference I get seqfault, since it might point into list where I erase from...
	std::list< std::string >::iterator iter = find( blocked_files_.begin(), blocked_files_.end(), filename );
	if ( iter != blocked_files_.end() ) {
		blocked_files_.erase( iter );
		runtime_assert( buffer_rank_ != my_rank_ );
#ifdef USEMPI
		int buf[ 4 ];
		buf[ 0 ] = my_rank_;
		buf[ 1 ] = filename.size();
		buf[ 2 ] = MPI_RELEASE_FILE;
		buf[ 3 ] = 0;
//#ifdef USEMPI
		tr.Debug << "release file " << filename << std::endl;
		MPI_Send(buf, 4, MPI_INT, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
		//don't send filename again....
		//MPI_Send(const_cast<char*> (filename.data()), filename.size(), MPI_CHAR, buffer_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
#endif
	} else {
		tr.Debug << filename << " is not blocked... release ignored " << std::endl;
	}
}


bool MpiFileBuffer::is_open_channel( std::string const& filename ) {
	Filenames::const_iterator iter = open_files_.find( filename );
	//Size channel;
	return ( iter != open_files_.end() );
}

void MpiFileBuffer::clear_channel_from_garbage_collector( core::Size channel ) {
	//remove from garbage_collector if its in there already
	GarbageList::iterator giter = garbage_collector_.find( channel );
	if ( giter != garbage_collector_.end() ) {
		tr.Debug << "remove channel " << channel << " from garbage-collector list " << std::endl;
		garbage_collector_.erase( giter );
	}
}

#ifdef USEMPI
void MpiFileBuffer::open_channel( Size slave, std::string const& filename, bool append, Size& status ) {
	//find filename in our file-map...
	tr.Debug << "open mpi-channel from slave-node " << slave << " for file: " << filename << std::endl;
	Filenames::const_iterator iter = open_files_.find( filename );
	Size channel;
	if ( iter != open_files_.end() ) {
		channel = iter->second;
		SingleFileBufferOP buf = open_buffers_[ channel ];
		runtime_assert( buf != 0 ); //consistent?
		//buf->flush( slave );
		tr.Debug << "channel exists already: " << channel << std::endl;
		status = MPI_SUCCESS_APPEND;

		clear_channel_from_garbage_collector( channel );
	} else {
		//new file
		channel = ++last_channel_; //this might overrun eventually
		runtime_assert( channel < 2147483647 /* 2^32 -1 */ );
		open_files_.insert( Filenames::value_type( filename, channel ) );
		open_buffers_.insert( Buffers::value_type( channel, generate_new_channel( filename, channel, append, status )) );
		tr.Debug << "new channel established: " << channel << std::endl;
	}
  int send_buf[ 2 ];
	send_buf[ 0 ] = channel;
	send_buf[ 1 ] = status;
	MPI_Send( &send_buf, 2, MPI_INT, slave, MPI_STREAM_TAG, MPI_COMM_WORLD );
	if ( status == MPI_SUCCESS_NEW ) {
			int buf[ 4 ];
			MPI_Status stat;
			MPI_Recv( buf, 4, MPI_INT, slave, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );
			Size const slave_id( buf[ 0 ] );
			Size const size( buf[ 1 ] );
			Size const msg_type( buf[ 2 ] );
			Size const channel_id( buf[ 3 ] );
			tr.Debug << "header? : received: " << slave_id << " " << size << " " << msg_type << " " << channel_id << std::endl;
			runtime_assert( msg_type == MPI_STREAM_SEND && slave_id == slave );
			std::string header;
			receive_str( slave, size, header );
			tr.Debug << header << std::endl;
			open_buffers_[ channel ]->store_line( slave, channel, header );
			open_buffers_[ channel ]->flush( slave ); //flush to write header first
	}
	//send file-descriptor out
}
#else
void MpiFileBuffer::open_channel( Size , std::string const&, bool, Size& status ) {
	status = 0;
	(void) last_channel_;
}
#endif


void MpiFileBuffer::store_to_channel( Size slave, Size channel, std::string const& line ) {
	//tr.Debug << "store channel for slave " << slave << " channel: " << channel << " length: " << line.length() << std::endl;
	SingleFileBufferOP buf = open_buffers_[ channel ];
	runtime_assert( buf != 0 ); //writing to open file ?
	buf->store_line( slave, channel, line );
	if ( buf->length(slave) > 5e6 ) {
		tr.Info << "autoflush threshold (5 MB) triggered for slave " << slave << " channel: " << channel << std::endl;
		flush_channel(slave, channel);
	}
}

void MpiFileBuffer::flush_channel( Size slave, Size channel_id ) {
	tr.Debug << "flush channel for slave " << slave << " channel: " << channel_id << std::endl;
	SingleFileBufferOP buf = open_buffers_[ channel_id ];
	runtime_assert( buf != 0 ); //writing to open file ?
	buf->flush( slave );
}

void MpiFileBuffer::close_channel( Size slave, Size channel_id ) {
	SingleFileBufferOP buf = open_buffers_[ channel_id ];
	runtime_assert( buf != 0 ); //writing to open file ?
	buf->close( slave );
	tr.Debug << "close channel "<< channel_id <<" for slave " << slave
		<< " currently " << buf->nr_open_slaves() << " open slave buffers; open files: " << open_buffers_.size() << std::endl;
	if ( !buf->has_open_slaves() && bSlaveCanOpenFile_ && !bKeepFilesAlive_ ) {
		tr.Debug << "channel has no more open slaves... close completely now" << std::endl;
		close_file( channel_id );
		mem_tr << "closed_channel" << std::endl;
	} else {
		if ( !buf->has_open_slaves() ) {
			tr.Debug << "mark channel " << channel_id << " for closing in " << seconds_to_keep_files_alive_ << std::endl;
			garbage_collector_[ channel_id ]=time(NULL);
		}
	}
	tr.Debug << "currently " << open_buffers_.size() << " open buffers" << std::endl;
}

bool MpiFileBuffer::close_file( std::string filename ) {
	if ( my_rank_ != buffer_rank_ ) {
		tr.Debug << "remote close file " << filename << " on node " << my_rank_ << std::endl;
		return remote_close_file( filename );
	} else {
		if ( is_open_channel( filename ) ) {
			Size channel_id = open_files_[ filename ];
			tr.Debug << "close file " << filename << " with channel_id " << channel_id << std::endl;
			close_file( channel_id );
			return true;
		}
	}
	return false;
}

void MpiFileBuffer::close_file( Size channel_id ) {
	Buffers::iterator iter = open_buffers_.find( channel_id );
	if ( iter != open_buffers_.end() ) {
		std::string filename( iter->second->filename() );
		Filenames::iterator file_iter = open_files_.find( filename );
		runtime_assert( file_iter != open_files_.end() ); //if this is not true we have inconsistency between open_buffers and open_files_
		open_files_.erase( file_iter );
		open_buffers_.erase( iter );
	} else {
		utility_exit_with_message( "illegal attempt to delete file with non-existent channel_id " );
	}
}

SingleFileBufferOP WriteOut_MpiFileBuffer::generate_new_channel( std::string const& filename, Size channel, bool append, Size& status ) {
	return SingleFileBufferOP( new WriteFileSFB( filename, channel, append, status ) );
}

SingleFileBufferOP DebugOut_MpiFileBuffer::generate_new_channel( std::string const& filename, Size channel, bool /*append*/, Size& status ) {
	return SingleFileBufferOP( new SingleFileBuffer( filename, channel, status ) );
}

}
}
