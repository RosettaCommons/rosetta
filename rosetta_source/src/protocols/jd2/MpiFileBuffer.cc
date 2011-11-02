// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MpiFileBuffer.hh
/// @brief  header file for MPISilentFileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @detail this outputter will send silentstructs via MPI to dedicated node that will collect all structures
/// @author Oliver Lange olange@u.washington.edu

#include <protocols/jd2/MpiFileBuffer.hh>
#include <iostream>
#include <utility/io/mpistream.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
// AUTO-REMOVED #include <iterator>

#include <protocols/jd2/SingleFileBuffer.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

using namespace core;
using namespace utility::io::mpi_stream;

static basic::Tracer tr("protocols.jd2.MpiFileBuffer");
using basic::mem_tr;

MpiFileBuffer::MpiFileBuffer( Size file_buf_rank )
	: buffer_rank_( file_buf_rank ), last_channel_( 0 ), bSlaveCanOpenFile_( true ),bKeepFilesAlive_( true )  {
#ifdef USEMPI
	int my_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank );/* get current process id */
	my_rank_ = my_rank;
#endif
}

Size const MPI_WIND_DOWN( 1001 );
Size const MPI_BLOCK_FILE( 1002 );
Size const MPI_RELEASE_FILE( 1003 );
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
}

void MpiFileBuffer::run() {
#ifdef USEMPI
	tr.Debug << "MpiFileBuffer " << (( my_rank_ != buffer_rank_ ) ? "not" : "" ) <<  " started on node " << my_rank_ << std::endl;
  if (my_rank_ != buffer_rank_ ) return;
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
		runtime_assert( buf ); //consistent?
		buf->block( from_node );
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
//if this is a reference i get seqfault, since it might point into list where I erease  from...
	std::list< std::string >::iterator iter = find( blocked_files_.begin(), blocked_files_.end(), filename );
	if ( iter != blocked_files_.end() ) {
		blocked_files_.erase( iter );
		runtime_assert( buffer_rank_ != my_rank_ );
		int buf[ 4 ];
		buf[ 0 ] = my_rank_;
		buf[ 1 ] = filename.size();
		buf[ 2 ] = MPI_RELEASE_FILE;
		buf[ 3 ] = 0;
#ifdef USEMPI
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


#ifdef USEMPI
void MpiFileBuffer::open_channel( Size slave, std::string const& filename, bool append, Size& status ) {
	//find filename in our file-map...
	tr.Debug << "open mpi-channel from slave-node " << slave << " for file: " << filename << std::endl;
	Filenames::const_iterator iter = open_files_.find( filename );
	Size channel;
	if ( iter != open_files_.end() ) {
		channel = iter->second;
		SingleFileBufferOP buf = open_buffers_[ channel ];
		runtime_assert( buf ); //consistent?
		buf->flush( slave );
		tr.Debug << "channel exists already: " << channel << std::endl;
		status = MPI_SUCCESS_APPEND;
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
			Size const msg_type( buf[ 2 ] );
			Size const size( buf[ 1 ] );
			Size const slave_id( buf[ 0 ] );
			Size const channel_id( buf[ 3 ] );
			tr.Debug << "header? : received: " << buf[ 0 ] << " " << buf[ 1 ] << " " << buf[ 2 ] << " " << buf[ 3 ] << std::endl;
			runtime_assert( msg_type == MPI_STREAM_SEND && slave_id == slave );
			std::string header;
			receive_str( slave, size, header );
			tr.Debug << header << std::endl;
			open_buffers_[ channel ]->store_line( slave, channel, header );
			open_buffers_[ channel ]->flush( slave );
	}
	//send file-descriptor out
}
#else
void MpiFileBuffer::open_channel( Size , std::string const&, bool, Size& status ) {
	status = 0;
}
#endif


void MpiFileBuffer::store_to_channel( Size slave, Size channel, std::string const& line ) {
	SingleFileBufferOP buf = open_buffers_[ channel ];
	runtime_assert( buf ); //writing to open file ?
	buf->store_line( slave, channel, line );
}

void MpiFileBuffer::flush_channel( Size slave, Size channel_id ) {
	tr.Debug << "flush channel for slave " << slave << " channel: " << channel_id << std::endl;
	SingleFileBufferOP buf = open_buffers_[ channel_id ];
	runtime_assert( buf ); //writing to open file ?
	buf->flush( slave );
}

void MpiFileBuffer::close_channel( Size slave, Size channel_id ) {
	tr.Debug << "close channel for slave " << slave << " channel: " << channel_id << std::endl;
	SingleFileBufferOP buf = open_buffers_[ channel_id ];
	runtime_assert( buf ); //writing to open file ?
	buf->close( slave );
	if ( !buf->has_open_slaves() && bSlaveCanOpenFile_ && !bKeepFilesAlive_ ) {
		tr.Debug << "channel has no more open slaves... close completely now" << std::endl;
		close_file( channel_id );
		mem_tr << "closed_channel" << std::endl;
	}
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
		utility_exit_with_message( "illegal attempt to delete file with non-existant channel_id " + channel_id );
	}
}

SingleFileBufferOP WriteOut_MpiFileBuffer::generate_new_channel( std::string const& filename, Size channel, bool append, Size& status ) {
	return new WriteFileSFB( filename, channel, append, status );
}

SingleFileBufferOP DebugOut_MpiFileBuffer::generate_new_channel( std::string const& filename, Size channel, bool /*append*/, Size& status ) {
	return new SingleFileBuffer( filename, channel, status );
}

}
}
