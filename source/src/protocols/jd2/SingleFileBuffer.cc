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
/// @detail this outputter will send silentstructs via MPI to dedicated node that will collect all structures
/// @author Oliver Lange olange@u.washington.edu

#include <protocols/jd2/SingleFileBuffer.hh>
//#include <protocols/jd2/MpiFileBuffer.hh> //only needed for the runtime_assert grumpf.
#include <iostream>
#include <utility/io/mpistream.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>
#include <iterator>

#include <utility/vector1.hh>

#ifdef WIN32
#include <windows.h>  // for sleep()
#include <ctime>      // for clock()
#else
#include <unistd.h>
#endif

namespace protocols {
namespace jd2 {

using namespace core;
using namespace utility::io::mpi_stream;

static thread_local basic::Tracer tr( "protocols.jd2.MpiFileBuffer" );

/// @details this is a implementation of Buffer for silent-file-based output.
std::string const START_BLOCK( "MPI_FILE_BUFFER_BLOCK_START" );
std::string const END_BLOCK( "MPI_FILE_BUFFER_BLOCK_END" );

SingleFileBuffer::~SingleFileBuffer() { runtime_assert( !has_open_slaves() ); }

void SingleFileBuffer::flush( Size slave ) {
	// std::cout << "flush channel: " << filename() << " for slave : " << slave << std::endl;
	if ( unfinished_blocks_[ slave ].size() ) {
		write_lines( unfinished_blocks_[ slave ] );
		unfinished_blocks_[ slave ].clear();
	}
}

bool SingleFileBuffer::has_open_slaves() const {
	for ( BufferMap::const_iterator it = unfinished_blocks_.begin(); it != unfinished_blocks_.end(); ++it ) {
		if ( it->second.size() > 0 ) {
			return true;
		}
	}
	return false;
}

Size SingleFileBuffer::nr_open_slaves() const {
	return unfinished_blocks_.size();
}

void SingleFileBuffer::close( Size slave ) {
	if ( unfinished_blocks_[ slave ].size() == 0 ) {
		tr.Trace << "EMPTY OPEN/CLOSE Operation from slave " << slave << std::endl;
	}
	tr.Debug << "close " << filename_ << " from slave " << slave << std::endl;
	flush( slave );
	//BufferMap::iterator iter = unfinished_blocks_.find( slave );
	// if ( iter!=unfinished_blocks_.end() ) {
	// unfinished_blocks_.erase( iter );
	//} else {
	// tr.Warning << "tried to close non-existant channel to slave-node " << slave << " for file " << filename_ << std::endl;
	//}
}

void SingleFileBuffer::store_line( Size slave, Size channel, std::string const& line ) {
	runtime_assert( channel == mpi_channel_ );
	unfinished_blocks_[ slave ].push_back( line );
	// std::cout << "channel: " << mpi_channel_ << " slave: " << slave <<std::endl;// << "line: " << line << std::endl;
}

core::Size SingleFileBuffer::length( core::Size slave ) {

	core::Size length = 0;
	LineBuffer & buf( unfinished_blocks_[ slave ] );
	for ( LineBuffer::iterator iter = buf.begin(); iter != buf.end(); ++iter ) {
		length += iter->length();
	}

	return length;
}

void SingleFileBuffer::write_lines( LineBuffer const& buf ) {
	std::cout << "START_BLOCK" << std::endl;
	copy( buf.begin(), buf.end(), std::ostream_iterator< std::string>( std::cout ) );
	std::cout << "END_BLOCK" << std::endl;
}


void SingleFileBuffer::block( core::Size MPI_ONLY( slave ) ) {
#ifdef USEMPI
	int status = 1;
	tr.Debug << "send blocking confirmation... " << filename() << std::endl;
	MPI_Send( &status, 1, MPI_INT,  slave, MPI_STREAM_TAG, MPI_COMM_WORLD );
	tr.Debug << "blocked..." << std::endl;
	int buf[ 4 ];
	MPI_Status stat;
	MPI_Recv( &buf, 4, MPI_INT, slave, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );
	tr.Debug << "release file? : received: " << buf[ 0 ] << " " << buf[ 1 ] << " " << buf[ 2 ] << " " << buf[ 3 ] << std::endl;
	//	runtime_assert( (Size) buf[ 2 ] == MPI_RELEASE_FILE && (Size) buf[ 1 ] == filename().size() ); //check of sizes ..
#endif
}

WriteFileSFB::WriteFileSFB( std::string const& filename, core::Size channel, bool append, core::Size& status ) :
	SingleFileBuffer( filename, channel, status ) {//, out_( filename.c_str() )
	//if ( append ) out_.open_append(filename); //still problems with this ???
	status = MPI_FAIL;
	int trials = 5;
	while ( status == MPI_FAIL && trials > 0 ) {
		--trials;
		if ( append ) {
			if ( !utility::file::file_exists( filename ) ) {
				//   out_.open( filename.c_str() );
				tr.Debug << "open file " << filename << " ... " << std::endl;
				out_.open( filename.c_str() );
				if ( out_.good() ) status = MPI_SUCCESS_NEW;
			} else {

				tr.Debug << "open file (append) " << filename << " ... " << std::endl;
				//   out_.open_append( filename.c_str() );
				out_.open( filename.c_str(), std::ios::app );
				if ( out_.good() ) status = MPI_SUCCESS_APPEND;
			}
		} else {
			//  out_.open( filename.c_str() );
			tr.Debug << "open file " << filename << " ... " << std::endl;
			out_.open( filename.c_str() );
			if ( out_.good() ) status = MPI_SUCCESS_NEW;
		}
		if ( status == MPI_FAIL ) {
			std::cerr << "failing to write file " << filename << " try again after 1 second of sleep ... " << std::endl;

#ifdef WIN32
			Sleep(1000);
#else
			sleep( 1 );
#endif
		}
	}
}

WriteFileSFB::~WriteFileSFB() {
	tr.Debug << "close file " << filename() << std::endl;
	out_.close();
}

void WriteFileSFB::write_lines( LineBuffer const& buf ) {
	static time_t const all_start_time = clock();
	static time_t ntime;
	static time_t last_time;
	ntime = clock();
	if ( buf.size()==1 ) tr.Debug << -1.0*(all_start_time-ntime)/CLOCKS_PER_SEC << " " << 1.0*( ntime-last_time )/CLOCKS_PER_SEC << " seconds: write block of " << buf.begin()->size() << " characters to file " << filename() << std::endl;
	else tr.Debug << -1.0*(all_start_time-ntime)/CLOCKS_PER_SEC << " " <<1.0*( ntime-last_time )/CLOCKS_PER_SEC << " seconds: write " << buf.size() << " blocks of data to file " << filename() << std::endl;
	last_time = ntime;

	copy( buf.begin(), buf.end(), std::ostream_iterator< std::string>( out_ ) );
	if ( tr.Trace.visible() ) {
		copy( buf.begin(), buf.end(), std::ostream_iterator< std::string>( tr.Trace ) );
		tr.Trace << std::endl;
	}
}

void WriteFileSFB::block( core::Size slave ) {
	out_.close();
	tr.Debug << "block file " << filename() << std::endl;
	//out_.flush();
	Base::block(slave);
	tr.Debug << "open file (append): " << filename() << std::endl;
	out_.open( filename().c_str() , std::ios::app );
}

}
}
