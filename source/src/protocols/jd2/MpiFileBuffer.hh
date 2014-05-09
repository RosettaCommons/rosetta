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


#ifndef INCLUDED_protocols_jd2_MpiFileBuffer_hh
#define INCLUDED_protocols_jd2_MpiFileBuffer_hh

#ifdef USEMPI
#include <mpi.h>
#endif

// unit headers
#include <protocols/jd2/MpiFileBuffer.fwd.hh>
// AUTO-REMOVED #include <protocols/jd2/SingleFileBuffer.hh>

//project headers
#include <core/types.hh>

//utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
//#include <utility/io/mpistream.hh>

//C++ headers

// AUTO-REMOVED #include <string>
#include <map>
#include <list>

#include <protocols/jd2/SingleFileBuffer.fwd.hh>


namespace protocols {
namespace jd2 {

class MpiFileBuffer { //this has to be singleton...
	typedef std::map< std::string, core::Size > Filenames;
	typedef std::map< int, SingleFileBufferOP > Buffers;
	typedef std::pair< time_t, int > TimeStampedChannel;
	typedef std::map< int, time_t > GarbageList;
public:
	MpiFileBuffer( core::Size file_buf_rank_ );
	virtual ~MpiFileBuffer();
	//core::Size receive_line( std::string& line );
private:
	void receive_str( core::Size slave, core::Size size, std::string& line );
	void open_channel( core::Size slave, std::string const& filename, bool append, core::Size& status ); //status to ret
	bool is_open_channel( std::string const& filename );
	void store_to_channel( core::Size slave, core::Size channel, std::string const& line );
	void show_status( std::ostream& ) const;
	void flush_channel( core::Size slave, core::Size channel_id );
	void close_channel( core::Size slave, core::Size channel_id );
	void block_file( core::Size slave, std::string const& filename ); //don't write to this file, and close-reopen_append streams
	void close_file( core::Size channel );
	bool remote_close_file( std::string const& filename ); //this file is no longer used... close
	//helper routine -- call when manually closing file or re-opening channel
	void clear_channel_from_garbage_collector( core::Size channel );

public:
	void release_file( std::string filename );
	void block_file( std::string const& filename );
	bool close_file( std::string fname );
  void run();
	void stop();
	void set_SlaveCanOpenFile( bool setting = true ) {
		bSlaveCanOpenFile_ = setting ;
	}
	void garbage_collection();

protected:
	virtual SingleFileBufferOP generate_new_channel( std::string const& filename, core::Size channel, bool append, core::Size& status ) = 0;
private:
  core::Size buffer_rank_;
	core::Size my_rank_;
	Filenames open_files_;
	Buffers open_buffers_;
	core::Size last_channel_;

	//if false slaves can only connect to already opened files.
	// dormant files will never be closed...
	bool bSlaveCanOpenFile_;
	bool bKeepFilesAlive_; //don't close files when no slaves want to write... probably speed up because not always reading from start
	time_t seconds_to_keep_files_alive_;
	bool bStop_;
	GarbageList garbage_collector_;
	time_t last_garbage_collection_;
	std::list< std::string > blocked_files_;
};


class WriteOut_MpiFileBuffer : public  MpiFileBuffer {
public:
	WriteOut_MpiFileBuffer( core::Size rank ) : MpiFileBuffer ( rank ) { };
protected:
	virtual SingleFileBufferOP generate_new_channel( std::string const& filename, core::Size channel, bool append, core::Size &status );
};

class DebugOut_MpiFileBuffer : public  MpiFileBuffer {
public:
	DebugOut_MpiFileBuffer( core::Size rank ) : MpiFileBuffer ( rank ) { };
protected:
	virtual SingleFileBufferOP generate_new_channel( std::string const& filename, core::Size channel, bool append, core::Size &status );
};

}
}

#endif
