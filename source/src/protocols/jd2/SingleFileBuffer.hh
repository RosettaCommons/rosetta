// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MpiFileBuffer.hh
/// @brief  header file for MPISilentFileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @detail this outputter will send silentstructs via MPI to dedicated node that will collect all structures
/// @author Oliver Lange olange@u.washington.edu


#ifndef INCLUDED_protocols_jd2_SingleFileBuffer_hh
#define INCLUDED_protocols_jd2_SingleFileBuffer_hh

#ifdef USEMPI
#include <mpi.h>
#endif

//unit headers
#include <protocols/jd2/SingleFileBuffer.fwd.hh>

//project headers
#include <core/types.hh>

//utility headers
#include <utility>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <string>
#include <map>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace jd2 {

/// @details this is a implementation of Buffer for silent-file-based output.

class SingleFileBuffer : public utility::pointer::ReferenceCount {
protected:
	typedef utility::vector1< std::string > LineBuffer;
	typedef std::map< int, LineBuffer> BufferMap;
public:
	SingleFileBuffer( std::string  filename, core::Size channel, core::Size&  status ) : filename_(std::move( filename )), mpi_channel_( channel ) { status = 0; };
	~SingleFileBuffer() override;
	void flush( core::Size slave );
	void store_line( core::Size slave, core::Size channel, std::string const & line );
	core::Size length( core::Size slave );
	void close( core::Size slave );
	bool has_open_slaves() const;
	core::Size nr_open_slaves() const;
	//a sequential block of lines hat come from a single slave node (i.e., a full-silent struct )
	// called when flush() or close() is called on the stream on the slave node
	virtual void write_lines( LineBuffer const & );
	virtual void block( core::Size slave );
	std::string const & filename() { return filename_; };
private:
	std::string filename_;
	core::Size mpi_channel_;
	BufferMap unfinished_blocks_;
};

class WriteFileSFB : public SingleFileBuffer {
public:
	typedef SingleFileBuffer Base;
	~WriteFileSFB() override;
	WriteFileSFB( std::string const & filename, core::Size channel, bool append, core::Size& status );
	void write_lines( LineBuffer const & ) override;
	void block( core::Size slave ) override;
private:
	/// @brief PyRosetta workaround, make copy constructor private so it will not try to create copy methods
	WriteFileSFB( WriteFileSFB const & src );

	std::ofstream out_;
};

}
}

#endif
