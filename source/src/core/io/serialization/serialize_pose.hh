// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author will sheffler

#ifndef INCLUDED_core_io_serialization_serialize_pose_hh
#define INCLUDED_core_io_serialization_serialize_pose_hh

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

#include <cstring>
#include <vector>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace serialization {

struct BUFFER
{
	BUFFER(size_t size) : start_(0),end_(0),size_(size),ownbuf_(true) {
		buf_ = new char[size_];
	}
	BUFFER(char * buf, std::size_t size) :
		buf_(buf),start_(0),end_(0),size_(size),ownbuf_(false) { }
	~BUFFER() {
		if ( ownbuf_ ) delete [] buf_;
	}
	int write(char * x, std::size_t nchar) {
		if ( end_+nchar >= size_ ) return -1;
		// strncpy(buf_+end_,x,nchar);
		memcpy(buf_+end_,x,nchar);
		end_ += nchar;
		return 1;
	}
	int read(char * out_buf, std::size_t nchar) {
		//if ( start_+nchar > end_ ) return -1;
		if ( start_+nchar >= size_ ) return -2;
		//  strncpy(out_buf,buf_+start_,nchar);
		memcpy(out_buf,buf_+start_,nchar);
		start_ += nchar;
		return 1;
	}
private:
	char *buf_;
	size_t start_,end_,size_;
	bool ownbuf_;
};

// stuff from interactive/util/binary_file.hh/cc
void write_binary(char           x, BUFFER & buf);
void  read_binary(char         & x, BUFFER & buf);
void write_binary(bool           x, BUFFER & buf);
void  read_binary(bool         & x, BUFFER & buf);
void write_binary(float          x, BUFFER & buf);
void  read_binary(float        & x, BUFFER & buf);
void write_binary(double         x, BUFFER & buf);
void  read_binary(double       & x, BUFFER & buf);
void write_binary(unsigned int   x, BUFFER & buf);
void  read_binary(unsigned int & x, BUFFER & buf);

/// Read/write simple structure to a file.
void write_binary(const utility::vector1_bool & x, BUFFER & buf);
void  read_binary(      utility::vector1_bool & x, BUFFER & buf);
void write_binary(const std::vector<std::string> & x, BUFFER & buf);
void  read_binary(      std::vector<std::string> & x, BUFFER & buf);
void write_binary(          const std::string & x, BUFFER & buf);
void  read_binary(                std::string & x, BUFFER & buf);
void write_binary(         const core::Vector & x, BUFFER & buf);
void  read_binary(               core::Vector & x, BUFFER & buf);

/// Utility read/write.
void check_binary_unsigned_int(unsigned int x, BUFFER & buf);
void write_binary_chars(const char *x, BUFFER & buf);
void check_binary_chars(const char *x, BUFFER & buf);

/// Read/Write a pose to a file
void write_binary(const core::pose::Pose & pose, BUFFER & buf);
void read_binary(core::pose::Pose & pose, BUFFER & buf);

} // serialization
} // io
} // core

#endif // INCLUDED_core_io_serialization_serialize_pose_HH
