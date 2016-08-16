// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#ifndef INCLUDED_utility_io_mpistream_IPP
#define INCLUDED_utility_io_mpistream_IPP

// Unit headers
#include <utility/io/zipstream.hpp>

#include <utility/assert.hh> //for MPI_ONLY macro
// C++ headers
#include <sstream>

namespace utility {
namespace io {
namespace mpi_stream {

template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::basic_mpi_streambuf(
	std::string MPI_ONLY( filename ),
	size_t buffer_size_,
	int master_rank,
	bool MPI_ONLY( append )
) : m_buffer( buffer_size_, 0 ), master_rank_( master_rank )
{
	//find channel_id_;
	this->setp( &(m_buffer[0]), &(m_buffer[m_buffer.size()-1]) );

#ifdef USEMPI
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank_);/* get current process id */
	//	std::cerr << "open " << filename << " from client " << my_rank_ << std::endl;
	//establish connection with master
	// akin to "open file"
	int buf[ 4 ];
  buf[ 0 ] = my_rank_;
  buf[ 1 ] = filename.size();
	buf[ 2 ] = append ? MPI_STREAM_OPEN_APPEND : MPI_STREAM_OPEN;
	buf[ 3 ] = 0;

  MPI_Send(buf, 4, MPI_INT, master_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
	MPI_Send(const_cast<char*> (filename.data()), filename.size(), MPI_CHAR, master_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
	MPI_Status stat;

	MPI_Recv(&buf, 2, MPI_INT, master_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD, &stat );
	channel_id_ = buf[ 0 ];
	file_status_ = buf[ 1 ];
	if ( file_status_ == MPI_FAIL ) {
		std::cerr << "ERROR when opening mpistream to write to " << filename << std::endl;
	}

#endif
}

template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
void
basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::print_header( std::string const& MPI_ONLY(header) ) {
#ifdef USEMPI
		int buf[ 4 ];
		buf[ 0 ] = my_rank_;
		buf[ 1 ] = header.size();
		buf[ 2 ] = MPI_STREAM_SEND;
		buf[ 3 ] = channel_id_;
		int master_rank = 0;

		//		std::cerr << "sending from client " << my_rank_ << std::endl;
		MPI_Send(buf, 4, MPI_INT, master_rank, MPI_STREAM_TAG, MPI_COMM_WORLD );
		MPI_Send(const_cast<char*> (header.data()), header.size(), MPI_CHAR, master_rank, MPI_STREAM_TAG, MPI_COMM_WORLD ); //MPI_CHAR OR MPI_INT ?
#endif
}

template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::~basic_mpi_streambuf()
{
	//  std::cerr << "destruct mpi_streambuf on node " << my_rank_ << std::endl;
	flush_final();//communicate closing
}

////SYNC
template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
int
basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::sync()
{
	if ( pptr() && pptr() > pbase() ) {
		if ( traits_type::eq_int_type( overflow( traits_type::eof() ), traits_type::eof() ) ) return -1;
	}

	return 0;
}


///OVERFLOW
template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
typename basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::int_type
basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::overflow(
	typename basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::int_type c
)
{
	bool const test_eof = traits_type::eq_int_type( c, traits_type::eof() );
	int w = static_cast< int >( pptr() - pbase() );
	if ( !test_eof ) {
		*pptr() = c;
		++w;
	}
	if ( send_to_master( pbase(), w ) ) {
		this->setp( pbase(), epptr() - 1 );
		return traits_type::not_eof( c );
	} else {
		return traits_type::eof();
	}
}

///+++ SEND_TO_MASTER ++++
template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
bool
basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::send_to_master(
	typename basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::char_type *
#ifdef USEMPI
buffer_
#endif
	,std::streamsize
#ifdef USEMPI
 buffer_size_
#endif
)
{
	//  std::streamsize written_byte_size = 0, total_written_byte_size = 0;


#ifdef USEMPI
		byte_buffer_type next_out = (byte_buffer_type)buffer_;
		uInt avail_out = static_cast< uInt >( buffer_size_ * sizeof(char_type) );

		//size_t remainder = 0;

		int buf[ 4 ];
		buf[ 0 ] = my_rank_;
		buf[ 1 ] = avail_out;
		buf[ 2 ] = MPI_STREAM_SEND;
		buf[ 3 ] = channel_id_;
		int master_rank = 0;

		//		std::cerr << "sending from client " << my_rank_ << std::endl;
		MPI_Send(buf, 4, MPI_INT, master_rank, MPI_STREAM_TAG, MPI_COMM_WORLD );
		MPI_Send(next_out, avail_out, MPI_CHAR, master_rank, MPI_STREAM_TAG, MPI_COMM_WORLD ); //MPI_CHAR OR MPI_INT ?
 		return true; // success detection ?

#else
	return false;
#endif
}

template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
std::streamsize basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::flush( bool
#ifdef USEMPI
 final
#endif
)
{
	//  std::streamsize written_byte_size = 0,


	int const buffer_size = static_cast< int >( pptr() - pbase() ); // amount of data currently in buffer
	send_to_master( pbase(), buffer_size );
	std::streamsize total_written_byte_size = buffer_size;
#ifdef USEMPI
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank_);/* get current process id */
	//establish connection with master
	// akin to "open file"
	int buf[ 4 ];
  buf[ 0 ] = my_rank_;
  buf[ 1 ] = 0;
	buf[ 2 ] = final ? MPI_STREAM_CLOSE : MPI_STREAM_FLUSH;
	buf[ 3 ] = channel_id_;

  MPI_Send(buf, 4, MPI_INT, master_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );

#endif

	return total_written_byte_size;
}


//  template<
//   typename Elem,
//   typename Tr,
//   typename ElemA,
//   typename ByteT,
//   typename ByteAT
//  >
//  void basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::release_file()
//  {
//   #ifdef USEMPI
//   int buf[ 4 ];
//   buf[ 0 ] = my_rank_;
//   buf[ 1 ] = filename_.size();
//   buf[ 2 ] = MPI_RELEASE_FILE;
//   buf[ 3 ] = 0;
//   MPI_Send(buf, 4, MPI_INT, master_rank_, MPI_STREAM_TAG, MPI_COMM_WORLD );
//   #endif
//  }


template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>
void basic_mpi_streambuf< Elem, Tr, ElemA, ByteT, ByteAT >::reset_state()
{
	std::cout << "call to reset_state" << std::endl;
}


template<
typename Elem,
typename Tr,
typename ElemA,
typename ByteT,
typename ByteAT
>

void basic_mpi_ostream< Elem, Tr, ElemA, ByteT, ByteAT >::put_long_as_uint32(
	typename basic_mpi_ostream< Elem, Tr, ElemA, ByteT, ByteAT >::ostream_reference out_,
	unsigned long x_
)
{
	// yab: 20090414, modified to conform to gmpi standard where
	// trailer crc and length must both be 32-bit, otherwise there
	// is breakage in systems where 'unsigned long' is not 32-bit
	// and external archiving programs end up complaining.
	char b1, b2, b3, b4; // assuming char is 8 bits
	b1 = 0xFF & x_;
	b2 = 0xFF & ( x_ >> 8 );
	b3 = 0xFF & ( x_ >> 16 );
	b4 = 0xFF & ( x_ >> 24 );

	out_.write( &b1, 1 );
	out_.write( &b2, 1 );
	out_.write( &b3, 1 );
	out_.write( &b4, 1 );
}

} // namespace mpi_stream
}
}

#endif // INCLUDED_utility_io_mpistream_IPP
