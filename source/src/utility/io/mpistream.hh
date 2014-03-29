// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/mpistream.hpp
/// @detail this is a stream class to communicate with MpiFileBuffer which should be in his .run() loop on the node given by master_rank_
/// all file-opening, file-closing and writing is rerouted via MPI to the dedicated FileBuf process
/// this source is copy-pasted from zipstream.hpp...


#ifndef INCLUDED_utility_io_mpistream_hh
#define INCLUDED_utility_io_mpistream_hh

#ifdef USEMPI
#include <mpi.h>
#endif

// C++ headers
#include <algorithm>
#include <iostream>
#include <vector>
#include <sstream>
#include <ios>

namespace utility {
namespace io {
namespace mpi_stream {


/// @brief Default gzip buffer size, change this to suite your needs
const std::size_t default_buffer_size = 921600; // Was 102400; Was 4096;
const int MPI_STREAM_TAG = 42; //should be a unique number...

/// messages to send to MpiFileBuffer
enum MPI_STREAM_MSG {
	MPI_STREAM_OPEN = 1,
	MPI_STREAM_OPEN_APPEND,
	MPI_STREAM_SEND,
	MPI_STREAM_FLUSH,
	MPI_STREAM_CLOSE,
	MPI_STREAM_FILE_EXIST
};

///reported file status after opening
enum MPI_FILE_STATUS {
	MPI_SUCCESS_NEW = 1, //file didn't exist and was opened
	MPI_SUCCESS_APPEND,//file existed and will be appended
	MPI_FAIL //file didn't exist
};
/// @brief A stream decorator that takes raw input and zips it to a ostream.
/// @note  The class wraps up the inflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/
template<
  typename Elem,
  typename Tr = std::char_traits< Elem >,
  typename ElemA = std::allocator< Elem >,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator< ByteT >
  >
class basic_mpi_streambuf :
  public std::basic_streambuf< Elem, Tr >
{

public:

  typedef  std::basic_streambuf< Elem, Tr >  basic_streambuf_type;
  typedef  std::basic_ostream< Elem, Tr > &  ostream_reference;
  typedef  Elem  char_type;
  typedef  ElemA  char_allocator_type;
  typedef  ByteT  byte_type;
  typedef  ByteAT  byte_allocator_type;
  typedef  byte_type *  byte_buffer_type;
  typedef  std::vector< byte_type, byte_allocator_type >  byte_vector_type;
  typedef  std::vector< char_type, char_allocator_type >  char_vector_type;
  typedef  Tr  traits_type;
  typedef  typename Tr::int_type  int_type;

  using basic_streambuf_type::epptr;
  using basic_streambuf_type::pbase;
  using basic_streambuf_type::pptr;

  /// @brief Construct a mpi stream
  /// @note  More info on the following parameters can be found in the zlib documentation
  basic_mpi_streambuf(
			std::string filename,
      std::size_t buffer_size_,
		  int master_rank_,
			bool append
  );

  virtual ~basic_mpi_streambuf();

  int sync();
  int_type overflow( int_type c );

  /// @brief flushes the mpi buffer and output buffer
  virtual std::streamsize flush() {
		return flush( false );
	}


  virtual std::streamsize flush_final() {
		return flush( true );
	}

  /// @brief resets the mpi stream and zeros the crc
  /// @details This method should be called after flush_finalize()
  /// @deatils to allow future writes
  void reset_state();
	int file_status() const { return file_status_; }
	//void release_file();
	void print_header( std::string const& );
private:
	virtual std::streamsize flush( bool final );
  bool send_to_master( char_type*, std::streamsize );
  std::size_t fill_input_buffer();

	char_vector_type m_buffer;

	int channel_id_;
	int my_rank_;
	int master_rank_;
	int file_status_;
	//std::string filename_;
}; // basic_mpi_streambuf


/// @brief Base class for mpi ostreams
/// @note  Contains a basic_mpi_streambuf
template<
  typename Elem,
  typename Tr = std::char_traits< Elem >,
  typename ElemA = std::allocator< Elem >,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator< ByteT >
  >
class basic_mpi_ostreambase :
  virtual public std::basic_ios< Elem, Tr >
{

public:

  typedef  std::basic_ostream<Elem, Tr> &  ostream_reference;
  typedef  basic_mpi_streambuf<
    Elem,
    Tr,
    ElemA,
    ByteT,
    ByteAT
    >  mpi_streambuf_type;

  /// @brief Construct a mpi stream
  /// @note  More info on the following parameters can be found in the zlib documentation.
  basic_mpi_ostreambase(
				std::string filename,
				size_t buffer_size_,
				int master_rank,
				bool append
  ) :
    m_buf( filename, buffer_size_, master_rank, append )
  {
    this->init( &m_buf );
  }

  /// @brief returns the underlying mpi ostream object
  mpi_streambuf_type * rdbuf() { return &m_buf; }
	int file_status() const { return m_buf.file_status(); };
	void release_file() {
		m_buf.release_file();
	}
	void print_header( std::string const& header ) {
		m_buf.print_header( header );
	}
private:

  mpi_streambuf_type m_buf;

}; // basic_mpi_ostreambase

// basic_mpi_istreambase


/// @brief A mpiper ostream
///
/// @remarks
///
/// This class is a ostream decorator that behaves 'almost' like any other ostream.
///
/// At construction, it takes any ostream that shall be used to output of the compressed data.
///
/// When finished, you need to call the special method zflush or call the destructor
/// to flush all the intermidiate streams.
///
/// Example:
/// \code
/// // creating the target mpi string, could be a fstream
/// ostringstream ostringstream_;
/// // creating the mpi layer
/// mpi_ostream mpiper(ostringstream_);
///
///
/// // writing data
/// mpiper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
/// // mpi ostream needs special flushing...
/// mpiper.zflush();
/// \endcode
template<
  typename Elem,
  typename Tr = std::char_traits< Elem >,
  typename ElemA = std::allocator< Elem >,
  typename ByteT = unsigned char,
  typename ByteAT = std::allocator< ByteT >
  >
class basic_mpi_ostream :
  public basic_mpi_ostreambase< Elem, Tr, ElemA, ByteT, ByteAT >,
  public std::basic_ostream< Elem, Tr >
{

public:

  typedef  basic_mpi_ostreambase< Elem, Tr, ElemA, ByteT, ByteAT >  mpi_ostreambase_type;
  typedef  std::basic_ostream< Elem, Tr >  ostream_type;
  typedef  std::basic_ostream< Elem, Tr > &  ostream_reference;
  typedef  Elem  char_type;

  using ostream_type::flush;
  using mpi_ostreambase_type::rdbuf;

  /// @brief Constructs a mpiper ostream decorator
  ///
  /// @param ostream_ ostream where the compressed output is written
  /// @param is_gmpi_ true if gmpi header and footer have to be added
  /// @param level_ level of compression 0, bad and fast, 9, good and slower,
  /// @param strategy_ compression strategy
  /// @param window_size_ see zlib doc
  /// @param memory_level_ see zlib doc
  /// @param buffer_size_ the buffer size used to mpi data
  ///
  /// @note  When is_gmpi_ is true, a gmpi header and footer is automatically added
  basic_mpi_ostream(
				std::string filename,						// int open_mode = std::ios::out,
				int master_rank,
				std::stringstream& header,
				bool append = false,
		    std::size_t buffer_size_ = default_buffer_size
  ) :
    mpi_ostreambase_type(
			 filename,
			 buffer_size_,
			 master_rank,
			 append
    ),
    ostream_type( rdbuf() ),
    m_mpi_stream_finalized( false )
  {
		if ( mpi_ostreambase_type::file_status() == MPI_SUCCESS_NEW ) {
			//print_header();
			//			(*this) << header.str();
			//			mpi_ostreambase_type::release_file();
			mpi_ostreambase_type::print_header( header.str() );
		} else if ( mpi_ostreambase_type::file_status() == MPI_FAIL ) {
			// Set failbit so failure can be detected
			ostream_type::setstate( std::ios_base::failbit );
		}
	};

	void close() {
	};

  ~basic_mpi_ostream()
  {
    // adding a footer is not necessary here, as it will be
    // taken care of during the last zflush_finalize()
    // called by the higher level close() routines
  }

  /// @brief stream output
  /// @details if mpi stream has been finalized, will reset
  /// @details the stream and add header if necessary
  template< typename T >
  inline
  basic_mpi_ostream &
  operator <<( T const & t )
  {
    static_cast< std::ostream & >( *this ) << t;
    return *this;
  }

  /// @brief write char
  /// @details if mpi stream has been finalized, will reset
  /// @details the stream and add header if necessary
  inline
  basic_mpi_ostream &
  put( char const c )
  {
		static_cast< std::ostream & >( *this ).put( c );
    return *this;
  }

  /// @brief write a string
  /// @details if mpi stream has been finalized, will reset
  /// @details the stream and add header if necessary
  inline
  basic_mpi_ostream &
  write( char const * str, std::streamsize const count )
  {
		static_cast< std::ostream & >( *this ).write( str, count );
    return *this;
  }

	inline
	basic_mpi_ostream &
	flush() {
		ostream_type::flush();
		// std::basic_ostream::flush() only calls rdbuf()->sync() not rdbuf()->flush()
		// therefore, include explicit call below to trigger MpiFileBuffer::flush_channel()
		rdbuf()->flush();
		return *this;
	}

private:


  static void put_long_as_uint32( ostream_reference out_, unsigned long x_ );

  /// @brief tracks to see if mpi stream was finalized
  /// @details set to true during zflush_finalize()
  /// @details set to false during reset_state()
  bool m_mpi_stream_finalized;

}; // basic_mpi_ostream


// Types
typedef  basic_mpi_ostream< char >  mpi_ostream;
typedef  basic_mpi_ostream< wchar_t >  mpi_wostream;


} // mpistream
} //io
} //utility

// Implementation [adding extra space before #include so PyRosetta skipp it...]
 #include <utility/io/mpistream.ipp>


#endif // INCLUDED_utility_io_mpistream_HPP

