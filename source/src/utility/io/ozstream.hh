// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/ozstream.hh
/// @brief  Output file stream wrapper for uncompressed and compressed files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author David Kim (dekim@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


#ifndef INCLUDED_utility_io_ozstream_hh
#define INCLUDED_utility_io_ozstream_hh


#ifdef USEMPI
#include <mpi.h>
#endif
// Unit headers
#include <utility/io/ozstream.fwd.hh>

// Package headers
#include <utility/io/orstream.hh>

// Project headers
#include <utility/file/gzip_util.hh>
#if defined( USE_FILE_PROVIDER )
#include <utility/inline_file_provider.hh>
#endif

// Utility headers
#include <utility/io/zipstream.hpp>
#include <utility/io/mpistream.hh>

// C++ headers
#include <algorithm>
#include <fstream>
#include <ios>
#include <sstream>

namespace utility {
namespace io {


/// @brief default buffer size for ozstreams (900KB)
/// @note  this must be at least 4KB, otherwise zipstream will break
std::streamsize const OZSTREAM_DEFAULT_BUFFER_SIZE = 921600;


/// @brief ozstream: Output file stream wrapper for uncompressed and compressed files
class ozstream :
	public orstream
{


private: // Friends


	friend long utility::file::gzip( std::string const & uncompressedfile, bool overwrite );
	friend long utility::file::gunzip( std::string const & compressedfile, bool overwrite );


private: // Types


	enum Compression { NONE, UNCOMPRESSED, GZIP };


public: // Creation


	/// @brief Default constructor
	inline
	ozstream() :
		compression_( NONE ),
		buffer_size_( OZSTREAM_DEFAULT_BUFFER_SIZE ),
		char_buffer_p_( NULL ),
		zip_stream_p_( 0 ),
		mpi_stream_p_( 0 )
#if defined( USE_FILE_PROVIDER )	
		,file_provider_stream( &bad_stream ) 
#endif

	{}


	/// @brief Filename constructor
	/// @param [in] filename_a filename
	/// @param [in] open_mode  opening mode bitmask, use std::ios_base::out for gzip files
	/// @param [in] buf_size   buffer size (bytes), enforced lower bound of 4KB
	inline
	explicit
	ozstream(
		std::string const & filename_a,
		std::ios_base::openmode open_mode = std::ios_base::out, // Ignored for gzip files
		std::streamsize buf_size = OZSTREAM_DEFAULT_BUFFER_SIZE
	) :
		compression_( NONE ),
		char_buffer_p_( NULL ),
		zip_stream_p_( 0 ),
		mpi_stream_p_( 0 )
#if defined( USE_FILE_PROVIDER )		
		,file_provider_stream( &bad_stream )
#endif

	{
		buffer_size( buf_size ); // set buffer
		open( filename_a, open_mode );
	}


	/// @brief Destructor
	inline
	virtual
	~ozstream()
	{
		close();
	}


public: // Methods: conversion


	/// @brief bool conversion
	inline
	operator bool() const
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        return true;
		  }
    #endif
		return ( zip_stream_p_ ? !zip_stream_p_->fail() : ( mpi_stream_p_ ? !mpi_stream_p_->fail() : !!of_stream_ ));
	}


	/// @brief Stream conversion
	inline
	operator std::ostream const &() const
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        return *file_provider_stream;
		  }
    #endif
		return ( zip_stream_p_
		 ? static_cast< std::ostream const & >( *zip_stream_p_ )
			: ( mpi_stream_p_ ? static_cast< std::ostream const& > ( *mpi_stream_p_ )
				: static_cast< std::ostream const & >( of_stream_ ) )
		);
	}


	/// @brief Stream conversion
	inline
	operator std::ostream &()
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        return *file_provider_stream;
		  }
    #endif
		return ( zip_stream_p_
		 ? static_cast< std::ostream & >( *zip_stream_p_ )
			: ( mpi_stream_p_ ? static_cast< std::ostream & > ( *mpi_stream_p_ )
				: static_cast< std::ostream & >( of_stream_ ) )
		);
	}


public: // Methods: formatting


	/// @brief Stream output: override to preserve type of return value
	template< typename T >
	inline
	ozstream &
	operator <<( T const & t )
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        (*file_provider_stream) << t;
			  return *this;
		  }
    #endif
		if ( zip_stream_p_ ) {
			(*zip_stream_p_) << t;
		} else if ( mpi_stream_p_ ) {
			(*mpi_stream_p_) << t;
		} else {
			of_stream_ << t;
		}
		return *this;
	}


	/// @brief Stream output overload to replace std::endl with \n to avoid flushing
	/// @brief and call ozstream::flush() when passed std::flush
	inline
	ozstream &
	operator <<( manipulator m )
	{
		static manipulator const std_endl = std::endl;
		static manipulator const std_flush = std::flush;
		if ( m == std_endl && ( mpi_stream_p_ || zip_stream_p_ ) )  {
			#if defined( USE_FILE_PROVIDER )
				if (file_provider_stream->good() ){
          (*file_provider_stream) << '\n';
				  return *this;
			  }
      #endif
			if ( zip_stream_p_ ) { // Output newline instead
				(*zip_stream_p_) << '\n';
			} else if ( mpi_stream_p_ ) {
				(*mpi_stream_p_) << '\n';
			}
		} else if ( ( m == std_flush ) && ( zip_stream_p_ || mpi_stream_p_) ) {
			#if defined( USE_FILE_PROVIDER )
				if (file_provider_stream->good() ){
          file_provider_stream->flush();
				  return *this;
			  }
      #endif
			flush(); // ozstream::flush()
		} else {
			#if defined( USE_FILE_PROVIDER )
				if (file_provider_stream->good() ){
          (*file_provider_stream) << m;
				  return *this;
			  }
      #endif
			of_stream_ << m;
		}
		return *this;
	}


public: // Methods: i/o


	/// @brief Open a file
	void
	open(
		std::string const & filename_a,
		std::ios_base::openmode open_mode = std::ios_base::out
	);


	/// @brief Open a text file for appending
	void
	open_append( std::string const & filename_a );

	/// @brief open file as append if it exists, return true if file existed before, false if it is new.
	void
	open_append_if_existed( std::string const& filename_a, std::stringstream& preprinted_header );

	/// @brief Write a char
	inline
	ozstream &
	put( char const c )
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        file_provider_stream->put( c );
			  return *this;
		  }
    #endif

		if ( zip_stream_p_ ) {
			zip_stream_p_->put( c );
		} else if ( mpi_stream_p_ ) {
			mpi_stream_p_->put( c );
		} else {
			of_stream_.put( c );
		}
		return *this;
	}


	/// @brief Write a string
	inline
	ozstream &
	write( char const * str, std::streamsize const count )
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        file_provider_stream->write( str, count );
			  return *this;
		  }
    #endif
		if ( zip_stream_p_ ) {
			zip_stream_p_->write( str, count );
		} else if ( mpi_stream_p_ ) {
			mpi_stream_p_->write( str, count );
		} else {
			of_stream_.write( str, count );
		}
		return *this;
	}


	/// @brief Write a string
	inline
	ozstream &
	write( std::string const & str, std::streamsize const count )
	{
		stream().write( str.c_str(), count );
		return *this;
	}


	/// @brief Flush the stream -- currently alias to flush_finalize()
	/// @details Instead doing a regular flush, we currently force a
	///          completion of the zip stream.  We do this to pre-empt
	///          against several things: (1) the expensive operation
	///          of closing and re-opening a stream (2) confusion and
	///          inconvenience that may result from users calling flushes
	///          and ending upon with broken or corrupted files if a
	///          job is somehow interrupted (e.g. on a cluster).
	///          Please note that calling flush_finalize() too often
	///          can seriously degrade compression.  zlib documentation
	///          appears to imply that every 1MB or so is a reasonable
	///          rule of thumb.
	inline
	ozstream &
	flush()
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        file_provider_stream->flush();
			  return *this;
		  }
    #endif
		// comment out the zflush_finalize() containing line and uncomment
		// the flush() containing line to switch to "regular" flush behavior
		if ( zip_stream_p_ ) zip_stream_p_->zflush_finalize();
//		if ( zip_stream_p_ ) zip_stream_p_->zflush();
		if ( mpi_stream_p_ ) {
			mpi_stream_p_->flush();
			return *this; //no of_stream_ if mpi_stream
		}
		of_stream_.flush();
		return *this;
	}


	/// @brief Flush the streams and finalize the zip stream
	/// @details Calling this will complete the zip stream.
	///          Upon the next write, a new zip stream will be started.
	///          Please note that calling flush_finalize() too often
	///          can seriously degrade compression.  zlib documentation
	///          appears to imply that every 1MB or so is a reasonable
	///          rule of thumb.
	inline
	ozstream &
	flush_finalize()
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        file_provider_stream->flush();
			  return *this;
      }
		#endif
		if ( zip_stream_p_ ) zip_stream_p_->zflush_finalize();
		if ( mpi_stream_p_ ) {
			mpi_stream_p_->flush();
			return *this; //no of_stream_ if mpi_stream
		}
		of_stream_.flush();
		return *this;
	}


	/// @brief Flush the zip_ostream
	/// @details this will flush but *not* finalize the zip stream
	inline
	void
	zflush()
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        return;
		  }
    #endif
		if ( zip_stream_p_ ) zip_stream_p_->zflush();
	}


	/// @brief Flush and finalize the zip_stream
	/// @details Calling this will complete the zip stream.
	///          Upon the next write, a new zip stream will be started
	///          Please note that calling zflush_finalize() too often
	///          can seriously degrade compression.  zlib documentation
	///          appears to imply that every 1MB or so is a reasonable
	///          rule of thumb.
	inline
	void
	zflush_finalize()
	{
		#if defined( USE_FILE_PROVIDER )
      if (file_provider_stream->good() ){
        return;
      }
		#endif
		if ( zip_stream_p_ ) zip_stream_p_->zflush_finalize();
	}

	/// @brief Clear the stream
	inline
	void
	clear()
	{
		#if defined( USE_FILE_PROVIDER )
			if (file_provider_stream->good() ){
        file_provider_stream->clear();
			  return;
		  }
    #endif
		of_stream_.clear();
		if ( zip_stream_p_ ) zip_stream_p_->clear();
		if ( mpi_stream_p_ ) mpi_stream_p_->clear();
	}


	/// @brief Close the ofstream and reset the state
	inline
	void
	close()
	{
		#if defined( USE_FILE_PROVIDER )
			 if (file_provider_stream->good() ){
         return;
		   }
    #endif
		if ( zip_stream_p_ ) {
			zip_stream_p_->zflush_finalize();
			delete zip_stream_p_; zip_stream_p_ = 0;
		}
		of_stream_.close();
		of_stream_.clear();
		if ( mpi_stream_p_ ) {
			mpi_stream_p_->close();
			mpi_stream_p_->clear();
			delete mpi_stream_p_;
			mpi_stream_p_ = 0;
		}
		compression_ = NONE;
		filename_.clear();

		destroy_char_buffer();
	}


public: // Properties


	/// @brief Stream access
	inline
	std::ostream const &
	operator ()() const
	{
		#if defined( USE_FILE_PROVIDER )
		  return (*file_provider_stream);
		#endif
		return ( zip_stream_p_
			? static_cast< std::ostream const & >( *zip_stream_p_ )
			: ( mpi_stream_p_ ? static_cast< std::ostream const & >( *mpi_stream_p_ )
				: static_cast< std::ostream const & >( of_stream_ ) ) );
	}


	/// @brief Stream access
	inline
	std::ostream &
	operator ()()
	{
		#if defined( USE_FILE_PROVIDER )
		   if (file_provider_stream->good() ){
         return (*file_provider_stream);
		   }
    #endif
		return ( zip_stream_p_
		 ? static_cast< std::ostream & >( *zip_stream_p_ )
			: ( mpi_stream_p_ ? static_cast< std::ostream & >( *mpi_stream_p_ )
				: static_cast< std::ostream & >( of_stream_ ) ));
	}


	/// @brief Stream access
	inline
	std::ostream const &
	stream() const
	{
		#if defined( USE_FILE_PROVIDER )
		  if (file_provider_stream->good() ){
        return (*file_provider_stream);
      }
		#endif
		return ( zip_stream_p_
		 ? static_cast< std::ostream const & >( *zip_stream_p_ )
			: ( mpi_stream_p_ ? static_cast< std::ostream const & >( *mpi_stream_p_ )
				: static_cast< std::ostream const & >( of_stream_ ) ));
	}


	/// @brief Stream access
	inline
	std::ostream &
	stream()
	{
		#if defined( USE_FILE_PROVIDER )
      if (file_provider_stream->good() ){
        return (*file_provider_stream);
      }
    #endif
		return ( zip_stream_p_
		 ? static_cast< std::ostream & >( *zip_stream_p_ )
			: ( mpi_stream_p_ ? static_cast< std::ostream & >( *mpi_stream_p_ )
				: static_cast< std::ostream & >( of_stream_ ) ));
	}


	/// @brief Pointer to the stream buffer
	inline
	std::streambuf *
	rdbuf() const
	{
		return stream().rdbuf();
	}


	/// @brief File name
	inline
	std::string const &
	filename() const
	{
		return filename_;
	}


public: // Properties: predicate


	/// @brief Good?
	inline
	bool
	good() const
	{
		return stream().good();
	}


	/// @brief End of file?
	inline
	bool
	eof() const
	{
		return stream().eof();
	}


	/// @brief Fail?
	inline
	bool
	fail() const
	{
		return stream().fail();
	}


	/// @brief Bad?
	inline
	bool
	bad() const
	{
		return stream().bad();
	}


	/// @brief Compressed?
	inline
	bool
	compressed() const
	{
		return ( compression_ == GZIP );
	}


	/// @brief Uncompressed?
	inline
	bool
	uncompressed() const
	{
		return ( compression_ == UNCOMPRESSED );
	}


	/// @brief gzipped?
	inline
	bool
	gzipped() const
	{
		return ( compression_ == GZIP );
	}


public: // Properties: buffer


	/// @brief get buffer size (bytes)
	/// @details In uncompressed mode this is the size of the character buffer.
	///          In compressed mode this is the size of the zip buffer.
	inline
	std::streamsize
	buffer_size() const
	{
		return buffer_size_;
	}


	/// @brief set buffer size (bytes)
	/// @details In uncompressed mode this is the size of the character buffer.
	///          In compressed mode this is the size of the zip buffer.
	///          Lower bound of 4KB is enforced due to zipstream requirements.
	///          Operation is skipped if the file is currently open, as the
	///          buffer is considered locked.
	void
	buffer_size(
		std::streamsize const & buf_size
	)
	{
		if ( !of_stream_.is_open() ) {
			buffer_size_ = std::max( buf_size, static_cast< std::streamsize >( 4096 ) );
		}
	}


private: // Properties: internal zip stream predicates


	/// @brief Is stream attached to a gzip file?
	inline
	bool
	is_gzip() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->is_gzip() : false );
	}


	/// @brief CRC of the uncompressed data (see zipstream documentation)
	inline
	long
	get_crc() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->get_crc() : 0 );
	}


	/// @brief Uncompressed data size
	inline
	long
	get_in_size() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->get_in_size() : 0 );
	}


	/// @brief Compressed data size
	inline
	long
	get_out_size() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->get_out_size() : 0 );
	}


private: // buffer management


	/// @brief if character buffer does not exist, create it and assign it to
	///        internal ofstream
	/// @details File must be closed for this operation to be successful,
	///          otherwise we can run into implementation dependent behavior
	///          of std::basic_filebuf.
	/// @return true if successful, false otherwise
	inline
	bool
	allocate_assign_char_buffer()
	{
		if ( !char_buffer_p_ && !of_stream_.is_open() ) {
			char_buffer_p_ = new char[ buffer_size_ ];
			of_stream_.rdbuf()->pubsetbuf( char_buffer_p_, buffer_size_ );

			return true;
		}

		return false;
	}

	/// @brief if character buffer exists, destroy it
	/// @details File must be closed for this operation to be successful,
	///          otherwise we are deleting a buffer that's still in use.
	/// @return true if successful, false otherwise
	inline
	bool
	destroy_char_buffer()
	{
		if ( char_buffer_p_ && !of_stream_.is_open() ) {
			delete [] char_buffer_p_;
			char_buffer_p_ = NULL;

			return true;
		}

		return false;
	}

public:

	static void enable_MPI_reroute( int min_client_rank, int master_rank );

	//return master( buffer ) rank or -1
	static int MPI_reroute_rank() {
#ifdef USEMPI
		return bMPI_reroute_stream_ ? mpi_FileBuf_master_rank_ : -1;
#else
		return -1;
#endif
	}
private: // Fields


	/// @brief Compression state
	Compression compression_;

	/// @brief File stream
	std::ofstream of_stream_;

	/// @brief File name
	std::string filename_;

	/// @brief size of buffer (in bytes)
	/// @details In uncompressed mode this is the size of the character buffer.
	///          In compressed mode this is the size of the zip buffer.
	///          Must be at least 4KB otherwise zipstream will break.
	std::streamsize buffer_size_;

	/// @brief character buffer pointer (owning)
	char * char_buffer_p_;

	/// @brief Zip file stream pointer (owning)
	zlib_stream::zip_ostream *zip_stream_p_;


	mpi_stream::mpi_ostream *mpi_stream_p_;

	static bool bMPI_reroute_stream_;
	static int mpi_FileBuf_master_rank_;

#if defined( USE_FILE_PROVIDER ) 
	std::ostream *file_provider_stream;
	std::stringstream bad_stream;
#endif

}; // ozstream


} // namespace io
} // namespace utility


#endif // INCLUDED_utility_io_ozstream_HH
