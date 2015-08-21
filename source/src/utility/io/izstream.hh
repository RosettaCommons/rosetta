// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/izstream.hh
/// @brief  Input file stream wrapper for uncompressed and compressed files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author David Kim (dekim@u.washington.edu)


#ifndef INCLUDED_utility_io_izstream_hh
#define INCLUDED_utility_io_izstream_hh


// Unit headers
#include <utility/io/izstream.fwd.hh>

// Package headers
#include <utility/io/irstream.hh>
#include <utility/io/zipstream.hpp>

// Project headers
#include <utility/file/gzip_util.hh>
#include <utility/vector1.hh>
#if defined( USE_FILE_PROVIDER )
#include <utility/inline_file_provider.hh>
#endif

// C++ headers
#include <fstream>
#include <istream>
#include <limits>
#include <string>
#include <iostream>

namespace utility {
namespace io {


/// @brief izstream: Input file stream wrapper for uncompressed and compressed files
class izstream :
	public irstream
{


private: // Friends


	friend long utility::file::gzip( std::string const & uncompressedfile, bool overwrite );
	friend long utility::file::gunzip( std::string const & compressedfile, bool overwrite );


private: // Types


	typedef  std::istream & (*manipulator)( izstream & );
	typedef  std::istream & (*std_manipulator)( std::istream & );

	enum Compression { NONE, UNCOMPRESSED, GZIP };


public: // Creation


	/// @brief Default constructor
	inline
	izstream() :
		compression_( NONE ),
		zip_stream_p_( 0 )
#if defined( USE_FILE_PROVIDER )
		,file_provider_stream( &bad_stream )
#endif

	{}


	/// @brief Filename constructor
	inline
	explicit
	izstream(
		std::string const & filename_a,
		std::ios_base::openmode open_mode = std::ios_base::in // Ignored for gzip files
	) :
		compression_( NONE ),
		zip_stream_p_( 0 )
#if defined( USE_FILE_PROVIDER )
		,file_provider_stream( &bad_stream )
#endif

	{
		open( filename_a, open_mode );
	}


	/// @brief Destructor
	inline
	virtual
	~izstream()
	{
		delete zip_stream_p_;
		if_stream_.close();
		if_stream_.clear();
	}


public: // Methods: conversion


	/// @brief bool conversion
	inline
	operator bool() const
	{
#if defined( USE_FILE_PROVIDER )
		// if file is present in inline file provider, return true
		if ( file_provider_stream->good() ) return true;
#endif
		//return ( zip_stream_p_ ? zip_stream_p_->good() : if_stream_.good() );
		// proper behavior is actually ( ! fail() )
		return ( zip_stream_p_ ? !zip_stream_p_->fail() : !!if_stream_ );
	}


	/// @brief Stream conversion
	inline
	operator std::istream const &() const
	{
#if defined( USE_FILE_PROVIDER )
		// if file is present in inline file provider, return that - otherwise return an actual istream
		if ( file_provider_stream->good() ) {
			return *file_provider_stream;
		}
#endif
		return ( zip_stream_p_
			? static_cast< std::istream const & >( *zip_stream_p_ )
			: static_cast< std::istream const & >( if_stream_ ) );
	}


	/// @brief Stream conversion
	inline
	operator std::istream &()
	{
#if defined( USE_FILE_PROVIDER )
		if ( file_provider_stream->good() ) {
			return *file_provider_stream;
		}
#endif
		return ( zip_stream_p_
			? static_cast< std::istream & >( *zip_stream_p_ )
			: static_cast< std::istream & >( if_stream_ ) );
	}


public: // Methods: formatting


	/// @brief Stream input
	template< typename T >
	inline
	std::istream &
	operator >>( T & t )
	{
		return stream() >> t;
	}


	/// @brief Stream manipulator input
	inline
	std::istream &
	operator >>( manipulator m )
	{
		return m( *this );
	}


	/// @brief Stream manipulator input
	inline
	std::istream &
	operator >>( std_manipulator m )
	{
		return m( *this );
	}


public: // Methods: i/o


	/// @brief Open a file
	void
	open(
		std::string const & filename_a,
		std::ios_base::openmode open_mode = std::ios_base::in
	);


	/// @brief Clear the stream(s)
	inline
	void
	clear()
	{
#if defined( USE_FILE_PROVIDER )
		// if the file is coming from the file provider then clear that stream. Otherwise go on to clear the actual file streams.
		if ( file_provider_stream->good() ) {
			file_provider_stream->clear();
			return;
		}
#endif
		if_stream_.clear();
		if ( zip_stream_p_ ) zip_stream_p_->clear();
	}


	/// @brief Close the ifstream and reset the state
	inline
	void
	close()
	{
#if defined( USE_FILE_PROVIDER )
		// no need to do anything if file is comign from file provider and not from disk
		if ( file_provider_stream->good() ) {
			return;
		}
#endif
		compression_ = NONE;
		if_stream_.close();
		if_stream_.clear();
		filename_.clear();
		delete zip_stream_p_; zip_stream_p_ = 0;
	}


	/// @brief Seek to the beginning
	inline
	void
	seek_beg()
	{
#if defined( USE_FILE_PROVIDER )
		if ( file_provider_stream->good() ) {
			file_provider_stream->clear();
			file_provider_stream->seekg( std::ios_base::beg );
			file_provider_stream->clear();
			return;
		}
#endif
		if_stream_.clear();
		if_stream_.seekg( std::ios_base::beg );
		if_stream_.clear();
		if ( zip_stream_p_ ) {
			delete zip_stream_p_; zip_stream_p_ = new zlib_stream::zip_istream( if_stream_ );
			if ( ( !zip_stream_p_ ) || ( !( *zip_stream_p_ ) ) || ( !zip_stream_p_->is_gzip() ) ) {
				delete zip_stream_p_; zip_stream_p_ = 0;
				if_stream_.close();
				if_stream_.setstate( std::ios_base::failbit ); // set ios state to failbit
			}
		}
	}


	/// @brief Get the next character
	inline
	int
	get()
	{
		return stream().get();
	}


	/// @brief Get the next character
	inline
	izstream &
	get( char & c )
	{
		stream().get( c );
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	izstream &
	get( char * str, std::streamsize const count )
	{
		stream().get( str, count );
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	izstream &
	get( char * str, std::streamsize const count, char const delim )
	{
		stream().get( str, count, delim );
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	izstream &
	get( std::string & str, std::streamsize const count )
	{
		char * cp = new char[ count ];
		stream().get( cp, count );
		str = cp;
		delete[] cp;
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	izstream &
	get( std::string & str, std::streamsize const count, char const delim )
	{
		char * cp = new char[ count ];
		stream().get( cp, count, delim );
		str = cp;
		delete[] cp;
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	izstream &
	getline( char * line, std::streamsize const count )
	{
		stream().getline( line, count );
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	izstream &
	getline( char * line, std::streamsize const count, char const delim )
	{
		stream().getline( line, count, delim );
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	izstream &
	getline( std::string & line )
	{
		std::getline( stream(), line );
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	izstream &
	getline( std::string & line, char const delim )
	{
		std::getline( stream(), line, delim );
		return *this;
	}


	/// @brief Read the next specified number of characters
	inline
	izstream &
	read( char * str, std::streamsize const count )
	{
		stream().read( str, count );
		return *this;
	}


	/// @brief Read the next specified number of characters
	inline
	izstream &
	read( std::string & str, std::streamsize const count )
	{
		char * cp = new char[ count ];
		stream().read( cp, count );
		str = cp;
		delete[] cp;
		return *this;
	}


	/// @brief Read the next available specified number of characters
	inline
	std::streamsize
	readsome( char * str, std::streamsize const count )
	{
		return stream().readsome( str, count );
	}


	/// @brief Read the next available specified number of characters
	inline
	std::streamsize
	readsome( std::string & str, std::streamsize const count )
	{
		char * cp = new char[ count ];
		std::streamsize const n_chars = stream().readsome( cp, count );
		str = cp;
		delete[] cp;
		return n_chars;
	}


	/// @brief Skip over the next character
	inline
	izstream &
	ignore()
	{
		stream().ignore();
		return *this;
	}


	/// @brief Skip over the next specified number of characters
	inline
	izstream &
	ignore( std::streamsize const count )
	{
		stream().ignore( count );
		return *this;
	}


	/// @brief Skip over the next specified number of characters
	inline
	izstream &
	ignore( std::streamsize const count, char const delim )
	{
		stream().ignore( count, delim );
		return *this;
	}


	/// @brief Returns the next character without extracting it
	inline
	int
	peek()
	{
		return stream().peek();
	}


	/// @brief Put the last character read back into the stream
	inline
	izstream &
	unget()
	{
		stream().unget();
		return *this;
	}


	/// @brief Put the last character read back into the stream and check
	///        that passed character is correct
	inline
	izstream &
	putback( char c )
	{
		stream().putback( c );
		return *this;
	}


public: // Properties


	/// @brief Stream access
	inline
	std::istream const &
	operator ()() const
	{
#if defined( USE_FILE_PROVIDER )
		if ( file_provider_stream->good() ) {
			return *file_provider_stream;
		}
#endif
		return ( zip_stream_p_
			? static_cast< std::istream const & >( *zip_stream_p_ )
			: static_cast< std::istream const & >( if_stream_ ) );
	}


	/// @brief Stream access
	inline
	std::istream &
	operator ()()
	{
#if defined( USE_FILE_PROVIDER )
		if ( file_provider_stream->good() ) {
			return *file_provider_stream;
		}
#endif
		return ( zip_stream_p_
			? static_cast< std::istream & >( *zip_stream_p_ )
			: static_cast< std::istream & >( if_stream_ ) );
	}


	/// @brief Stream access
	inline
	std::istream const &
	stream() const
	{
#if defined( USE_FILE_PROVIDER )
		if ( file_provider_stream->good() ) {
			return *file_provider_stream;
		}
#endif
		return ( zip_stream_p_
			? static_cast< std::istream const & >( *zip_stream_p_ )
			: static_cast< std::istream const & >( if_stream_ ) );
	}


	/// @brief Stream access
	inline
	std::istream &
	stream()
	{
#if defined( USE_FILE_PROVIDER )
		if ( file_provider_stream->good() ) {
			return *file_provider_stream;
		}
#endif
		return ( zip_stream_p_
			? static_cast< std::istream & >( *zip_stream_p_ )
			: static_cast< std::istream & >( if_stream_ ) );
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


	/// @brief Get the number of characters read by the last unformatted read
	inline
	std::streamsize
	gcount() const
	{
		return stream().gcount();
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


private: // Properties: internal zip stream predicates


	/// @brief Is stream attached to a gzip file?
	inline
	bool
	is_gzip() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->is_gzip() : false );
	}


	/// @brief CRC of gzip file valid?
	inline
	bool
	check_crc() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->check_crc() : false );
	}


	/// @brief CRC of the uncompressed data (see zipstream documentation)
	inline
	long
	get_crc() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->get_crc() : 0 );
	}


	/// @brief Compressed data size
	inline
	long
	get_in_size() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->get_in_size() : 0 );
	}


	/// @brief Uncompressed data size
	inline
	long
	get_out_size() const
	{
		return ( zip_stream_p_ ? zip_stream_p_->get_out_size() : 0 );
	}

private:

	/// @brief Helper function for opening files with alternative search paths
	void
	open_ifstream(
		std::string const & name,
		std::ios_base::openmode open_mode);

public:

	static
	void
	set_alternative_search_paths(
		vector1< std::string > alternative_search_paths
	){
		alternative_search_paths_ = alternative_search_paths;
	}


	static
	vector1< std::string >
	get_alternative_search_paths() {
		return alternative_search_paths_;
	}


private: // Fields


	/// @brief Compression state
	Compression compression_;

	/// @brief File stream
	std::ifstream if_stream_;

	/// @brief File name
	std::string filename_;

	/// @brief Zip file stream pointer (owning)
	zlib_stream::zip_istream * zip_stream_p_;

	/// @brief Alternative search paths
	/// This initialized by the option system -in:path:path Notice that
	/// izstream cannot access the option system (because the utility
	/// library comes before the basic library), so setting the
	/// alternate search paths is it the responsibility of core::init::init()
	static vector1< std::string > alternative_search_paths_;
#if defined( USE_FILE_PROVIDER )
	std::istream *file_provider_stream;
	std::stringstream bad_stream;
#endif

}; // izstream


// Non-member izstream functions


/// @brief std::getline( std::istream, std::string ) wrapper
inline
izstream &
getline( izstream & stream, std::string & line )
{
	std::getline( stream(), line );
	return stream;
}


/// @brief std::getline( std::istream, std::string, char ) wrapper
inline
izstream &
getline( izstream & stream, std::string & line, char const delim )
{
	std::getline( stream(), line, delim );
	return stream;
}


/// @brief Skip rest of line and line terminator (manipulator)
inline
std::istream &
skip( izstream & stream )
{
	return stream().ignore( std::numeric_limits< std::streamsize >::max(), '\n' );
}


/// @brief Skip rest of line and line terminator (manipulator)
inline
std::istream &
skip( std::istream & stream )
{
	return stream.ignore( std::numeric_limits< std::streamsize >::max(), '\n' );
}


} // namespace io
} // namespace utility


#endif // INCLUDED_utility_io_izstream_HH
