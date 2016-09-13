// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/icstream.hh
/// @brief  Input channel stream wrapper class
/// @author Mike Tyka (M.Tyka@bristol.ac.uk) (based on ocstream.hh)


#ifndef INCLUDED_utility_io_icstream_hh
#define INCLUDED_utility_io_icstream_hh


// Unit headers
#include <utility/io/icstream.fwd.hh>

// Package headers
#include <utility/io/irstream.hh>

// C++ headers
#include <istream>


namespace utility {
namespace io {


/// @brief icstream: Input channel stream wrapper class
class icstream :
	public irstream
{


public: // Creation


	/// @brief Constructor
	inline
	icstream( std::istream & i_stream_a ) :
		i_stream_( i_stream_a )
	{}


	/// @brief Destructor
	inline

	~icstream() {}


public: // Methods: conversion


	/// @brief bool conversion
	inline
	operator bool() const override
	{
		//return i_stream_.good();
		// proper behavior is actually ( ! fail() )
		return !!i_stream_;
	}


	/// @brief Stream conversion
	inline
	operator std::istream const &() const override
	{
		return i_stream_;
	}


	/// @brief Stream conversion
	inline
	operator std::istream &() override
	{
		return i_stream_;
	}


public: // Methods: formatting


	/// @brief Stream input: Override to preserve type of return value
	template< typename T >
	inline
	icstream &
	operator >>( T & t )
	{
		i_stream_ >> t;
		return *this;
	}


	/// @brief Stream manipulator input
	inline
	std::istream &
	operator >>( std_manipulator m ) override
	{
		return m( *this );
	}


public: // Methods: i/o


	/// @brief Seek to the beginning
	inline
	void
	seek_beg() override
	{
		stream().clear();
		stream().seekg( std::ios_base::beg );
		stream().clear();
	}


	/// @brief Get the next character
	inline
	int
	get() override
	{
		return stream().get();
	}


	/// @brief Get the next character
	inline
	icstream &
	get( char & c ) override
	{
		stream().get( c );
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	icstream &
	get( char * str, std::streamsize const count ) override
	{
		stream().get( str, count );
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	icstream &
	get( char * str, std::streamsize const count, char const delim ) override
	{
		stream().get( str, count, delim );
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	icstream &
	get( std::string & str, std::streamsize const count ) override
	{
		auto * cp = new char[ count ];
		stream().get( cp, count );
		str = cp;
		delete[] cp;
		return *this;
	}


	/// @brief Get the next specified number of characters
	inline
	icstream &
	get( std::string & str, std::streamsize const count, char const delim ) override
	{
		auto * cp = new char[ count ];
		stream().get( cp, count, delim );
		str = cp;
		delete[] cp;
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	icstream &
	getline( char * line, std::streamsize const count ) override
	{
		stream().getline( line, count );
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	icstream &
	getline( char * line, std::streamsize const count, char const delim ) override
	{
		stream().getline( line, count, delim );
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	icstream &
	getline( std::string & line ) override
	{
		std::getline( stream(), line );
		return *this;
	}


	/// @brief Get the rest of the line
	inline
	icstream &
	getline( std::string & line, char const delim ) override
	{
		std::getline( stream(), line, delim );
		return *this;
	}


	/// @brief Read the next specified number of characters
	inline
	icstream &
	read( char * str, std::streamsize const count ) override
	{
		stream().read( str, count );
		return *this;
	}


	/// @brief Read the next specified number of characters
	inline
	icstream &
	read( std::string & str, std::streamsize const count ) override
	{
		auto * cp = new char[ count ];
		stream().read( cp, count );
		str = cp;
		delete[] cp;
		return *this;
	}


	/// @brief Read the next available specified number of characters
	inline
	std::streamsize
	readsome( char * str, std::streamsize const count ) override
	{
		return stream().readsome( str, count );
	}


	/// @brief Read the next available specified number of characters
	inline
	std::streamsize
	readsome( std::string & str, std::streamsize const count ) override
	{
		auto * cp = new char[ count ];
		std::streamsize const n_chars = stream().readsome( cp, count );
		str = cp;
		delete[] cp;
		return n_chars;
	}


	/// @brief Skip over the next character
	inline
	icstream &
	ignore() override
	{
		stream().ignore();
		return *this;
	}


	/// @brief Skip over the next specified number of characters
	inline
	icstream &
	ignore( std::streamsize const count ) override
	{
		stream().ignore( count );
		return *this;
	}


	/// @brief Skip over the next specified number of characters
	inline
	icstream &
	ignore( std::streamsize const count, char const delim ) override
	{
		stream().ignore( count, delim );
		return *this;
	}


	/// @brief Returns the next character without extracting it
	inline
	int
	peek() override
	{
		return stream().peek();
	}


	/// @brief Put the last character read back into the stream
	inline
	icstream &
	unget() override
	{
		stream().unget();
		return *this;
	}


	/// @brief Put the last character read back into the stream and check that passed character is correct
	inline
	icstream &
	putback( char c ) override
	{
		stream().putback( c );
		return *this;
	}


	/// @brief Clear the stream
	inline
	void
	clear() override
	{
		i_stream_.clear();
	}


public: // Properties


	/// @brief Stream access
	inline
	std::istream const &
	operator ()() const override
	{
		return i_stream_;
	}


	/// @brief Stream access
	inline
	std::istream &
	operator ()() override
	{
		return i_stream_;
	}


	/// @brief Stream access
	inline
	std::istream const &
	stream() const override
	{
		return i_stream_;
	}


	/// @brief Stream access
	inline
	std::istream &
	stream() override
	{
		return i_stream_;
	}


	/// @brief Pointer to the stream buffer
	inline
	std::streambuf *
	rdbuf() const override
	{
		return i_stream_.rdbuf();
	}


public: // Properties: predicate


	/// @brief Good?
	inline
	bool
	good() const override
	{
		return i_stream_.good();
	}


	/// @brief End of file?
	inline
	bool
	eof() const override
	{
		return i_stream_.eof();
	}


	/// @brief Fail?
	inline
	bool
	fail() const override
	{
		return i_stream_.fail();
	}


	/// @brief Bad?
	inline
	bool
	bad() const override
	{
		return i_stream_.bad();
	}


	/// @brief Compressed?
	inline
	bool
	compressed() const override
	{
		return false;
	}


	/// @brief Uncompressed?
	inline
	bool
	uncompressed() const override
	{
		return true;
	}


	/// @brief gzipped?
	inline
	bool
	gzipped() const override
	{
		return false;
	}


private: // Fields


	/// @brief Input stream reference
	std::istream & i_stream_;


}; // icstream


namespace ic { // Predefined icstreams


/// @brief Wrapper around std::cin
extern icstream cin;


} // namespace ic


} // namespace io
} // namespace utility


#endif // INCLUDED_utility_io_icstream_HH
