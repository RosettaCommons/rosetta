// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/ocstream.hh
/// @brief  Output channel stream wrapper class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_io_ocstream_hh
#define INCLUDED_utility_io_ocstream_hh


// Unit headers
#include <utility/io/ocstream.fwd.hh>

// Package headers
#include <utility/io/orstream.hh>

// C++ headers
#include <ostream>


namespace utility {
namespace io {


/// @brief ocstream: Output channel stream wrapper class
class ocstream :
	public orstream
{


public: // Creation


	/// @brief Constructor
	inline
	ocstream( std::ostream & o_stream_a ) :
		o_stream_( o_stream_a )
	{}


	/// @brief Destructor
	inline

	~ocstream() override
	= default;


public: // Methods: conversion


	/// @brief bool conversion
	inline
	operator bool() const override
	{
		return !!o_stream_;
	}


	/// @brief Stream conversion
	inline
	operator std::ostream const &() const override
	{
		return o_stream_;
	}


	/// @brief Stream conversion
	inline
	operator std::ostream &() override
	{
		return o_stream_;
	}


public: // Methods: formatting


	/// @brief Stream output: override to preserve type of return value
	template< typename T >
	inline
	ocstream &
	operator <<( T const & t )
	{
		o_stream_ << t;
		return *this;
	}


	/// @brief Stream manipulator output
	inline
	ocstream &
	operator <<( manipulator m ) override
	{
		o_stream_ << m;
		return *this;
	}


public: // Methods: i/o


	/// @brief Write a char
	inline
	ocstream &
	put( char const c ) override
	{
		o_stream_.put( c );
		return *this;
	}


	/// @brief Write a string
	inline
	ocstream &
	write( char const * str, std::streamsize const count ) override
	{
		o_stream_.write( str, count );
		return *this;
	}


	/// @brief Write a string
	inline
	ocstream &
	write( std::string const & str, std::streamsize const count ) override
	{
		o_stream_.write( str.c_str(), count );
		return *this;
	}


	/// @brief Flush the stream
	inline
	ocstream &
	flush() override
	{
		o_stream_.flush();
		return *this;
	}


	/// @brief Clear the stream
	inline
	void
	clear() override
	{
		o_stream_.clear();
	}


public: // Properties


	/// @brief Stream access
	inline
	std::ostream const &
	operator ()() const override
	{
		return o_stream_;
	}


	/// @brief Stream access
	inline
	std::ostream &
	operator ()() override
	{
		return o_stream_;
	}


	/// @brief Stream access
	inline
	std::ostream const &
	stream() const override
	{
		return o_stream_;
	}


	/// @brief Stream access
	inline
	std::ostream &
	stream() override
	{
		return o_stream_;
	}


	/// @brief Pointer to the stream buffer
	inline
	std::streambuf *
	rdbuf() const override
	{
		return o_stream_.rdbuf();
	}


public: // Properties: predicate


	/// @brief Good?
	inline
	bool
	good() const override
	{
		return o_stream_.good();
	}


	/// @brief End of file?
	inline
	bool
	eof() const override
	{
		return o_stream_.eof();
	}


	/// @brief Fail?
	inline
	bool
	fail() const override
	{
		return o_stream_.fail();
	}


	/// @brief Bad?
	inline
	bool
	bad() const override
	{
		return o_stream_.bad();
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


	/// @brief Output stream reference
	std::ostream & o_stream_;


}; // ocstream


namespace oc { // Predefined ocstreams


/// @brief Wrapper around std::cout
extern ocstream cout;

/// @brief Wrapper around std::cerr
extern ocstream cerr;

/// @brief Wrapper around std::clog
extern ocstream clog;


} // namespace oc


} // namespace io
} // namespace utility


#endif // INCLUDED_utility_io_ocstream_HH
