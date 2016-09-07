// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/irstream.hh
/// @brief  Input stream wrapper abstract base class
/// @author Mike Tyka (M.Tyka@bristol.ac.uk) (based on orstream.hh)


#ifndef INCLUDED_utility_io_irstream_hh
#define INCLUDED_utility_io_irstream_hh


// Unit headers
#include <utility/io/irstream.fwd.hh>

// C++ headers
#include <iosfwd>
#include <string>
#include <iostream>


namespace utility {
namespace io {


/// @brief orstream: Input stream wrapper abstract base class
class irstream
{


protected: // Types


	typedef  std::istream & (*std_manipulator)( std::istream & );


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~irstream()
	= default;


protected: // Creation


	/// @brief Default constructor
	inline
	irstream()
	= default;


private: // Creation


	/// @brief Copy constructor: Undefined
	irstream( irstream const & );


private: // Methods: assignment


	/// @brief Copy assignment: Undefined
	irstream &
	operator =( irstream const & );


public: // Methods: conversion


	/// @brief bool conversion
	virtual
	operator bool() const = 0;


	/// @brief Stream conversion
	virtual
	operator std::istream const &() const = 0;


	/// @brief Stream conversion
	virtual
	operator std::istream &() = 0;


public: // Methods: formatting


	/// @brief Stream input
	template< typename T >
	inline
	irstream &
	operator >>( T & t )
	{
		stream() >> t;
		return *this;
	}


	/// @brief Stream manipulator input
	virtual
	std::istream &
	operator >>( std_manipulator m ) = 0;


public: // Methods: i/o


	/// @brief Clear the stream
	virtual
	void
	clear() = 0;


	/// @brief Seek to the beginning
	virtual
	void
	seek_beg() = 0;


	/// @brief Get the next character
	virtual
	int
	get() = 0;


	/// @brief Get the next character
	virtual
	irstream &
	get( char & c ) = 0;


	/// @brief Get the next specified number of characters
	virtual
	irstream &
	get( char * str, std::streamsize const count ) = 0;


	/// @brief Get the next specified number of characters
	virtual
	irstream &
	get( char * str, std::streamsize const count, char const delim ) = 0;


	/// @brief Get the next specified number of characters
	virtual
	irstream &
	get( std::string & str, std::streamsize const count ) = 0;


	/// @brief Get the next specified number of characters
	virtual
	irstream &
	get( std::string & str, std::streamsize const count, char const delim ) = 0;


	/// @brief Get the rest of the line
	virtual
	irstream &
	getline( char * line, std::streamsize const count ) = 0;


	/// @brief Get the rest of the line
	virtual
	irstream &
	getline( char * line, std::streamsize const count, char const delim ) = 0;

	/// @brief Get the rest of the line
	virtual
	irstream &
	getline( std::string & line ) = 0;


	/// @brief Get the rest of the line
	virtual
	irstream &
	getline( std::string & line, char const delim ) = 0;


	/// @brief Read the next specified number of characters
	virtual
	irstream &
	read( char * str, std::streamsize const count ) = 0;


	/// @brief Read the next specified number of characters
	virtual
	irstream &
	read( std::string & str, std::streamsize const count ) = 0;


	/// @brief Read the next available specified number of characters
	virtual
	std::streamsize
	readsome( char * str, std::streamsize const count ) = 0;


	/// @brief Read the next available specified number of characters
	virtual
	std::streamsize
	readsome( std::string & str, std::streamsize const count ) = 0;


	/// @brief Skip over the next character
	virtual
	irstream &
	ignore() = 0;


	/// @brief Skip over the next specified number of characters
	virtual
	irstream &
	ignore( std::streamsize const count ) = 0;


	/// @brief Skip over the next specified number of characters
	virtual
	irstream &
	ignore( std::streamsize const count, char const delim ) = 0;


	/// @brief Returns the next character without extracting it
	virtual
	int
	peek() = 0;


	/// @brief Put the last character read back into the stream
	virtual
	irstream &
	unget() = 0;


	/// @brief Put the last character read back into the stream and check that passed character is correct
	virtual
	irstream &
	putback( char c ) = 0;


public: // Properties


	/// @brief Stream access
	virtual
	std::istream const &
	operator ()() const = 0;


	/// @brief Stream access
	virtual
	std::istream &
	operator ()() = 0;


	/// @brief Stream access
	virtual
	std::istream const &
	stream() const = 0;


	/// @brief Stream access
	virtual
	std::istream &
	stream() = 0;


	/// @brief Pointer to the stream buffer
	virtual
	std::streambuf *
	rdbuf() const = 0;


public: // Properties: predicate


	/// @brief Good?
	virtual
	bool
	good() const = 0;


	/// @brief End of file?
	virtual
	bool
	eof() const = 0;


	/// @brief Fail?
	virtual
	bool
	fail() const = 0;


	/// @brief Bad?
	virtual
	bool
	bad() const = 0;


	/// @brief Compressed?
	virtual
	bool
	compressed() const = 0;


	/// @brief Uncompressed?
	virtual
	bool
	uncompressed() const = 0;


	/// @brief gzipped?
	virtual
	bool
	gzipped() const = 0;


}; // irstream


} // namespace io
} // namespace utility


#endif // INCLUDED_utility_io_irstream_HH
