// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/orstream.hh
/// @brief  Output stream wrapper abstract base class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_io_orstream_hh
#define INCLUDED_utility_io_orstream_hh


// Unit headers
#include <utility/io/orstream.fwd.hh>

// C++ headers
#include <iosfwd>
#include <string>
#include <iostream>


namespace utility {
namespace io {


/// @brief orstream: Output stream wrapper base class
class orstream
{


protected: // Types


	typedef  std::ostream & (*manipulator)( std::ostream & );


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~orstream()
	{}


protected: // Creation


	/// @brief Default constructor
	inline
	orstream()
	{}


private: // Creation


	/// @brief Copy constructor: Undefined
	orstream( orstream const & );


private: // Methods: assignment


	/// @brief Copy assignment: Undefined
	orstream &
	operator =( orstream const & );


public: // Methods: conversion


	/// @brief bool conversion
	virtual
	operator bool() const = 0;


	/// @brief Stream conversion
	virtual
	operator std::ostream const &() const = 0;


	/// @brief Stream conversion
	virtual
	operator std::ostream &() = 0;


public: // Methods: formatting


	/// @brief Stream output
	template< typename T >
	inline
	orstream &
	operator <<( T const & t )
	{
		stream() << t;
		return *this;
	}


	/// @brief Stream manipulator output
	virtual
	orstream &
	operator <<( manipulator m ) = 0;


public: // Methods: i/o


	/// @brief Flush the stream
	virtual
	orstream &
	flush() = 0;


	/// @brief Clear the stream
	virtual
	void
	clear() = 0;


	/// @brief Write a char
	virtual
	orstream &
	put( char const c ) = 0;


	/// @brief Write a string
	virtual
	orstream &
	write( char const * str, std::streamsize const count ) = 0;


	/// @brief Write a string
	virtual
	orstream &
	write( std::string const & str, std::streamsize const count ) = 0;


public: // Properties


	/// @brief Stream access
	virtual
	std::ostream const &
	operator ()() const = 0;


	/// @brief Stream access
	virtual
	std::ostream &
	operator ()() = 0;


	/// @brief Stream access
	virtual
	std::ostream const &
	stream() const = 0;


	/// @brief Stream access
	virtual
	std::ostream &
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


}; // orstream


} // namespace io
} // namespace utility


#endif // INCLUDED_utility_io_orstream_HH
