// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/excn/Exceptions.hh
/// @brief  common derived classes for thrown exceptions
/// @author Oliver Lange
/// @author Sergey Lyskov


#ifndef INCLUDED_utility_excn_Exceptions_HH
#define INCLUDED_utility_excn_Exceptions_HH


// Unit Headers
#include <utility/excn/Exceptions.fwd.hh>

// Package Headers
#include <utility>
#include <string>
#include <ostream>
#include <exception>

namespace utility {
namespace excn {

/*
generally:
include files will be found in
/<namespace>/Exceptions.hh
for specialized Exceptions e.g. a EXCN_InvalidFoldTree

all-purpose exceptions are all bundled together in this header.
if this gets to big we will have extra forward declarations in Exceptions.fwd.hh
*/

#define CREATE_EXCEPTION(type, ...) type(__FILE__, __LINE__, __VA_ARGS__)

class Exception : public std::exception
{
public:
	Exception() = delete;  // NO DEFAULT CONSTRUCTOR, do not remove this line!!! This is needed to ensure that appropriate message specifying origin of the exception is provided. This will allow identification of origin of exceptions without need to dive into debugger.

	//
	// intended usage is CREATE_EXCEPTION(Exception, "something went wrong: ...");
	// or if you need to avoid macro Exception(__FILE__, __LINE__, "something went wrong: ...")
	Exception(char const *file, int line, std::string const &msg);

	virtual ~Exception() {}


	void show( std::ostream& ) const;

	std::string msg() const { return msg_; };

	void msg(std::string const &m) { msg_ = m; };

	void add_msg( std::string const& str ) {
		msg_ = msg_ + '\n' + str;
	}

private:
	// disabling access to what in favor of msg()
	char const * what() const noexcept override {
		return msg_.c_str();
	}

private:
	std::string msg_;
};

inline std::ostream& operator << ( std::ostream& os, Exception const & excn ) {
	excn.show( os );
	return os;
}


class IOError : public Exception {
protected:
	using Exception::Exception;
};

/// @brief EXCN_BadInput, as an IO error, should only be used for bad *user* input.
/// Do not use for something which is just bad function input.
class BadInput : public IOError {
public:
	using IOError::IOError;

};

class FileNotFound : public IOError {
public:
	FileNotFound(char const * file, int line, std::string const& not_found_file_name ) : IOError(file, line, std::string("unable to open file ") + file ), file_( not_found_file_name ) {};
private:
	std::string file_;
};

class RangeError : public Exception {
public:
	using Exception::Exception;
};

class KeyError : public Exception {
public:

	using Exception::Exception;
};

class NullPointerError: public RangeError {
public:
	using RangeError::RangeError;
};

class RosettaScriptsOptionError : public Exception {
public:
	using Exception::Exception;
};

class JD2Failure : public Exception {
public:
	using Exception::Exception;
};

}
}

#endif
