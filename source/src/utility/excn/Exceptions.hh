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
	// Avoid the latter, though, as it circumvents the crash log handling
	Exception(char const *file, int line, std::string const &msg);

	virtual ~Exception() {}

	/// @brief Present this exception to the user.
	/// Will invoke crash log reporting, if applicable
	virtual void display() const;

	/// @brief Invoke a crash log for throwing this exception.
	/// If your exception is one which is a "non-error" exception,
	/// override this function to do nothing.
	virtual void crash_log() const;

	void show( std::ostream& ) const;

	/// @brief Will return a formatted message (with file/line information)
	std::string msg() const;

	/// @brief Will return the raw message (without file/line information)
	std::string raw_msg() const { return msg_; };

	/// @brief Will set the *raw* message.
	void msg(std::string const &m) { msg_ = m; };

	void add_msg( std::string const& str ) {
		msg_ = msg_ + '\n' + str;
	}

	void prepend_to_msg( std::string const& str ) {
		msg_ = str + '\n' + msg_;
	}

	std::string const & file() { return file_; }
	int line() { return line_; }
	std::string const & traceback() { return traceback_; }

private:
	// disabling access to what in favor of msg()
	char const * what() const noexcept override {
		// As msg returns a temp string, we need to keep a reference to it
		// in order to return a pointer to the contents.
		whatstring_ = msg();
		return whatstring_.c_str();
	}

private:
	std::string msg_;
	std::string file_;
	int line_;
	std::string traceback_;
	mutable std::string whatstring_; // only used by what()
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

	// Should be a user error, but I'm not yet confident the error messages are all that revealing yet -- wait until they're better to turn off crash reporting - RM
	// void crash_log() const override {}
};

class FileNotFound : public IOError {
public:
	FileNotFound(char const * file, int line, std::string const& not_found_file_name ) : IOError(file, line, std::string("unable to open file ") + file ), file_( not_found_file_name ) {};
	void crash_log() const override {} // User error, not crash, & fix should be clear from error message.
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

class JD2Failure : public Exception {
public:
	using Exception::Exception;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//
class UserCorrectableIssue: public Exception {
public:
	using Exception::Exception;

	void display() const override;

	void crash_log() const override;
};

class RosettaScriptsOptionError : public UserCorrectableIssue {
public:
	using UserCorrectableIssue::UserCorrectableIssue;
};

}
}

#endif
