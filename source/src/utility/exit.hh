// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/exit.hh
/// @brief  Program exit functions and macros
/// @author David Kim (dekim@u.washington.edu)
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note The point of these is:
///  @li  Show the file and line number where the exit originates
///  @li  Optionally show a message about why the exit occurred
///  @li  Provide a core dump on Linux/UNIX so a post-mortem backtrace
///       can be performed
///  @li  Provide macro functions to add the file and line for you
/// @note Break on utility::exit when debugging to allow you to
///       get a backtrace from the point of exit


#ifndef INCLUDED_utility_exit_hh
#define INCLUDED_utility_exit_hh


// C++ headers
#include <string>


/// @brief Macro function wrappers for utility::exit
///
/// @note Convenience macros that fills in the file and line
/// @note These have to be macros to get the file and line from the point of call
/// @note Don't use variadic macros to reduce the "overloads": They aren't standard C++

/// @brief Exit with file + line
#define utility_exit() utility::exit( __FILE__, __LINE__ )


/// @brief Exit with file + line + message
///
/// @note The m argument is a message string
#define utility_exit_with_message(m) utility::exit( __FILE__, __LINE__, m )


/// @brief Exit with file + line + status
///
/// @note The s argument is a status value
#define utility_exit_with_status(s) utility::exit( __FILE__, __LINE__, s )


/// @brief Exit with file + line + message + status
///
/// @note The m argument is a message string
/// @note The s argument is a status value
#define utility_exit_with_message_status(m,s) utility::exit( __FILE__, __LINE__, m, s )

/// @brief Assert that the condition holds. Evaluated for both debug and release builds
#define runtime_assert(_Expression) if ( !(_Expression) ) utility::exit(__FILE__, __LINE__, #_Expression)

/// @brief Assert that the condition holds. Evaluated for both debug and release builds
#define runtime_assert_msg(_Expression, msg) \
	if ( !(_Expression) ) utility::exit(__FILE__, __LINE__, #_Expression " MSG:" msg )

/// @brief Assert that the condition holds. Evaluated for both debug and release builds
// Does the same thing as runtime_assert_msg but allows C++ strings to be used for the message.
#define runtime_assert_string_msg(_Expression, msg) \
	if ( !(_Expression) ) utility::exit(__FILE__, __LINE__, msg )

namespace utility {


#ifdef __GNUC__
#  define NORETURN __attribute__ ((noreturn))
#elif __clang__
#  define NORETURN __attribute__ ((noreturn))
#else
#  define NORETURN
#endif

/// @brief Exit with file + line + message + optional status
void
exit(
	std::string const & file,
	int const line,
	std::string const & message,
	int const status = 1
) NORETURN;

/// @brief Conditional Exit with file + line + message + optional status. WIll exit if the condition is not met!
int
cond_exit(
	bool condition,
	std::string const & file,
	int const line,
	std::string const & message,
	int const status = 1
);


/// @brief Exit with file + line + optional status
inline
void
exit(
	std::string const & file,
	int const line,
	int const status = 1
) NORETURN;


/// @brief Exit with file + line + optional status
inline
void
exit(
	std::string const & file,
	int const line,
	int const status
)
{
	utility::exit( file, line, std::string(), status );
}


/// @brief Exit with file + line + status
///
/// @note  Deprecated: For backwards compatibility with earlier version
inline
void
exit(
	int const status,
	std::string const & file,
	int const line
)
{
	utility::exit( file, line, std::string(), status );
}

typedef void (* UtilityExitCallBack)(void);

/// @brief Set call back funtion that will be called on utility::exit.
///        Use this function to overload default behavior of sys.exit to more appropriate to your application
///        Defaut value for callback function is 0, whicth mean no sys exit is called.
void set_main_exit_callback( UtilityExitCallBack = 0 );

/// @brief Add additional callback function that will be called *before* standard exit(…) is executed.
///        [Note: do not confuse this function with 'set_main_exit_callback' which is replacing the end behavior of exit(…)]
void add_exit_callback( UtilityExitCallBack );

/// @brief Remove additional callback function that was previously added by using add_exit_callback.
void remove_exit_callback( UtilityExitCallBack );

} // namespace utility


#endif // INCLUDED_utility_exit_HH
