// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/crashreport.hh
/// @brief  Save crash reporting information to a file
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_crashreport_hh
#define INCLUDED_utility_crashreport_hh

#include <string>

#define do_crash(message) utility::save_crash_report( message, __FILE__, __LINE__ )

namespace utility {

/// @brief Set the crash handler for non-standard crashes (e.g. segfaults)
/// Should only be called once at startup.
void install_crash_handler();

/// @brief Set the application name that this was launched as.
/// Should only be called once on startup.
void set_application_name( char const * appname);

/// @brief Set the string representation of the options to be used in the crash report
/// Should only be called once on startup.
void set_options_string(std::string const & options);

/// @brief Save a crash report to the crash reporter file.
void save_crash_report(char const * message = "(none)", std::string const & file = "(none)", int const line = 0);

void save_crash_report(char const * message, std::string const & file, int const line, std::string const & traceback);

/// @brief Save a crash report to the crash reporter file.
inline
void save_crash_report(std::string const & message, std::string const & file = "(none)", int const line = 0) {
	save_crash_report( message.c_str(), file, line );
}

inline
void save_crash_report(std::string const & message, std::string const & file, int const line, std::string const & traceback) {
	save_crash_report( message.c_str(), file, line, traceback );
}

} // namespace utility

#endif // INCLUDED_utility_crashreport_HH
