// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/crashreport.cc
/// @brief  Save crash reporting information to a file
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <utility/crash_report.hh>

#include <utility/version.hh>
#include <utility/backtrace.hh>
#include <utility/excn/Exceptions.hh>

#include <sstream>
#include <fstream>
#include <iostream>
#include <exception>
#include <typeinfo>
#ifndef __native_client__
#include <csignal>
#endif

// For the MPI version macros used below
#ifdef USEMPI
	#include <mpi.h>
#endif

// This dance is because it's the best way to turn a macro into a string
// see http://www.decompile.com/cpp/faq/file_and_line_error_string.htm

#define STRINGIFY_IMPL(x) #x
#define STRINGIFY(x) STRINGIFY_IMPL(x)

//////////////////////////////////
// Compilation/system settings
// These come primarily from https://sourceforge.net/p/predef/wiki/Home/

// Complier & version

#if defined(__xlC__)
static std::string const COMPILER("IBM XLC version " STRINGIFY(__xlc__));
#elif defined(_MSC_VER)
static std::string const COMPILER("Microsoft Visual C++ version " STRINGIFY(_MSC_VER));
#elif defined(__MINGW32__)
//need stdlib.h for the version macros
#include <stdlib.h>
static std::string const COMPILER("MinGW version " STRINGIFY(__MINGW32_MAJOR_VERSION) "." STRINGIFY(__MINGW32_MINOR_VERSION));
#elif defined(__ICC) || defined(__INTEL_COMPILER)
static std::string const COMPILER("Intel compiler version " STRINGIFY(__INTEL_COMPILER));
#elif defined(__clang__)
#ifdef __clang_version__
		static std::string const COMPILER("Clang version " STRINGIFY(__clang_version__));
#else
static std::string const COMPILER("Clang version (unknown)");
#endif
#elif defined(__GNUC__)
// GCC comes at the very end, as other compilers often define __GNUC__ if they're GCC compatible
#ifdef __VERSION__
		static std::string const COMPILER("GCC version " STRINGIFY(__VERSION__) );
#else
static std::string const COMPILER("GCC version (unknown)");
#endif
#else
	static std::string const COMPILER("Unknown compiler")
#endif

// C++ standard library details

#if defined(__INTEL_CXXLIB_ICC)
static std::string const STDLIB("Intel C++ Run-Time Libraries");
#elif defined(_LIBCPP_VERSION)
static std::string const STDLIB("libc++ version " STRINGIFY(_LIBCPP_VERSION));
#elif defined(__GLIBCXX__)
static std::string const STDLIB("libstdc++ version " STRINGIFY(__GLIBCXX__));
#elif defined(_CPPLIB_VER)
static std::string const STDLIB("Dinkumware libraries version " STRINGIFY(_CPPLIB_VER));
#else
	static std::string const STDLIB("Unknown C++ standard library");
#endif

// OS details
#if defined(__ANDROID__)
static std::string const OS_VER("Android");
#elif defined(__bg__)
static std::string const OS_VER("BlueGene");
#elif defined(__CYGWIN__)
static std::string const OS_VER("Cygwin");
#elif defined(_WIN32)
static std::string const OS_VER("Microsoft Windows");
#elif defined(__APPLE__) && defined(__MACH__)
#include <TargetConditionals.h>
#if TARGET_OS_IPHONE
static std::string const OS_VER("Apple iOS");
#else
		static std::string const OS_VER("Apple Mac OS X");
#endif
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__bsdi__) || defined(__DragonFly__)
static std::string const OS_VER("BSD");
#elif defined(__gnu_linux__)
static std::string const OS_VER("GNU/Linux");
#elif defined(__linux__)
static std::string const OS_VER("Linux");
#elif defined(__unix__) || defined(__unix)
static std::string const OS_VER("Unix, not otherwise specified");
#else
	static std::string const OS_VER("Unknown Operating System");
#endif

// Mode
// release & debug are the two big ones - the others aren't important for a crashback sense

#ifdef NDEBUG
	static std::string const MODE("Release");
#else
static std::string const MODE("Debug");
#endif


// Extras
// Have to do in a function because we can't build it all in one shot
std::string
get_extras() {
	std::string extras("");
#ifdef ANDROID
	extras += "android_arm ";
#endif
#ifdef LINK_APBS_LIB
	extras += "apbs ";
#endif
#ifdef BOINC
	extras += "boinc ";
#endif
#ifdef USEBOOSTMPI
	extras += "boost_mpi ";
#endif
#ifdef USE_BOOST_THREAD
	extras += "boost_thread ";
#endif
#ifdef USECUDA
	extras += "cuda ";
#endif
#ifdef MULTI_THREADED
	extras += "cxx11thread "; // Though also defined for omp
#endif
#ifdef CXX14
	extras += "cxx14 ";
#endif
#ifdef CXX17
	extras += "cxx17 ";
#endif
#ifdef CXX20
	extras += "cxx20 ";
#endif
#ifdef GL_GRAPHICS
	extras += "graphics ";
#endif
#ifdef USEHDF5
	extras += "hdf5 ";
#endif
#ifdef USEMPI
	extras += "mpi ";
	// Add MPI version information
#if defined(OMPI_MAJOR_VERSION)
		extras += "(OpenMPI " STRINGIFY(OMPI_MAJOR_VERSION) "." STRINGIFY(OMPI_MINOR_VERSION) "." STRINGIFY(OMPI_RELEASE_VERSION) ") ";
#elif defined(MVAPICH2_VERSION)
		extras += "(MVAPICH2 "  MVAPICH2_VERSION ")";
#elif defined(MPICH2_VERSION)
		extras += "(MPICH2 " MPICH2_VERSION ")";
#elif defined(MPICH_VERSION)
		extras += "(MPICH " MPICH_VERSION ")";
#else
		extras += "(MPI other)"
#endif
#endif
#ifdef USEMYSQL
	extras += "mysql ";
#endif
#ifdef USEOPENCL
	extras += "opencl ";
#endif
#ifdef USE_OPENMP
	extras += "omp ";
#endif
#ifdef USEPOSTGRES
	extras += "postgres ";
#endif
#ifdef PTR_BOOST
	extras += "ptr_boost ";
#endif
#ifdef WITH_PYTHON
	extras += "python ";
#endif
#ifdef SERIALIZATION
	extras += "serialization "; // also for zeromq
#endif
#ifdef USE_TENSORFLOW_CPU
	extras += "tensorflow ";
#endif
#ifdef USE_TENSORFLOW_GPU
	extras += "tensorflow_gpu ";
#endif
#ifdef ZEROMQ
	extras += "zeromq ";
#endif
	if ( extras.empty() ) {
		return "default";
	} else {
		return extras;
	}
}


namespace utility {

// Various variables and parameters for the crash reporting.

static std::string const CRASH_FILE("ROSETTA_CRASH.log");

static std::string APPNAME("UNKNOWN APPLICATION");
static std::string OPTIONS("NO OPTIONS SPECIFIED");

static std::string const HEADER(
	"##############################################################################################################\n"
	"#\n"
	"# Rosetta crash log. \n" // " Please submit the contents of this file to http://test.rosettacommons.org/crash_reports\n"
	"#\n\n"
);

#ifndef __native_client__
void signal_handler(int signal) {

	// First thing, reset the signal handlers such that if we trigger them again, we abort
	// This makes sure we're not re-entrant here
	std::signal(SIGFPE, SIG_DFL);
	std::signal(SIGILL, SIG_DFL);
	std::signal(SIGSEGV, SIG_DFL);
#ifndef _WIN32
	std::signal(SIGBUS, SIG_DFL);
#endif

	// Technically, doing basically anything here is undefined behavior
	// That's not really much of a problem, as we're dying in a fiery crash anyway.

	std::string message, name;
	switch(signal) {
	case SIGFPE :
		name = "SIGFPE";
		message = "Fatal arithmetic error encountered";
		break;
	case SIGILL :
		name = "SIGILL";
		message = "Illegal CPU instruction encountered";
		break;
	case SIGSEGV :
		name = "SIGSEGV";
		message = "Segmentation Fault";
		break;
#ifndef _WIN32
	case SIGBUS :
		name = "SIGBUS";
		message = "Bus Error";
		break;
#endif
	default :
		name = "SIGUNK";
		message = "Unknown signal encountered";
		break;
	}

	save_crash_report( message, name, signal ); // Reuse line for signal number

	raise(signal); // Re-raise signal, to let the default error handler run.
}

// The original terminate handler - for chaining purposes;
std::terminate_handler old_terminate_handler_ = nullptr;

void terminate_handler() {
	try {
		// Attempt to get some introspection about any pending exception
		try {
			auto exptr( std::current_exception() );
			if ( exptr ) {
				std::rethrow_exception(exptr);
			} else {
				// Non-exception reason to call terminate.
				save_crash_report("Unknown terminate","Unknown terminate");
			}
		} catch( utility::excn::Exception const & e ) {
			// Handle uncaught Rosetta error.
			e.crash_log();
		} catch( std::exception const & e ) {
			// Handle uncaught stdlib error
			std::string const & name( typeid( e ).name() );
			save_crash_report(e.what(), name);
		} catch ( ... ) {
			// Rethrowing something that's not a standard exception ... punt
			save_crash_report("Uncaught non-std exception","Uncaught non-std exception");
		}

	} catch ( ... ) {
		// Do nothing - The try/catch here is to be robust to secondary exceptions which are thrown in crash handling
	}
	old_terminate_handler_();
}

#else // __native_client__
void signal_handler(int) {}
void terminate_handler() {}
#endif // __native_client__

void install_crash_handler() {
#ifndef NOCRASHREPORT
#ifndef __native_client__
	if ( old_terminate_handler_ == nullptr ) { // Don't crush the terminate handler if we're re-entrant
		old_terminate_handler_ = std::set_terminate( terminate_handler ); // To catch unhandled (stdlib) exceptions
	}
	// SIGABRT SIGINT and SIGTERM are installed by JD2
	// We don't bother here, as they aren't internal crash errors
	std::signal(SIGFPE, signal_handler);
	std::signal(SIGILL, signal_handler);
	std::signal(SIGSEGV, signal_handler);
#ifndef _WIN32
	// If your platform doesn't recognize SIGBUS, or results in an error on this call, add to the define guard above
	std::signal(SIGBUS, signal_handler);
#endif // _WIN32
#endif // __native_client__
#endif // NOCRASHREPORT
}

void set_application_name( char const * appname) {
	APPNAME = appname;
}
void set_options_string(std::string const & options) {
	OPTIONS = options;
}

void save_crash_report(char const * message, std::string const & file, int line) {
	save_crash_report(message, file, line, backtrace_string(1)); // Capture current backtrace.
}

#ifndef NOCRASHREPORT
void save_crash_report(char const * message, std::string const & file, int line, std::string const & traceback) {
	// We create a single string with all the contents here, as
	// we want to save it all in one go to avoid interleaving multiple simultaneous crash reports.
	std::stringstream crash_log;
	crash_log << HEADER;
	crash_log << "[START_CRASH_REPORT]\n";

	crash_log << "[ROSETTA_VERSION]: " << utility::Version::version() << '\n';
	crash_log << "[COMMIT_DATE]: " << utility::Version::date() << '\n';
	crash_log << "[APPLICATION]: " << APPNAME << '\n';
	crash_log << "[MODE]: " << MODE << '\n';
	crash_log << "[EXTRAS]: " << get_extras() << '\n';
	crash_log << "[OS]: " << OS_VER << '\n';
	crash_log << "[COMPILER]: " << COMPILER << '\n';
	crash_log << "[STDLIB]: " << STDLIB << '\n';
	crash_log << "[START_OPTIONS]\n";
	crash_log << OPTIONS << '\n';
	crash_log << "\n[END_OPTIONS]\n";
	crash_log << '\n';

	// Add additional infomation here, like:
	// option settings

	crash_log << "[START_BACKTRACE]: " << "RAW_LIBC" << "\n"; // Add actual details about backtrace format here
	crash_log << traceback;
	crash_log << "\n[END_BACKTRACE]\n\n";

	crash_log << "[FILE]: " << file << '\n';
	crash_log << "[LINE]: " << line << '\n';
	crash_log << "[START_MESSAGE]\n";
	crash_log << message << '\n';
	crash_log << "\n[END_MESSAGE]\n";

	crash_log << "[END_CRASH_REPORT]\n";
	crash_log << '\n';

	// TODO: Try to open in the home directory?

	// Don't need izstream, because we always want to open this unzipped
	// NOTE: I hope opening as `app` here means that interleaving behavior from multiple processes is sane
	// if not, we may have to drop to the plain C interface level, which has better gurantees
	std::fstream fileout( CRASH_FILE, std::fstream::out | std::fstream::app);
	if ( fileout ) {
		fileout << crash_log.str();
		fileout.flush();
	}
	fileout.close();

	//std::cerr << "\n\nAN INTERNAL ERROR HAS OCCURED. PLEASE SUBMIT THE CONTENTS OF " << CRASH_FILE << " TO http://test.rosettacommons.org/crash_reports TO REPORT.\n\n" << std::endl;
	std::cerr << "\n\nAN INTERNAL ERROR HAS OCCURED. PLEASE SEE THE CONTENTS OF " << CRASH_FILE << " FOR DETAILS.\n\n" << std::endl;
	// TODO: Have a script which does the submission for you (and then deletes)
}
#else

// If we're running NOCRASHREPORT, just print the backtrace.
void save_crash_report(char const * message, std::string const &, int, std::string const & traceback) {
	std::cerr << utility::CSI_Magenta(); // set color of cerr to magenta
	std::cerr << "BACKTRACE:\n";
	std::cerr << traceback;
	std::cerr << utility::CSI_Reset();
}
#endif

} // namespace utility
