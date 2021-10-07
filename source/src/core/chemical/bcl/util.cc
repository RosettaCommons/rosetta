// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/bcl/util.cc
/// @brief Utilities for interacting with the BCL.
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Benjamin P. Brown (benjamin.p.brown17@gmail.com) (edited 10/2021)

#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__)
#include <mach-o/dyld.h> // for _NSGetExecutablePath
#endif

// core includes
#include <core/chemical/bcl/util.hh>
#include <core/types.hh>

// utility includes
#include <utility/exit.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/OptionCollection.hh>

// basic includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/app/bcl_app_apps.h>
#include <bcl/include/command/bcl_command_command.h>
#include <bcl/include/command/bcl_command_app_default_flags.h>
#include <bcl/include/command/bcl_command_command_state.h>
#include <bcl/include/command/bcl_command_parameter_check_file_in_search_path.h>
#include <bcl/include/model/bcl_model.h>
#include <bcl/include/random/bcl_random_uniform_distribution.h>
#include <bcl/include/io/bcl_io_directory.h>
#include <bcl/include/storage/bcl_storage_vector.h>
#include <bcl/include/util/bcl_util_logger_interface.h>
#include <bcl/include/util/bcl_util_message.h>
#endif

// C++
#include <iostream>

namespace core {
namespace chemical {
namespace bcl {

static const std::string BCL_VERSION = "BCL v4.1.0";

static basic::Tracer TR("core.chemical.bcl.util");
// The following tracer is not meant for actual output, but is instead intended to be used
// to control the BCL output levels via the Rosetta tracer controls.
static basic::Tracer TR_BCL("BCL");

/// @brief Utility function to print BCL usage details.
void require_bcl() {
#ifndef USEBCL
	TR.Fatal << std::endl;
	TR.Fatal << "ATTENTION: This functionality requires the BCL. (An extras=bcl build.)" << std::endl;
	TR.Fatal << "ATTENTION: Please re-run the protocol with a version of Rosetta compiled against the BCL." << std::endl;
	TR.Fatal << "ATTENTION: The use of BCL version " << BCL_VERSION << " is highly recommended." << std::endl;
	TR.Fatal << "ATTENTION: The BCL can be obtained from https:/www.meilerlab.org" << std::endl;
	TR.Fatal << "ATTENTION: See https://www.rosettacommons.org/docs/latest/BCL.html for more details." << std::endl;
	TR.Fatal << std::endl;
	utility_exit_with_message("Use of extras=bcl build required.");
#endif
#ifdef MULTI_THREADED
	TR.Fatal << std::endl;
	utility_exit_with_message("Use of extras=bcl is not currently threadsafe");
	TR.Fatal << std::endl;
#endif
}

/// @brief Initialize core components of the bcl apps/apps.cpp
/// @details e.g. GetExecutablePath, etc.
void initialize_bcl_main()
{
#ifdef USEBCL
	// update whether this is static initialization time
	::bcl::command::CommandState::IsInStaticInitialization() = false;

	// set bool that we're done with Main's command line parsing module
	::bcl::command::CommandState::GetInMainCommandLineParsing() = false;

	//! Set some important paths; reference the Rosetta database location

	// set BCL executable path
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string bcl_exe_dir( std::string( option[ in::path::bcl ](1).name()) + "../");
	TR << "Externals directory: " + bcl_exe_dir << std::endl;
	::bcl::app::Apps::GetExecutablePath() = bcl_exe_dir;

	// set path to the rotamer library
	std::string bcl_rotamer_lib( std::string( option[ in::path::bcl ](1).name()) + "rotamer_library/");
	TR << "BCL rotamer library directory: " + bcl_rotamer_lib << std::endl;
	::bcl::command::ParameterCheckFileInSearchPath("rotamer_library", bcl_rotamer_lib, ::bcl::io::Directory::e_Dir);

	// set path to model directory
	std::string bcl_model_dir( std::string( option[ in::path::bcl ](1).name()) + "model/");
	TR << "BCL model directory: " + bcl_model_dir << std::endl;
	::bcl::command::ParameterCheckFileInSearchPath("model", bcl_model_dir, ::bcl::io::Directory::e_Dir);

	// set BCL command-line required arguments;
	// get the arguments as a vector
	const ::bcl::storage::Vector< std::string> arguments
		(
		::bcl::storage::Vector< std::string>::Create
		(
		"-message_level", "Standard",
		"-model_path", bcl_model_dir,
		"-random_seed"
		)
	);

	// print model path
	TR << "BCL model directory passed to commandline: " + ::bcl::model::Model::GetModelPathFlag()->GetFirstParameter()->GetValue() << std::endl;

	// create the command state object; parse all arguments except the 1st, which is the bcl executable command
	::bcl::command::CommandState &cmd_state( ::bcl::command::CommandState::GetGlobalCommandState());

	// create a string stream to catch all error output
	std::ostringstream errors;

	// default BCL command line flags
	::bcl::command::Command def_cmd;
	::bcl::command::GetAppDefaultFlags().AddDefaultCommandlineFlags( def_cmd);
	def_cmd.SetFlags( cmd_state, ::bcl::util::GetLogger());
	::bcl::util::GetLogger() << errors.str() << std::endl;

#endif
}

/// @brief Initialize the BCL random number generator.
/// @details Note that seed is an int to match the seed generated in core/init.cc
void initialize_bcl_random( int const seed )
{
	TR.Trace << "Passing " << seed << " to BCL RNG, if enabled" << std::endl; // Silence compiler warning when USEBCL isn't defined.
#ifdef USEBCL
	::bcl::random::GetGlobalRandom().SetSeed( seed);
	TR << "Initializing BCL random number generator with seed: " << ::bcl::random::GetGlobalRandom().GetSeed() << std::endl;
#endif
	// No-op if we're not running BCL.
}

void initialize_bcl_tracers()
{
#ifdef USEBCL
	::bcl::util::Message::MessageLevel level( ::bcl::util::Message::e_Standard );
	// Have to go chattiest to most laconic, as visibility is nested.
	if ( TR_BCL.visible(basic::t_trace) ) {
		level = ::bcl::util::Message::e_Debug;
	} else if ( TR_BCL.visible(basic::t_debug) ) {
		level = ::bcl::util::Message::e_Verbose;
	} else if ( TR_BCL.visible(basic::t_info) ) {
		level = ::bcl::util::Message::e_Standard;
	} else if ( TR_BCL.visible(basic::t_warning) ) {
		level = ::bcl::util::Message::e_Critical;
	} else if ( TR_BCL.visible(basic::t_error) ) {
		level = ::bcl::util::Message::e_Silent;
	} else { // t_fatal or muted
		level = ::bcl::util::Message::e_Error;
	}
	::bcl::util::GetMessenger().SetMessageLevel(level);

	TR << "BCL message level set at " << ::bcl::util::Message::GetLevelString( ::bcl::util::GetMessenger().GetCurrentMessageLevel() ) << std::endl;
#endif
	// No-op if we're not running BCL
}

/// @brief Locate the BCL executable path within /main/source/externals/bcl
/// @details Required that we find the BCL submodule executable path so that
/// we can set paths to the rotamer library and model directories; this is based
/// obviously on the locate_database() function in init.cc with some extra logic
/// and minor additional differences
void locate_bcl(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// blatantly copied and modified from locate_database() in init.cc
#ifndef __native_client__
	if ( !option[ in::path::bcl ].user() ) {
		std::string bcl_path;

		char path[1024];
		uint32_t path_size = sizeof( path );
		TR.Debug << path_size << std::endl;
		std::string path_string;


#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__)
		uint32_t result = _NSGetExecutablePath(path, &path_size );

		// _NSGetExecutablePath returns 0 if the path was successfully copied.
		// Copies a null-terminated string
		if ( result == 0 ) {
			path_string = std::string( path );
		}
#endif

#if defined(linux) || defined(__linux__) || defined(__linux)
		path_size = readlink( "/proc/self/exe", path, path_size ); // This works on plain Linux, but not FreeBSD or Solaris.

		// readlink returns -1 on error, otherwise number of bytes written.
		// Does not append null to string
		if ( path_size > 0 ) {
			path_string = std::string( path, path_size );
		}
#endif
		// Rosetta executable path
		TR << "Resolved executable path: " << path_string << std::endl;

		// Attempt to resolve from source/external/bcl if 'source' is in path or external/bcl if not.
		// Note that here we must do source/external rather than source/.. as in locate_database()
		if ( path_string.length() > 0 ) {
			// get the substring "source/" from the executable path
			Size found = std::string::npos;
			if ( (found = path_string.find("source/")) && (found != std::string::npos) ) {
				// this should put us at main/
				std::string rosetta_exe_dir = path_string.substr(0,found);
				bcl_path = rosetta_exe_dir + "source/external/bcl/";
				TR << "Looking for bcl based on location of executable: " << bcl_path << std::endl;
			} else if ( (found = path_string.rfind("/")) && (found != std::string::npos) ) {
				// this should put us at main/source/
				std::string rosetta_exe_dir = path_string.substr(0,found);
				bcl_path = rosetta_exe_dir + "external/bcl/";
				TR << "Looking for bcl based on location of executable: " << bcl_path << std::endl;
			}
		} else if ( option[ in::path::database ]().size() ) {
			// Attempt to resolve from database (not original with locate_database())
			std::string database_path( std::string( option[ in::path::database ](1).name()));
			bcl_path = database_path + "../source/external/bcl";
		} else {
			TR << "Could not determine location of executable or database." << std::endl;
		}

		// set the final bcl path
		if ( bcl_path.size() > 0 ) {
			TR << "Found BCL submodule path located at: " + bcl_path << std::endl;
			option[ in::path::bcl ].value( bcl_path );
		} else {
			TR.Fatal << "Could not find the BCL submodule." << std::endl;
			utility_exit_with_message("Double check that a valid BCL submodule has been cloned with Rosetta");
		}
	}
#endif
}



} // namespace bcl
} // namespace chemical
} // namespace core
