// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/relax/RelaxScriptManager.cc
/// @brief A singleton class for managing relax scripts, to ensure that they are loaded once and only once from disk.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <protocols/relax/RelaxScriptManager.hh>


// Unit headers

// Project header
#include <core/types.hh>

// Utility headers
#include <utility/pointer/memory.hh>
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <fstream>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Construct tracer.
static basic::Tracer TR( "protocols.relax.RelaxScriptManager" );

namespace protocols {
namespace relax {


// RelaxScriptFileContents Public methods /////////////////////////////////////////////////////////////
/// @brief File contents constructor.
RelaxScriptFileContents::RelaxScriptFileContents( utility::vector1< std::string > const & file_lines_in ) :
	utility::pointer::ReferenceCount(),
	file_lines_( file_lines_in )
{}

/// @brief Destructor.
RelaxScriptFileContents::~RelaxScriptFileContents() {}

/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
RelaxScriptFileContentsOP
RelaxScriptFileContents::clone() const {
	return utility::pointer::make_shared< RelaxScriptFileContents >(*this);
}

// RelaxScriptManager Public methods /////////////////////////////////////////////////////////////
// Static constant data access

/// @brief Get a relax script.  Load it from disk if it has not already been loaded.
/// @details Threadsafe and lazily loaded.
RelaxScriptFileContents const &
RelaxScriptManager::get_relax_script(
	std::string const &filename
) const {
	boost::function< RelaxScriptFileContentsOP () > creator( boost::bind( &RelaxScriptManager::create_relax_script_instance, boost::cref( filename ) ) );
	return *( utility::thread::safely_check_map_for_key_and_insert_if_absent( creator, SAFELY_PASS_MUTEX( relax_script_mutex_ ), filename, filename_to_filecontents_map_ ) );
}

// RelaxScriptManager Private methods ////////////////////////////////////////////////////////////

/// @brief Empty constructor.
RelaxScriptManager::RelaxScriptManager() :
	SingletonBase< RelaxScriptManager >(),
	filename_to_filecontents_map_()
#ifdef MULTI_THREADED
	,
	relax_script_mutex_()
#endif //MULTI_THREADED
{}

/// @brief Create an instance of a RelaxScriptFileContents object, by owning pointer.
/// @details Needed for threadsafe creation.  Loads data from disk.  NOT for repeated calls!
/// @note Not intended for use outside of RelaxScriptManager.
RelaxScriptFileContentsOP
RelaxScriptManager::create_relax_script_instance(
	std::string const & filename
) {
	utility::vector1< std::string > filelines;
	std::ifstream infile( filename.c_str() );
	if ( !infile.good() ) {
		utility_exit_with_message( "[ERROR] Error opening relaxscript file '" + filename );
	}
	TR << "================== Reading script file: " << filename << " ==================" << std::endl;

	std::string line;
	while ( getline( infile, line ) ) {
		//store line
		filelines.push_back( line );
	}
	infile.close();

	for ( auto const & line : filelines ) {
		TR << line << std::endl;
	}

	return utility::pointer::make_shared< RelaxScriptFileContents >(filelines);
}


} //protocols
} //relax
