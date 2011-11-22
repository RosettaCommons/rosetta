// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#if (defined _WIN32) && (!defined PYROSETTA)
	#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

// Unit headers
#include <basic/options/util.hh>

#include <basic/options/option.hh>


// Utility headers
//#include <utility/options/OptionCollection.hh>
//#include <utility/options/keys/OptionKeys.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>




namespace basic {
namespace options {

/// kind of like the old -s -- just one file
std::string
start_file()
{
	utility::vector1< std::string > const files( start_files() );
	runtime_assert( files.size() == 1 );
	return files[1];
}


/// @details Retrieve all files specified using -l or -s
/// NOTE:: will die if neither option is used!
///
utility::vector1< std::string >
start_files()
{
	using namespace OptionKeys::in::file;

	utility::vector1< std::string > filenames;

	// -l files // this section should be deleted and this option removed
	if ( option[ l ].user() ) {
		utility::vector1< std::string > list_files( option[ l ]() );
		for ( size_t i=1; i<= list_files.size(); ++i ) {
			utility::io::izstream data( list_files[i].c_str() );
			if ( !data.good() ) {
				utility_exit_with_message("Unable to open list file: "+list_files[i]);
			}
			while ( data.good() ) {
				std::string name;
				data >> name;
				if ( data.good() ) filenames.push_back( name );
			}
			data.close();
		}
	}

	// -list files
	if ( option[ list ].user() ) {
		utility::vector1< std::string > list_files( option[ list ]() );
		for ( size_t i=1; i<= list_files.size(); ++i ) {
			utility::io::izstream data( list_files[i].c_str() );
			if ( !data.good() ) {
				utility_exit_with_message("Unable to open list file: "+list_files[i]);
			}
			while ( data.good() ) {
				std::string name;
				data.getline(name);
				if ( data.good() ) filenames.push_back( name );
			}
			data.close();
		}
	}


	// -s files
	if ( option[ s ].user() ) {
		utility::vector1< std::string > const names( option[ s ]() );
		for ( size_t i=1; i<= names.size(); ++i ) {
			filenames.push_back( names[i] );
		}
	}

	if ( filenames.empty() ) {
		basic::T("basic.options.util") << "Use either -s or -l to designate one or more start_files" << std::endl;
		utility_exit();
	}

	return filenames;
}


} // namespace options
} // namespace basic
