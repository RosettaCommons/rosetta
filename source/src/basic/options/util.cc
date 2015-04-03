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

	// -screening_list files
	if( option[ screening_list ].user() ) {
		utility::vector1<std::string> const list_files( option[ screening_list ]() );
		if(list_files.size() != 2)
		{
			utility_exit_with_message("-in:file:screening_list currently takes exactly 2 file lists as arguments");
		}

		utility::io::izstream outer_file(list_files[1].c_str());
		utility::vector1<std::string> outer_paths;
		if(!outer_file.good())
		{
			utility_exit_with_message("unable to open list file: "+list_files[1]);
		}
		while (outer_file.good())
		{
			std::string name;
			outer_file.getline(name);
			if(outer_file.good()) outer_paths.push_back(name);
		}
		outer_file.close();

		utility::io::izstream inner_file(list_files[2].c_str());
		utility::vector1<std::string> inner_paths;
		if(!inner_file.good())
		{
			utility_exit_with_message("unable to open list file: "+list_files[2]);
		}
		while(inner_file.good())
		{
			std::string name;
			inner_file.getline(name);
			if(inner_file.good()) inner_paths.push_back(name);
		}
		inner_file.close();

		for(size_t i = 1; i <= outer_paths.size();++i)
		{
			std::string outer_path(outer_paths[i]);
			for(size_t j = 1; j <= inner_paths.size();++j)
			{
				std::string inner_path(inner_paths[j]);
				filenames.push_back(outer_path+" "+inner_path);
			}
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
