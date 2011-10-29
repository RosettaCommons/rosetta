// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file inline_file_provider.cc 


#include <utility/inline_file_provider.hh>

// C/C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <utility/static_database.hh>

namespace utility {


Inline_File_Provider* Inline_File_Provider::get_instance(){
	if(!instance_) instance_ = new Inline_File_Provider();
	return instance_;
}

void Inline_File_Provider::show_contents(){
//	for( unsigned int i = 0; i < static_database_size; ++i ){
//		std::cout << "FILE:     " << static_database[i][0] << std::endl;
//		std::cout << "CONTENTS: " << static_database[i][1] << std::endl;
//		std::cout << "-----------------------" << std::endl;
//	}
}

std::string Inline_File_Provider::standardise_filename( std::string filename ){
	if( filename.substr(0,2) == "./" ){
		filename = filename.substr(2);
	}

	std::string no_dup_slash;
	for( unsigned int i = 0; i < filename.size(); ++i){
		if( (no_dup_slash.size() > 0) && // previous chars existed 
		    (filename[i] == '/') &&  // current char is a forward slash
				(no_dup_slash[no_dup_slash.size()-1] == '/') // and last added character is a forward slash 
			){
			continue;
		}
		no_dup_slash = no_dup_slash + filename[i];
	}
	return no_dup_slash;
}

bool Inline_File_Provider::file_exists( const std::string& filename )
{
	std::string filtered_filename( standardise_filename( filename ) );

	std::cout << "Looking for inline file: '" << filtered_filename << "'" << std::endl; 
	// first find the data stupid simple search
	unsigned int i = 0;
	for( i = 0; i < static_database_size; ++i ){
		if( std::string( static_database[i][0] ) == filtered_filename ){
			return true;	
		}
	}
	
	// now look through outfiles
	for( std::vector < std::pair < std::string, std::stringstream* > >::iterator it = output_files.begin();
	     it != output_files.end(); ++it ){
		if( it->first == filtered_filename ){
			return true;
		}
	}
	return false;
}

bool Inline_File_Provider::get_ostream( const std::string& filename, std::ostream **the_stream )
{
	std::cout << "Creating inline file: '" << filename << "'" << std::endl; 
	std::stringstream *newstream = new std::stringstream( );

	std::string filtered_filename( standardise_filename( filename ) );
	output_files.push_back( std::make_pair( filtered_filename, newstream ) ); 
	(*the_stream) = newstream;
	return true;
}




bool Inline_File_Provider::get_istream( const std::string& filename, std::istream **the_stream ){
	std::stringstream *the_sstream;
	bool result = get_sstream( filename, &the_sstream );
	(*the_stream) = the_sstream;
	return result;
}

bool Inline_File_Provider::get_sstream( const std::string& filename, std::stringstream **the_stream ){
	std::string filtered_filename( standardise_filename( filename ) );
	std::cout << "Looking for inline file: '" << filtered_filename << "'" << std::endl; 
	
	// first find the data stupid simple search
	unsigned int i = 0;
	for( i = 0; i < static_database_size; ++i ){
		if( std::string( static_database[i][0] ) == filtered_filename ){
			break;
		}
	}
	// did we find it ?
	if( i <  static_database_size ){
		int data = i;
		const char *the_data = static_database[i][1];
		std::stringstream *newstream = new std::stringstream( the_data );
		streambucket.push_back( newstream );
		(*the_stream) = streambucket[ streambucket.size()-1 ];
		return true;
	}

	// now look through outfiles
	std::vector < std::pair < std::string, std::stringstream* > >::iterator it = output_files.begin();
	for( ; it != output_files.end(); ++it ){
		if( it->first == filtered_filename ){
			break;	
		}
	}

	if( it != output_files.end() ){
		(*the_stream) = it->second; 
		return true;
	}
	
	return false;
}


Inline_File_Provider* Inline_File_Provider::instance_ = NULL;

}



