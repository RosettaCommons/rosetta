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
#include <utility/exit.hh>

// Singleton instance and mutex static data members
namespace utility {

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< Inline_File_Provider >::singleton_mutex_{};
template <> std::atomic< Inline_File_Provider * > utility::SingletonBase< Inline_File_Provider >::instance_( 0 );
#else
template <> Inline_File_Provider * utility::SingletonBase< Inline_File_Provider >::instance_( 0 );
#endif

}

namespace utility {

Inline_File_Provider * Inline_File_Provider::create_singleton_instance()
{
	return new Inline_File_Provider;
}


void Inline_File_Provider::show_contents(){

}

void Inline_File_Provider::init_static_inputs(){

}

void Inline_File_Provider::add_input_file( const std::string &filename, const std::string &contents ){
	std::stringstream *newstream = new std::stringstream(contents);
	std::string filtered_filename( standardise_filename( filename ) );
	std::cout << "Creating inline input file: '" << filename << "' = (" << filtered_filename << ")" << std::endl;
	input_files.push_back( std::make_pair( filtered_filename, newstream ) );
}

void Inline_File_Provider::add_black_listed_file( const std::string &filename ){
	black_listed_files_.push_back( filename );
}

class predicate_cmp_filename
{
public:
	predicate_cmp_filename( std::string filename ): filename_(filename) {}

	bool operator() ( std::pair < std::string, std::stringstream* > &a ) const
	{
		return a.first == filename_;
	}
private:
	std::string filename_;
};


void Inline_File_Provider::clear_input_files(){
	input_files.clear();
}

void Inline_File_Provider::remove_input_file( const std::string &filename ){

	std::vector < std::pair < std::string, std::stringstream* > >::iterator found;

	predicate_cmp_filename pred(filename);

	do{
		found = std::find_if( input_files.begin(), input_files.end(), pred );
		if ( found != input_files.end() ) {
			input_files.erase( found );
		} else {
			break;
		}
	}while(true);

}

std::string Inline_File_Provider::standardise_filename( std::string filename ){
	if ( filename.substr(0,2) == "./" ) {
		filename = filename.substr(2);
	}

	std::string no_dup_slash;
	for ( unsigned int i = 0; i < filename.size(); ++i ) {
		if ( (no_dup_slash.size() > 0) && // previous chars existed
				(filename[i] == '/') &&  // current char is a forward slash
				(no_dup_slash[no_dup_slash.size()-1] == '/') // and last added character is a forward slash
				) {
			continue;
		}
		no_dup_slash = no_dup_slash + filename[i];
	}
	return no_dup_slash;
}

bool Inline_File_Provider::file_exists( const std::string& filename )
{

	// std::string filtered_filename( standardise_filename( filename ) );
	// std::cout << "Looking for inline file: '" << filtered_filename << "'" << std::endl;
	// // first find the data stupid simple search
	//  // look through input_files
	// for( std::vector < std::pair < std::string, std::stringstream* > >::iterator it = input_files.begin();
	//      it != input_files.end(); ++it ){
	//  if( it->first == filtered_filename ){
	//   return true;
	//  }
	// }
	//
	// // now look through output_files
	// for( std::vector < std::pair < std::string, std::stringstream* > >::iterator it = output_files.begin();
	//      it != output_files.end(); ++it ){
	//  if( it->first == filtered_filename ){
	//   return true;
	//  }
	// }

	// determine if the file exists by opening it. On non-physical local filesystems there sometimes is not specific "does file exist"
	// mechanism. Trying to access the resource is the only way.
	std::istream *temp;
	return get_istream( filename, &temp );
}

bool Inline_File_Provider::get_ostream( const std::string& filename, std::ostream **the_stream )
{
	std::stringstream *newstream = new std::stringstream( );

	std::string filtered_filename( standardise_filename( filename ) );
	std::cout << "Creating inline output file: '" << filename << "' = (" << filtered_filename << ")" << std::endl;
	output_files.push_back( std::make_pair( filtered_filename, newstream ) );
	(*the_stream) = newstream;
	return true;
}


bool Inline_File_Provider::is_black_listed_file( const std::string &filename ){
	std::vector < std::string >::iterator found;
	found = std::find( black_listed_files_.begin(), black_listed_files_.end(), filename);
	return ( found != black_listed_files_.end() );
}


bool Inline_File_Provider::get_istream( const std::string& filename, std::istream **the_stream ){

	std::string standard_filename =  standardise_filename( filename );
	if ( is_black_listed_file(standard_filename) ) return false;


	if ( standard_filename.length() == 0 ) {
		return false;
	}

	// check last character
	if ( *(standard_filename.rbegin()) == '/' ) {
		// cannot download a directory - only files. Reject this
		return false;
	}

	std::stringstream *the_sstream(0);
	bool result = get_sstream( filename, &the_sstream );
	if ( result ) {
		(*the_stream) = the_sstream;
		return true;
	}

	// if here, we haven't found the file. Now go through all the hooks (if any) and ask each if they have the file. If so, take the file, add it to our file store and then
	// return the result.

	for ( std::vector<Inline_File_Provider_HookOP>::const_iterator it = file_provider_hooks_.begin();
			it != file_provider_hooks_.end(); ++it ) {
		std::string file_data;

		// try and request the file from that resource - these calls can take a long time to return(!) because the resource can be slow.
		// this call is blocking !
		if ( (*it)->request_file( standard_filename, file_data ) ) {
			// add the contents under the filename
			add_input_file( filename, file_data );
			// grab the stream
			if ( get_sstream( filename, &the_sstream ) ) {
				(*the_stream) = the_sstream;
				return true;
			}
		}
	}

	// if here we have failed to find the resource
	return false;
}

bool Inline_File_Provider::find_sstream( std::vector < std::pair < std::string, std::stringstream* > > &file_catalog, const std::string& filename, std::stringstream **the_stream ){
	std::vector < std::pair < std::string, std::stringstream* > >::iterator it = file_catalog.begin();
	for ( ; it != file_catalog.end(); ++it ) {
		if ( it->first == filename ) {
			break;
		}
	}
	if ( it != file_catalog.end() ) {
		(*the_stream) = it->second;
		return true;
	}
	return false;
}

bool Inline_File_Provider::get_sstream( const std::string& filename, std::stringstream **the_stream ){
	std::string filtered_filename( standardise_filename( filename ) );
	std::cout << "Looking for inline file: '" << filtered_filename << "'" << std::endl;

	if ( find_sstream( input_files, filtered_filename, the_stream) ) return true;
	if ( find_sstream( output_files, filtered_filename, the_stream) ) return true;


	return false;
}

void Inline_File_Provider::add_file_provider_hook( const Inline_File_Provider_HookOP &new_hook ){
	file_provider_hooks_.push_back( new_hook );
}

}
