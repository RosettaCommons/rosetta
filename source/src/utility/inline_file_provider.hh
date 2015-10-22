// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/inline_file_provider.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_utility_inline_file_provider_hh
#define INCLUDED_utility_inline_file_provider_hh

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <vector>


namespace utility {


// Forward
class Inline_File_Provider_Hook;

typedef utility::pointer::shared_ptr< Inline_File_Provider_Hook > Inline_File_Provider_HookOP;
typedef utility::pointer::shared_ptr< Inline_File_Provider_Hook const > Inline_File_Provider_HookCOP;

// base class for adding file-providing hooks into Inline_File_Provider
class Inline_File_Provider_Hook: public utility::pointer::ReferenceCount {
public:
	Inline_File_Provider_Hook(){}
	virtual bool request_file( const std::string & filename, std::string &result_data ) = 0;
	virtual void return_file_callback( const std::string &result_data, bool error ) = 0;
};

class Inline_File_Provider : public utility::SingletonBase< Inline_File_Provider > {
public:
	friend class utility::SingletonBase< Inline_File_Provider >;

private:
	Inline_File_Provider() {}
	static Inline_File_Provider * create_singleton_instance();

public:

	void init_static_inputs();
	void show_contents();

	void add_input_file( const std::string &filename, const std::string &contents );

	// files that will be rejected immediately as if they didn't exist on disk
	// this is used to reduce the number of necessary HTTP requests when we already know
	// certain files dont exist. Its essentailly a black list.
	void add_black_listed_file( const std::string &filename );

	void clear_input_files();

	void remove_input_file( const std::string &filename );

	bool file_exists( const std::string& filename );

	bool get_ostream( const std::string& filename, std::ostream **the_stream );

	bool get_istream( const std::string& filename, std::istream **the_stream );

	bool get_sstream( const std::string& filename, std::stringstream **the_stream );

	// Add a functor to the list of hooks. These will be called by file requests to allow
	// external addition of file sources.
	void add_file_provider_hook( const Inline_File_Provider_HookOP &new_hook );
private:

	bool is_black_listed_file( const std::string &filename );

	bool find_sstream( std::vector < std::pair < std::string, std::stringstream* > > &file_catalog, const std::string& filename, std::stringstream **the_stream );

	std::string standardise_filename( std::string filename );

	std::vector < std::pair < std::string, std::stringstream* > > input_files;

	std::vector < std::pair < std::string, std::stringstream* > > output_files;

	std::vector<Inline_File_Provider_HookOP> file_provider_hooks_;

	// files that will be rejected immediately as if they didn't exist on disk
	// this is used to reduce the number of necessary HTTP requests when we already know
	// certain files dont exist. Its essentailly a black list.
	std::vector < std::string > black_listed_files_;
};


}


#endif
