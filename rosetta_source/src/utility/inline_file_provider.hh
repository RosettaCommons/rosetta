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
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>
#include <sstream>
#include <vector>


namespace utility {


class Inline_File_Provider {

	private:
		Inline_File_Provider(){
		}

	public:
		static Inline_File_Provider* get_instance();
		void init_static_inputs();
    void show_contents();

    void add_input_file( const std::string &filename, const std::string &contents ); 
	  void clear_input_files();
    void remove_input_file( const std::string &filename );
	
    bool file_exists( const std::string& filename );
		
		bool get_ostream( const std::string& filename, std::ostream **the_stream );
		
		bool get_istream( const std::string& filename, std::istream **the_stream );
		bool get_sstream( const std::string& filename, std::stringstream **the_stream );
	private:
    
    bool find_sstream( std::vector < std::pair < std::string, std::stringstream* > > &file_catalog, const std::string& filename, std::stringstream **the_stream );

		std::string standardise_filename( std::string filename );
		static Inline_File_Provider* instance_;
		
    std::vector < std::pair < std::string, std::stringstream* > > input_files;

		std::vector < std::pair < std::string, std::stringstream* > > output_files;
};


}


#endif

