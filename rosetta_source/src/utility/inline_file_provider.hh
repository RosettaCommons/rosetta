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
#include <vector>


namespace utility {


class Inline_File_Provider {

	private:
		Inline_File_Provider(){
		}

	public:
		static Inline_File_Provider* get_instance();
		void show_contents();
		bool file_exists( const std::string& filename );
		bool get_istream( const std::string& filename, std::istream **the_stream );
	private:
		std::string standardise_filename( std::string filename );
		static Inline_File_Provider* instance_;
		std::vector < std::stringstream* > streambucket;
};


}


#endif

