// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/setup/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/setup/util.hh>
#include <utility/io/izstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.setup.util" );

namespace protocols {
namespace stepwise {
namespace setup {

//////////////////////////////////////////////////////////
utility::vector1< std::string > load_s_and_l()
{
	using basic::options::option;
	using utility::vector1;
	using namespace basic::options::OptionKeys;

	// concatenate -s and -l flags together to get total list of PDB files
	vector1< std::string > pdb_file_names;
	if ( option[ in::file::s ].active() ) {
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
	}

	vector1< std::string > list_file_names;
	if ( option[ in::file::l ].active() ) {
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)
	}
	if ( option[ in::file::list ].active() ) {
		vector1< std::string > better_list_file_names;
		better_list_file_names= option[in::file::list ]().vector(); // make a copy (-list)
		for ( auto & better_list_file_name : better_list_file_names ) {
			list_file_names.push_back(better_list_file_name); // make a copy (-l)
		}
	}

	for ( auto filename : list_file_names ) {
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while ( getline(data, line) ) {
			pdb_file_names.push_back( std::string(line) );
		}
		data.close();
	}

	return pdb_file_names;
}


} //setup
} //stepwise
} //protocols
