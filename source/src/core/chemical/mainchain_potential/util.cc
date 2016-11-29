// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/mainchain_potential/util.cc
/// @brief  Utility functions for mainchain torsional potentials.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <core/chemical/mainchain_potential/util.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.hh>

// Project Headers
//#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/string_util.hh>
//#include <utility/vector1.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//C++ includes

namespace core {
namespace chemical {
namespace mainchain_potential {

static basic::Tracer TR("core.chemical.mainchain_potential.util");

/// @brief Read a Shapovalov/Ramachandran-style mainchain torsion file, parse it for mainchain
/// potentials corresponding to a vector of ResidueType names, and return a map of (name->MainchainScoreTableCOP).
/// @param[in] filename The name of the file to read.
/// @param[in] res_type_names A vector of ResidueType names that this function will seek data for.
/// @param[in] use_polycubic_interpolation Should the MainchainScoreTables use polycubic interpolation?
/// @param[out] mainchain_score_table_map The output map of (name->MainchainScoreTableCOP).
void read_rama_map_file_shapovalov(
	std::string const &filename,
	utility::vector1< std::string > const &res_type_names,
	bool const use_polycubic_interpolation,
	std::map< std::string, MainchainScoreTableCOP > &mainchain_score_table_map
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	debug_assert( res_type_names.size() > 0 );

	// The input stream for the file:
	utility::io::izstream  iunit;
	// A string to hold the file contents:
	std::string filecontents;

	// Search in the three places:
	//   1) filename
	//   2) $DATABASE/[rama_pp_map]/filename
	//   3) $DATABASE/filename
	iunit.open( filename );
	if ( !iunit.good() ) {
		iunit.close();
		std::string full_filename = std::string(option[ corrections::score::rama_pp_map]())+"/"+filename;
		if ( !basic::database::open( iunit, full_filename ) ) {
			if ( !basic::database::open( iunit, filename ) ) {
				std::stringstream err_msg;
				err_msg << "Unable to open mainchain torsion score table '" << filename << "'.";
				utility_exit_with_message(err_msg.str());
			}
		}
	}
	utility::slurp( iunit, filecontents ); //Read the whole file into a string.
	iunit.close();

	//Now do the parsing:
	mainchain_score_table_map.clear();

	for ( core::Size i=1, imax=res_type_names.size(); i<=imax; ++i ) {
		runtime_assert_string_msg( !mainchain_score_table_map.count( res_type_names[i] ), "Error in core::chemical::mainchain_potential::read_rama_map_file_shapovalov(): The residue type " + res_type_names[i] + " has already been parsed.  (It was included multiple times in the list of types to parse.)" );
		MainchainScoreTableOP scoretable( new MainchainScoreTable );
		scoretable->parse_rama_map_file_shapovalov( filename, filecontents, res_type_names[i], use_polycubic_interpolation );
		mainchain_score_table_map[ res_type_names[i] ] = MainchainScoreTableCOP( scoretable );
	}
}

/// @brief Read a Shapovalov/Ramachandran-style mainchain torsion file, parse it for all the mainchain
/// potentials that it contains, and return a vector of pairs of (map name, MainchainScoreTableOP).
/// @param[in] filename The name of the file to read.
/// @param[in] use_polycubic_interpolation Should the MainchainScoreTables use polycubic interpolation?
/// @param[out] newtables The output vector of pairs of (name, MainchainScoreTableOP).
void read_rama_map_file_shapovalov(
	std::string const &filename,
	bool const use_polycubic_interpolation,
	utility::vector1< std::pair < std::string, MainchainScoreTableOP> > &newtables
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// The input stream for the file:
	utility::io::izstream  iunit;
	// A string to hold the file contents:
	std::string filecontents;

	// Search in the three places:
	//   1) filename
	//   2) $DATABASE/[rama_pp_map]/filename
	//   3) $DATABASE/filename
	iunit.open( filename );
	if ( !iunit.good() ) {
		iunit.close();
		std::string full_filename = std::string(option[ corrections::score::rama_pp_map]())+"/"+filename;

		//see if fullfilename exists (this is a waste but db::open throws an exception if file does not exist)
		utility::io::izstream temp( basic::database::full_name( full_filename, false ) );
		if (!temp.good()) {
			full_filename = filename;
		}

		if ( !basic::database::open( iunit, full_filename ) ) {
			std::stringstream err_msg;
			err_msg << "Unable to open mainchain torsion score table '" << filename << "'.";
			utility_exit_with_message(err_msg.str());
		}
	}
	utility::slurp( iunit, filecontents ); //Read the whole file into a string.
	iunit.close();

	//Go through the file and make a list of amino acids it contains.
	std::istringstream istream(filecontents);
	utility::vector1< std::string > namelist;
	do {
		std::string line;
		std::getline( istream, line );
		if ( istream.eof() ) break;
		utility::trim( line, " \n\t" ); //Strip terminal whitespace.
		line = line.substr( 0, line.find('#') ); //Strip anything following the number sign (comments).
		if ( line.empty() || line[0] == '\n' || line[0] == '\r' || line[0] == '@' ) continue;  //Skip blank lines and command lines.

		std::istringstream linestream(line); //Put the line into a stringstream for easy parsing.
		std::string curname;
		linestream >> curname;

		bool found(false);
		for ( core::Size i=1, imax=namelist.size(); i<=imax; ++i ) {
			if ( namelist[i].compare(curname) == 0 ) {
				found=true;
				break;
			}
		}
		if ( !found ) namelist.push_back(curname);
	} while(true);

	//Populate newtables
	newtables.clear();
	newtables.reserve( namelist.size() );
	for ( core::Size i=1, imax=namelist.size(); i<=imax; ++i ) {
		MainchainScoreTableOP scoretable( new MainchainScoreTable );
		scoretable->parse_rama_map_file_shapovalov( filename, filecontents, namelist[i], use_polycubic_interpolation );
		newtables.push_back( std::pair< std::string, MainchainScoreTableOP >(namelist[i], scoretable) );
	}
}

} //mainchain_potential
} //chemical
} //core
