// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chromophore/ChromophoreDataReader.cc
/// @brief implementation of reading information about correspondence between FP chromophore's
/// and standard Rosetta residue's atom names from a provided file
/// @author Nina Bozhanova (nbozhanova@gmail.com)

// Project headers:
#include <protocols/chromophore/ChromophoreDataReader.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <fstream>

static basic::Tracer TR( "protocols.chromophore.ChromophoreDataReader" );


namespace protocols {
namespace chromophore {

/// @brief Default constructor.
ChromophoreDataReader::ChromophoreDataReader():
	was_initialized_(false)
{}

/// @brief Destructor.
ChromophoreDataReader::~ChromophoreDataReader() = default;

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
ChromophoreDataReaderOP
ChromophoreDataReader::clone() const {
	return utility::pointer::make_shared< ChromophoreDataReader >( *this );
}

std::map <std::string, std::string>
ChromophoreDataReader::get_residue_names_map (core::Size const & residue_number) const {
	if ( residue_name_maps_.find(residue_number) != residue_name_maps_.end() ) {
		return residue_name_maps_.at(residue_number);
	} else {
		utility_exit_with_message( "The requested residue number is not found." );
	}
}


core::Size
ChromophoreDataReader::number_of_residues () const {
	return residue_name_maps_.size();
}


bool
ChromophoreDataReader::map_for_residue_exists (core::Size const & residue_number) const {
	if ( number_of_residues() != 0 ) {
		return (residue_name_maps_.find(residue_number) != residue_name_maps_.end() );
	} else {
		return false;
	}
}


void
ChromophoreDataReader::initialize_from_file (std::string const & filename) {

	if ( was_initialized_ ) {
		utility_exit_with_message("One file was already loaded by this object");
	}

	utility::io::izstream infile(filename);
	if ( !infile ) {
		utility_exit_with_message("Cannot open file " + filename);
	} else {
		create_residue_names_maps (infile);
	}

	was_initialized_ = true;
}


void
ChromophoreDataReader::create_residue_names_maps (std::istream & infile) {

	if ( was_initialized_ ) {
		utility_exit_with_message("One file was already loaded by this object");
	}

	// vector of tuples containing information about input (CRO) atom name, Rosetta atom name replacement, and the amino acid number in the CRO
	utility::vector1 <std::tuple <std::string, std::string, core::Size> > meaningful_file_content;
	parse_file( infile, meaningful_file_content );
	// Does the file look ok?
	if ( meaningful_file_content.size() == 0 ) {
		utility_exit_with_message( "No information in the input file. Please check your input file." );
	} else {
		for ( core::Size i(1); i <= meaningful_file_content.size(); i++ ) {
			// Add information about input atom name : Rosetta atom name correspondence to the appropriate map (based on the residue number)
			core::Size const residue_number = std::get<2>(meaningful_file_content[i]);
			std::string const & input_atom_name = std::get<0>(meaningful_file_content[i]);
			std::string const & rosetta_atom_name = std::get<1>(meaningful_file_content[i]);
			residue_name_maps_[residue_number][input_atom_name] = rosetta_atom_name;
		}
	}

	was_initialized_ = true;

}


void
ChromophoreDataReader::parse_file (std::istream & infile,
	utility::vector1 <std::tuple <std::string, std::string, core::Size> > & meaningful_file_content) {
	// Expected format:
	// # Input (CRO) atom name | Rosetta atom name replacement | aa number in the CRO
	// # residue #1 (Thr in the immature chromophore)
	// N1   N    1
	// CA1  CA   1
	std::string file_line;
	while ( std::getline(infile, file_line) ) {
		// Skip empty lines
		if ( file_line.empty() ) {
			continue;
		} else if ( file_line[0] == '#' ) {
			// Skip the line if it is a comment
			continue;
		} else {
			// Process the line
			std::tuple <std::string, std::string, core::Size> parsed_line;
			parse_line (file_line, parsed_line);
			meaningful_file_content.push_back( parsed_line );
		}
	}
}


void
ChromophoreDataReader::parse_line (std::string & file_line,
	std::tuple <std::string, std::string, core::Size> & parsed_line ) {

	// Split the line
	utility::vector1 <std::string> dataline_content = utility::split_whitespace(file_line);
	// Check the number of white space-delimited columns in the line
	// If not exactly 3, throw out an error
	if ( dataline_content.size() != 3 ) {
		utility_exit_with_message( "Wrong number of columns in the parsed line. Please check your input file." );
	} else {
		// Check whether the atom names look reasonable (no more than 4 characters)
		if ( dataline_content[1].size() >= 5 || dataline_content[2].size() >= 5 ) {
			utility_exit_with_message("Atom names are expected to be no longer than 4 characters: " + file_line);
		} else if ( !utility::is_string_numeric(dataline_content[3]) ) {
			// Check whether the last column is a number
			utility_exit_with_message( "The last column is expected to be a number but it is not. Please check your input file." );
		} else {
			core::Size const residue_number = std::stoi(dataline_content[3]);
			parsed_line = std::make_tuple(dataline_content[1], dataline_content[2], residue_number);
		}
	}
}


} // chromophore
} // protocols
