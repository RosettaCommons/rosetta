// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chromophore/ChromophoreDataReader.hh
/// @brief header of class for reading information about correspondence between FP chromophore's
/// and standard Rosetta residue's atom names
/// @note This class is not currently threadsafe
/// @author Nina Bozhanova (nbozhanova@gmail.com)


#ifndef INCLUDED_protocols_chromophore_ChromophoreDataReader_hh
#define INCLUDED_protocols_chromophore_ChromophoreDataReader_hh

#include <protocols/chromophore/ChromophoreDataReader.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <map>
#include <string>

namespace protocols {
namespace chromophore {

class ChromophoreDataReader : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	ChromophoreDataReader();

	/// @brief Destructor.
	~ChromophoreDataReader() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	ChromophoreDataReaderOP clone() const;

	// Get information for the particular residue
	std::map <std::string, std::string> get_residue_names_map (core::Size const & residue_number) const;

	// Open the file, read and store information form it
	void initialize_from_file (std::string const & filename);

	// The actual reading of the input file and storing information
	void create_residue_names_maps (std::istream & instream);

	// Do we have any information on the residue?
	bool map_for_residue_exists (core::Size const & residue_number) const;

	// How many residues we have information about?
	core::Size number_of_residues () const;

private:

	// Choose which lines to parse and parse them
	void parse_file (std::istream & infile, utility::vector1 <std::tuple <std::string, std::string, core::Size> > & meaningful_file_content);

	// Parsing of an individual line
	void parse_line (std::string & file_line, std::tuple <std::string, std::string, core::Size> & parsed_line);

private:

	// The outer map is a map over residue numbers, the inner map contains atom name mappings for that residue
	std::map < core::Size, std::map <std::string, std::string> > residue_name_maps_;

	// Has chromophore data been already loaded?
	bool was_initialized_;

};

} // chromophore
} // protocols

#endif //INCLUDED_protocols_chromophore_ChromophoreDataReader_hh
