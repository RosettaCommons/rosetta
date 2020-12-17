// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chromophore/ReadResidueCoordinatesFromPDB.hh
/// @brief re-reading residue's coordinates from a PDB file
/// @author Nina Bozhanova (nbozhanova@gmail.com)


#ifndef INCLUDED_protocols_chromophore_ReadResidueCoordinatesFromPDB_hh
#define INCLUDED_protocols_chromophore_ReadResidueCoordinatesFromPDB_hh

#include <protocols/chromophore/ReadResidueCoordinatesFromPDB.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <core/io/ResidueInformation.hh>
#include <core/io/StructFileRep.hh>
#include <numeric/xyzVector.hh>

namespace protocols {
namespace chromophore {

/// @brief re-reading residue's coordinates from a PDB file
/// @author Nina Bozhanova (nbozhanova@gmail.com)
class ReadResidueCoordinatesFromPDB : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	ReadResidueCoordinatesFromPDB();

	/// @brief Destructor.
	~ReadResidueCoordinatesFromPDB() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	ReadResidueCoordinatesFromPDBOP clone() const;

	// Get coordinates for the residue specified by the residue number and the chain
	utility::vector1 <std::tuple <std::string, core::Vector> > get_residue_coordinates(int const resnum, char const chain_name) const;

	// Open the PDB file and read the residue(s) coordinates (specified by the residue number and the chain)
	void read_coordinates_from_file (std::string const & filename, utility::vector1 <std::tuple <int, char> > const & residues_to_read);

	// Read the residue(s) coordinates (specified by the residue number and the chain)
	void read_coordinates (std::istream & instream, utility::vector1 <std::tuple <int, char> > const & residues_to_read);

	// Do we have any coordinates stored for the residue?
	bool coordinates_exist(int const resnum, char const chain_name) const;

	// How many residues we have information about?
	core::Size number_of_residues () const;

private:

	// Get sfr from the PDB file
	void parse_pdb (std::istream & instream, core::io::StructFileRep & sfr);

	// Get residue coordinates
	void save_residue_coordinates (core::io::ResidueInformation const & residue, int const resnum, char const chain_name);

private:

	// A map that contains all stored coordinates
	// The key is a tuple of the residue number and the residue chain
	// The value is a vector of tuples consisting of all residue's atom names and the corresponding xyz coordinates
	std::map < std::tuple < int, char >, utility::vector1 <std::tuple <std::string, core::Vector> > > coordinates_;

};

} //chromophore
} //protocols

#endif //INCLUDED_protocols_chromophore_ReadResidueCoordinatesFromPDB_hh
