// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chromophore/ReadResidueCoordinatesFromPDB.cc
/// @brief re-reading residue's coordinates from a PDB file
/// @author Nina Bozhanova (nbozhanova@gmail.com)

// Project headers:
#include <protocols/chromophore/ReadResidueCoordinatesFromPDB.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>
#include <utility/io/izstream.hh>

#include <core/io/ResidueInformation.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/pdb/pdb_reader.hh>

#include <core/io/StructFileRep.hh> // AUTO IWYU For StructFileRep


static basic::Tracer TR( "protocols.chromophore.ReadResidueCoordinatesFromPDB" );


namespace protocols {
namespace chromophore {

/// @brief Default constructor.
ReadResidueCoordinatesFromPDB::ReadResidueCoordinatesFromPDB() = default;

/// @brief Destructor.
ReadResidueCoordinatesFromPDB::~ReadResidueCoordinatesFromPDB() = default;

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
ReadResidueCoordinatesFromPDBOP
ReadResidueCoordinatesFromPDB::clone() const {
	return utility::pointer::make_shared< ReadResidueCoordinatesFromPDB >( *this );
}


void
ReadResidueCoordinatesFromPDB::read_coordinates_from_file (
	std::string const & filename,
	utility::vector1 <std::tuple <int, char> > const & residues_to_read) {

	utility::io::izstream infile(filename);
	if ( !infile ) {
		utility_exit_with_message("Cannot open file " + filename);
	} else {
		read_coordinates (infile, residues_to_read);
	}
}


void
ReadResidueCoordinatesFromPDB::read_coordinates (
	std::istream & instream,
	utility::vector1 <std::tuple <int, char> > const & residues_to_read) {

	core::io::StructFileRep sfr;
	parse_pdb (instream, sfr);

	// Reorganize information about the residues for which we need coordinates
	std::map <char, utility::vector1 <int> > residues_to_read_map;
	for ( core::Size i(1); i <= residues_to_read.size(); i++ ) {
		int resnum = std::get<0>(residues_to_read[i]);
		char chain_name = std::get<1>(residues_to_read[i]);
		residues_to_read_map[chain_name].push_back(resnum);
	}

	utility::vector1<core::io::ResidueInformation> rinfo;

	// Iterate all chains in the PDB
	for ( core::Size ch=0; ch < sfr.chains().size(); ++ch ) {
		// Get information about the first atom in the chain
		core::io::AtomInformation dummy_ai( sfr.chains()[ch][0] );
		// Is it the chain that we are interested in?
		if ( residues_to_read_map.find(dummy_ai.chainID) != residues_to_read_map.end() ) {
			utility::vector1 <int> resnums = residues_to_read_map[dummy_ai.chainID];
			// Iterate all atoms in the chain
			for ( core::Size i=0; i < sfr.chains()[ch].size(); ++i ) {
				// Get atom information
				core::io::AtomInformation ai( sfr.chains()[ch][i] );
				// Is this atom from the residue we need?
				if ( std::find(resnums.begin(), resnums.end(), ai.resSeq ) != resnums.end() ) {
					// Get residue information
					core::io::ResidueInformation residue = core::io::ResidueInformation (ai);
					// Have we already added this residue to the rinfo?
					if ( rinfo.size() == 0 || rinfo.back() != residue ) {
						rinfo.push_back(residue);
						rinfo.back().append_atom(ai);
					} else {
						rinfo.back().append_atom(ai);
					}
				}
			}
		}
	}

	// Throw an error if no residues were found
	if ( rinfo.size() == 0 ) {
		utility_exit_with_message( "No residues fulfill the specified criteria.");
	} else if ( rinfo.size() != residues_to_read.size() ) {
		utility_exit_with_message( "The number of found residues does not match the number of requested residues.");
	} else {
		for ( core::Size i(1); i <= rinfo.size(); i++ ) {
			save_residue_coordinates (rinfo[i], rinfo[i].resSeq(), rinfo[i].chainID());
		}
	}
}


void
ReadResidueCoordinatesFromPDB::parse_pdb (
	std::istream & instream,
	core::io::StructFileRep & sfr) {

	std::string file_content;
	utility::slurp( instream, file_content );
	sfr = core::io::pdb::create_sfr_from_pdb_file_contents( file_content);
}


void
ReadResidueCoordinatesFromPDB::save_residue_coordinates (
	core::io::ResidueInformation const & residue,
	int const resnum,
	char const chain_name) {

	// I guess we don't want to accidentally add coordinates of the residue twice
	if ( !coordinates_exist(resnum, chain_name) ) {
		std::tuple < int, char > const key (resnum, chain_name);
		for ( core::Size i(1); i <= residue.atoms().size(); ++i ) {
			core::io::AtomInformation const & atom( residue.atoms()[i] );
			std::string const name( utility::stripped_whitespace( atom.name ) );
			core::Vector const xyz_vector ( atom.x, atom.y, atom.z );
			coordinates_[key].push_back(std::tuple<std::string, core::Vector>(name, xyz_vector));
		}
	} else {
		utility_exit_with_message( "The coordinates for this residue seems to be already saved");
	}
}


utility::vector1 <std::tuple <std::string, core::Vector> >
ReadResidueCoordinatesFromPDB::get_residue_coordinates(
	int const resnum,
	char const chain_name) const {

	if ( coordinates_exist(resnum, chain_name) ) {
		std::tuple < int, char > const key (resnum, chain_name);
		return coordinates_.at(key);
	} else {
		// Throw an error if empty
		utility_exit_with_message( "No coordinates are available for the residue.");
	}
}


bool
ReadResidueCoordinatesFromPDB::coordinates_exist(
	int const resnum,
	char const chain_name) const {

	std::tuple < int, char > key (resnum, chain_name);
	return (coordinates_.find(key) != coordinates_.end() );
}


core::Size
ReadResidueCoordinatesFromPDB::number_of_residues () const {

	return coordinates_.size();
}


} //chromophore
} //protocols
