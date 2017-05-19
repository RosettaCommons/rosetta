// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/PreloadedPDB.cc
/// @author SC Li & YK Ng (kalngyk@gmail.com)
/// @author Andy Watkins (amw579@stanford.edu)

#include <iostream>
#include <sstream>
#include <fstream>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <map>
#include <vector>

#include <protocols/cluster/calibur/SimPDB.fwd.hh>
#include <protocols/cluster/calibur/PreloadedPDB.hh>
#include <protocols/cluster/calibur/pdb_util.hh>

#include <utility/exit.hh> // for runtime_assert

namespace protocols {
namespace cluster {
namespace calibur {

/**
* Reasons for choosing 36,000 as the threshold for using disk access
* 1. 36,000 decoys of 100 residues each is about 44MB, which is still
* acceptable for a workstation PC. If user has less capabled hardware,
* the "-d" switch may be used.
* 2. For more than 36,000 decoys, the bottleneck is probably not in disk.
*/
unsigned int PreloadedPDB::ADVISED_THRESHOLD = 36000;

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

PreloadedPDB::PreloadedPDB():
	num_residue_( 0 ),
	mNumDecoy( 0 )
{}


/**
* Populate the PreloadedPDB with a silentfile
*/
void
PreloadedPDB::loadSilentFile( std::string const & filename ) {
	char buf[400];
	std::ifstream input( filename );
	std::string line;
	if ( !input ) {
		std::cerr << "Can't open silent file \"" << filename << "\"" << std::endl;
		exit(0);
	}

	silentfilename_ = filename;

	/**
	* Determine the number of residues from the silent file
	*/
	input.getline(buf, 400);
	line=buf;
	if ( line.substr(0, 9)!="SEQUENCE:" ) {
		std::cerr << "Silent file \"" << filename << "\" began with:" << std::endl
			<< "\t" << line << std::endl;
		exit(0);
	}

	input.getline(buf, 400);
	line=buf;
	if ( line.substr(0, 6)!="SCORE:" ) {
		std::cerr << "Silent file \"" << filename << "\" began with:" << std::endl
			<< "\t" << line << std::endl;
		exit(0);
	}

	input.getline(buf, 400);
	line=buf;
	if ( line.substr(0, 6)!="SCORE:" ) {
		std::cerr << "Silent file \"" << filename << "\" began with:" << std::endl
			<< "\t" << line << std::endl;
		exit(0);
	}

	int residueID = 0;
	for ( int i=1; !input.eof(); i++ ) {
		input.getline(buf, 400);
		line = buf;
		if ( line.substr(0, 6) == "SCORE:" ) break;
		residueID = toInt(line.substr(1, 4));
		if ( residueID != i ) {
			std::cerr << "Residue id " << residueID
				<< " out of sequence in silent file" << std::endl;
			exit(0);
		}
	}

	// AMW: does this assume residues are numbered sequentially from 1?
	runtime_assert( num_residue_ != 0 );
	num_residue_ = residueID;
	input.seekg(0);
	input.getline(buf, 400);
	input.getline(buf, 400);
	input.getline(buf, 400);

	/**
	* Read PDBs into filename2PDB
	*/
	SimPDBOP pdb( new SimPDB(num_residue_) );
	bool isNewPDB = true;
	int numResidue = 0;
	int decoyCount = 1;
	std::string key;
	do
	{
		input.getline(buf, 400);
		line = buf;

		if ( line.substr(0, 6) == "SCORE:" ) { // Old PDB done
			// Check if PDB has num_residue_
			if ( numResidue != num_residue_ ) {
				std::cerr << "Insufficient residues in the " << decoyCount
					<< "-th decoy in silent file" << std::endl;
				exit(0);
			}

			// Insert the pdb
			filename2PDB[key] = pdb;
			pdb->protein_file_name_ = key;
			center_residues(pdb->calpha_vector_, pdb->num_residue_);

			// Start a new pdb
			pdb = SimPDBOP( new SimPDB(num_residue_) );
			isNewPDB = true;
			numResidue = 0;

			decoyCount++;
		} else if ( line != "" ) { // Sometimes an empty string is read at eof
			if ( isNewPDB ) {
				// Get the filename
				std::string filename = line.substr(62, 400);
				std::string ext = filename.substr(filename.length()-4, 4);

				// If filename does not end in .pdb, generate a filename
				if ( strcmp(ext.c_str(), ".pdb") ) {
					std::stringstream ss;
					ss << "decoy" << decoyCount;
					key = ss.str();
				} else {
					key = filename;
				}

				isNewPDB = false;
			}

			residueID = toInt(line.substr( 0, 4));
			double x = toFloat(line.substr(35, 8));
			double y = toFloat(line.substr(44, 9));
			double z = toFloat(line.substr(54, 8));
			pdb->calpha_vector_[numResidue*3]   = x;
			pdb->calpha_vector_[numResidue*3+1] = y;
			pdb->calpha_vector_[numResidue*3+2] = z;

			if ( residueID != numResidue+1 ) {
				std::cout << residueID << "," << numResidue << std::endl;
				std::cerr << "Residue ID out of sequence in the " << decoyCount
					<< "-th decoy in silent file" << std::endl;
				exit(0);
			}

			numResidue++;
		}

	}
	while (!input.eof());
	input.close();

	// Insert the final pdb
	filename2PDB[key] = pdb;
	pdb->protein_file_name_ = key;

	mNumDecoy = decoyCount;

	/**
	* Generate the vector of filenames.
	* We can insert filename at the same time as inserting the SimPDB into
	* filename2PDB. However, doing this here has the advantage that if the
	* same filename is entered twice into filename2PDB, the name will not
	* be duplicated here.
	*/
	names_ = StringVecOP( new StringVec );
	for ( auto const & elem : filename2PDB ) {
		names_->push_back( elem.first );
	}
}


/**
* Populate the PreloadedPDB with the PDB files specified in a list.
*
* Reading of the PDB files is through SimPDB.
*/
void
PreloadedPDB::loadPDBFromList( std::string const & filename )
{
	std::ifstream input(filename);
	if ( !input ) {
		std::cerr << "Cannot find file \"" << filename << "\"" << std::endl;
		exit(0);
	}

	pdblistfilename_ = filename;

	/**
	* Read in names of PDB files
	*/
	//mNames = new std::vector<char *>(0);
	std::string line;
	names_ = StringVecOP( new StringVec );
	while ( std::getline( input, line ) ) {
		// Line has only one entry in this format: the PDB file name.
		names_->push_back(line);
	}
	input.close();

	mNumDecoy = names_->size();

	SimPDBOP pdb(new SimPDB);
	pdb->protein_file_name_ = (*names_)[0];
	pdb->num_residue_ = LONGEST_CHAIN;
	pdb->calpha_vector_.resize(3*LONGEST_CHAIN);
	pdb->read();

	num_residue_ = pdb->num_residue_;

	filename2PDB[(*names_)[0]] = pdb;

	for ( unsigned int i=1; i < names_->size(); i++ ) {
		SimPDBOP pdb(new SimPDB(num_residue_));
		pdb->protein_file_name_ = (*names_)[i];
		pdb->read();
		filename2PDB[(*names_)[i]] = pdb;
	}
}


SimPDBOP
PreloadedPDB::getSimPDB( std::string const & filename )
{
	return filename2PDB[filename];
}


}
}
}
