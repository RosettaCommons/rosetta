// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/SimPDB.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)
/// @author Andy Watkins (amw579@stanford.edu)

#include <iostream>
#include <fstream>

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iomanip>

#include <map>
#include <vector>

#include <protocols/cluster/calibur/SimPDB.hh>
#include <protocols/cluster/calibur/PreloadedPDB.hh>
#include <protocols/cluster/calibur/pdb_util.hh>

#include <memory>

namespace protocols {
namespace cluster {
namespace calibur {

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Memory and default values for the static variables

// Modifies how SimPDB reads PDB files
int SimPDB::s_residue = 1;
int SimPDB::e_residue = LONGEST_CHAIN;
char * SimPDB::chains = strdup("AC ");

// Set these two fields to tell SimPDB to use the PreloadedPDB mechanism
PreloadedPDB * SimPDB::preloadedPDB = nullptr;
bool SimPDB::preloadPDB = true;


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// General purpose constructors which handle the preloadedPDB mechanism

SimPDB::SimPDB( std::string const & aFileName )
{
	if ( preloadPDB ) {
		/*for(int i=0; i<len; i++)
		{
		calpha_vector[i]=new double [3];
		}*/
		SimPDBOP pdb = preloadedPDB->filename2PDB[aFileName];
		protein_file_name_ = pdb->protein_file_name_; // copies
		num_residue_ = pdb->num_residue_;
		//calpha_vector.resize(3*num_residue_);
		// copy by copy ctor
		calpha_vector_ = pdb->calpha_vector_;
		//memcpy(calpha_vector, pdb->calpha_vector_, 3 * num_residue_ * sizeof(double));
	} else {
		protein_file_name_ = aFileName;
		num_residue_ = LONGEST_CHAIN;
		calpha_vector_.resize(3*num_residue_);
		read();
	}
}

SimPDB::SimPDB( std::string const & aFileName, int len )
{
	if ( preloadPDB ) {
		/*for(int i=0; i<len; i++)
		{
		calpha_vector[i]=new double [3];
		}*/
		SimPDBOP pdb = preloadedPDB->filename2PDB[aFileName];
		protein_file_name_ = pdb->protein_file_name_;
		num_residue_ = pdb->num_residue_;
		calpha_vector_ = pdb->calpha_vector_;
	} else {
		protein_file_name_ = aFileName;
		num_residue_ = len;
		calpha_vector_.resize(3*num_residue_);
		read();
	}
}


SimPDB::~SimPDB() = default;


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Constructors which DO NOT handle the preloadedPDB mechanism
//
// They should be called from PreloadedPDB only, since it will need to
// bypass the mechanism

SimPDB::SimPDB() = default;

SimPDB::SimPDB(int len)
{
	num_residue_ = len;
	calpha_vector_.resize(3 * len );
}



//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
* Reads a PDB file from disk. Do not read from PDB files anywhere else.
*/
void
SimPDB::read()
{
	std::ifstream input(protein_file_name_);
	if ( !input ) {
		std::cerr << "Cannot find protein file " << protein_file_name_.c_str() << std::endl;
		exit(0);
	}
	std::cout.flush();
	char buf[400];
	num_residue_ = 0;
	double x, y, z;
	//char c = 'a';
	bool read = false;

	int prevID = -10000;
	int count = 0;

	int CA_number = 1;

	//double cx = 0;
	//double cy = 0;
	//double cz = 0;
	int count3 = 0;
	//mSquaredSum=0;
	while ( !input.eof() ) {
		input.getline(buf, 400);
		std::string line=buf;
		if ( line.substr(0, 3) == "TER" && read == true ) break;
		if ( line.substr(0, 6) == "ENDMDL" ) break;

		if ( line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM" ) {
			continue;
		}

		if ( line.substr(13, 4) == "CA  " || line.substr(13, 4) == " CA "
				|| line.substr(13, 4) == "  CA" || line.substr(13, 2) == "CA" ) {
			// At this point a CA atom has been discovered
			// We want to further filter it based on the following two
			// criteria: chain and region.

			// Check if the chain which this atom belongs to is to be included
			bool include_chain = false;
			for ( char * c = SimPDB::chains; *c; c++ ) {
				if ( toupper(line[21]) == toupper(*c) ) {
					include_chain = true;
					break;
				}
			}
			if ( !include_chain ) {
				CA_number++;
				continue;
			}

			// Check if the CA atom is within the region to analyze
			if ( CA_number < SimPDB::s_residue ) {
				CA_number++;
				continue;
			} else if ( CA_number > SimPDB::e_residue ) break;

			read = true;
			int residueID = toInt(line.substr(22, 6));
			if ( residueID == prevID ) continue;

			prevID = residueID;
			x = toFloat(line.substr(30, 8));
			y = toFloat(line.substr(38, 8));
			z = toFloat(line.substr(46, 8));
			//std::string AAType = line.substr();
			count3 = 3*count;
			calpha_vector_[count3]   = x;
			calpha_vector_[count3+1] = y;
			calpha_vector_[count3+2] = z;
			//mSquaredSum+=x*x+y*y+z*z;
			count++;

			CA_number++;
		}
	}//while

	num_residue_ = count;
	input.close();

	center_residues(calpha_vector_, num_residue_);
}


}
}
}
