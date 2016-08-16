// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/StructFileRep.hh
/// @brief  Class/structure definitions for StructFileRep
/// @author Sergey Lyskov
/// @author Andy Watkins
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_StructFileRep_HH
#define INCLUDED_core_io_pdb_StructFileRep_HH

// Unit headers
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/HeaderInformation.fwd.hh>
#include <core/io/Remarks.fwd.hh>
#include <core/io/CrystInfo.hh>
#include <core/io/AtomInformation.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <map>
#include <string>


namespace core {
namespace io {

typedef utility::vector0< AtomInformation > ChainAtoms;  // A representation of all atoms in a chain.


// A structure for storing information from .pdb MODRES records.
/// @author Labonte <JWLabonte@jhu.edu>
struct ModifiedResidueInformation {
	std::string resName;
	char chainID;
	int seqNum;
	char iCode;
	std::string stdRes;
	std::string comment;
};  // struct ModifiedResidueInformation


// A structure for storing information from .pdb LINK records.
/// @author Labonte <JWLabonte@jhu.edu>
struct LinkInformation {
	std::string name1;  // 1st atom name
	std::string resName1;
	char chainID1;
	int resSeq1;
	char iCode1;
	std::string resID1;  // a 6-character resID, as defined elsewhere in StructFileRep (not from the PDB)
	std::string name2;  // 2nd atom name
	std::string resName2;
	char chainID2;
	int resSeq2;
	char iCode2;
	std::string resID2;
	core::Distance length;
};  // struct LinkInformation


// A structure for storing information from .pdb SSBOND records.
/// @author Watkins
struct SSBondInformation {
	std::string name1;  // 1st atom name - usually not read from PDB
	std::string resName1;
	char chainID1;
	int resSeq1;
	char iCode1;
	std::string resID1;  // a 6-character resID, as defined elsewhere in StructFileRep (not from the PDB)
	std::string name2;  // 2nd atom name - usually not read from PDB
	std::string resName2;
	char chainID2;
	int resSeq2;
	char iCode2;
	std::string resID2;
	core::Distance length;
	std::string sym1;
	std::string sym2;
};  // struct SSBondInformation


// A structure for storing information from .pdb CISPEP records.
/// @author Labonte <JWLabonte@jhu.edu>
struct CisPeptideInformation {
	std::string pep1;
	char chainID1;
	int seqNum1;
	char icode1;
	std::string pep2;
	char chainID2;
	int seqNum2;
	char icode2;
	//int modNum;
	core::Angle measure;
};  // struct CisPeptideInformation


/// @details This class is designed to store information about a structure in a
/// "non-Rosetta" way, such that it can serve as an intermediary between
/// various file formats and a Pose.
/// @remarks This class is intended to be light, consisting primarily of simple
/// accessors.
class StructFileRep : public utility::pointer::ReferenceCount {
public:  // Typedefs //////////////////////////////////////////////////////////
	typedef std::map< std::string, double > ResidueTemps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef std::map< std::string, ResidueCoords > Coords;
	typedef utility::vector1< std::string > Strings;


public:  // Standard Methods //////////////////////////////////////////////////
	/// @brief Constructor.
	StructFileRep();

	/// @brief empty destructor in C++ file to reduce number of necessary includes.
	~StructFileRep();

	/// @brief Create a copy and return an owning pointer to the copy.
	StructFileRepOP clone() const;


public:  // Accessors /////////////////////////////////////////////////////////
	/// @brief Access the name of the input/output file, if any.
	std::string const & filename() const { return filename_; }

	/// @brief Access the name of the input/output file, if any.
	std::string & filename() { return filename_; }


	/// @brief Access the model tag.
	std::string const & modeltag() const { return modeltag_; }

	/// @brief Access the model tag.
	std::string & modeltag() { return modeltag_; }


	// PDB Title Section //////////////////////////////////////////////////////
	/// @brief    Access HEADER, TITLE, EXPDTA, KEYWDS, and COMPND records.
	/// @details "header" is a misnomer, as it actually stores HEADER, TITLE,
	/// EXPDTA, KEYWDS, and COMPND records.
	HeaderInformationCOP header() const { return header_; }
	HeaderInformationOP & header() { return header_; }

	// Accessors for OBSLTE, SPLT, CAVEAT, NUMMDL, MDLTYP, AUTHOR, REVDAT,
	// SPRSDE, and/or JRNL records data should be declared here if ever
	// implemented.

	/// @brief Access PDB remarks.
	RemarksCOP remarks() const { return remarks_; }
	RemarksOP & remarks() { return remarks_; }


	// PDB Primary Structure Section //////////////////////////////////////////
	/// @brief  Access the sequences for each chain.
	std::map< char, utility::vector1< std::string > > const & chain_sequences() const { return chain_sequences_; }
	std::map< char, utility::vector1< std::string > > & chain_sequences() { return chain_sequences_; }

	/// @brief    Access map for storing MODRES records.
	/// @details  Key is 6-character resID of the modified residue.
	std::map< std::string, ModifiedResidueInformation > const & modres_map() const { return modres_map_; }
	std::map< std::string, ModifiedResidueInformation > & modres_map() { return modres_map_; }


	// PDB Heterogen Section //////////////////////////////////////////////////
	// Accessors for HET data should be defined here if ever implemented,
	// though from what I can tell, HET gives no extra information.  ~Labonte

	/// @brief   Access map for storing HETNAM records.
	/// @details key is hetID\n
	/// value is the chemical name field
	std::map< std::string, std::string > const & heterogen_names() const { return heterogen_names_; }
	std::map< std::string, std::string > & heterogen_names() { return heterogen_names_; }

	/// @brief   Access map for storing HETSYN records.
	/// @details key is hetID\n
	/// value is the chemical synonym list field
	std::map< std::string, utility::vector1< std::string > > const & heterogen_synonyms() const { return heterogen_synonyms_; }
	std::map< std::string, utility::vector1< std::string > > & heterogen_synonyms() { return heterogen_synonyms_; }

	/// @brief   Access map for storing FORMUL records.
	/// @details key is hetID\n
	/// value is the chemical formula, including a potential asterisk character
	std::map< std::string, std::string > const & heterogen_formulae() const { return heterogen_formulae_; }
	std::map< std::string, std::string > & heterogen_formulae() { return heterogen_formulae_; }

	/// @brief   Access map for storing ResidueType base (non-variant) names; parsed from HETNAM records:
	/// @details key is 6-character resID\n
	/// first value of pair is 3-letter-code; second value of pair is the base_name.
	std::map< std::string, std::pair< std::string, std::string > > const & residue_type_base_names() const { return residue_type_base_names_; }
	std::map< std::string, std::pair< std::string, std::string > > & residue_type_base_names() { return residue_type_base_names_; }


	// PDB Secondary Structure Section ////////////////////////////////////////
	// Accessors for HELIX and/or SHEET records data should be declared here
	// if ever implemented.




	// PDB Connectivity Annotation Section ////////////////////////////////////
	/// @brief   Access map for storing SSBOND records.
	/// @details key is 6-character resID of 1st residue in ssbond
	/// @note    (A vector is needed because to futureproof if we ever handle
	/// weird disorder situations.)
	std::map< std::string, utility::vector1< SSBondInformation > > const & ssbond_map() const { return ssbond_map_; }
	std::map< std::string, utility::vector1< SSBondInformation > > & ssbond_map() { return ssbond_map_; }


	/// @brief   Access map for storing LINK records.
	/// @details Key is 6-character resID of 1st residue in link.
	/// @note    (A vector is needed because a single saccharide residue can
	/// have multiple branches.)
	std::map< std::string, utility::vector1< LinkInformation > > const & link_map() const { return link_map_; }
	std::map< std::string, utility::vector1< LinkInformation > > & link_map() { return link_map_; }

	/// @brief   Access map for storing CISPEP records.
	/// @details Key is 6-character resID of 1st residue in the peptide bond.
	std::map< std::string, CisPeptideInformation > const & cispep_map() const { return cispep_map_; }
	std::map< std::string, CisPeptideInformation > & cispep_map() { return cispep_map_; }


	// PDB Miscellaneous Features Section /////////////////////////////////////
	// Accessors for SITE records data should be declared here if ever
	// implemented.


	// PDB Crystallographic and Coordinate Transformation Section /////////////
	/// @brief  Access crystallographic information.
	CrystInfo const & crystinfo() const { return crystinfo_; }
	CrystInfo & crystinfo() { return crystinfo_; }

	// Accessors for MTRIX, ORIGX, and/or SCALE records data should be declared
	// here if ever implemented.


	// PDB Coordinate Section /////////////////////////////////////////////////
	/// @brief  Access the actual atomic coordinates, stored as chains.
	utility::vector0< ChainAtoms > const & chains() const { return chains_; }
	utility::vector0< ChainAtoms > & chains() { return chains_; }


	// PDB Connectivity Section ///////////////////////////////////////////////
	// Accessors for CONECT records data should be declared here for
	// consistency.


	// Non-PDB Stuff //////////////////////////////////////////////////////////
	/// @brief   Access the The foldtree, represented as a string.
	/// @details Each file outputter must figure out how to write this out in
	/// its output format.
	std::string const & foldtree_string() const { return foldtree_string_; }
	std::string & foldtree_string() { return foldtree_string_; }

	/// @brief   Access The PDB comments, represented as a map of string->string.
	/// @details Each file outputter must figure out how to write this out in
	/// its output format.
	std::map< std::string, std::string > const & pdb_comments() const { return pdb_comments_; }
	std::map< std::string, std::string > & pdb_comments() { return pdb_comments_; }

	/// @brief   Access a catch-all place to store additional data for output.
	/// @details Each file outputter must figure out how to write this out in
	/// its output format.
	std::string const & additional_string_output() const { return additional_string_output_; }
	std::string & additional_string_output() { return additional_string_output_; }


	/// @brief   Append more string data to the additional_string_output_ string in the SFR.
	void append_to_additional_string_output( std::string const & input_string );


private:
	std::string filename_;
	std::string modeltag_;
	HeaderInformationOP header_;
	RemarksOP remarks_;
	std::map< char, utility::vector1< std::string > > chain_sequences_;  // for storing SEQRES records
	std::map< std::string, ModifiedResidueInformation > modres_map_;  // key is 6-character resID
	std::map< std::string, std::string > heterogen_names_;  // key is hetID
	std::map< std::string, utility::vector1< std::string > > heterogen_synonyms_;  // key is hetID
	std::map< std::string, std::string > heterogen_formulae_;  // key is hetID
	std::map< std::string, std::pair< std::string, std::string > > residue_type_base_names_;  // key is 6-char. resID
	std::map< std::string, utility::vector1< SSBondInformation > > ssbond_map_;  // key is 6-character resID
	std::map< std::string, utility::vector1< LinkInformation > > link_map_;  // key is 6-character resID
	std::map< std::string, CisPeptideInformation > cispep_map_;  // key is 6-character resID
	CrystInfo crystinfo_;  //fpd:  CRYST1 line
	utility::vector0< ChainAtoms > chains_;
	std::string foldtree_string_;
	std::map < std::string, std::string > pdb_comments_;
	std::string additional_string_output_;
};  // class StructFileRep

/// @brief  Debug printing, serializing to Tracer like object
std::ostream & operator<<( std::ostream & os, AtomInformation const & ai );

/// @brief  Output StructFileRep object to TR-like stream in human-readable format.
std::ostream & operator<<( std::ostream & os, StructFileRep const & sfr );

} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_StructFileRep_HH
