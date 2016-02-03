// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/StructFileRep.cc
/// @brief  The representation of a structure file
/// @author Andy Watkins
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- conversion to class from struct
/// @author Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_core_io_pdb_StructFileRep_hh
#define INCLUDED_core_io_pdb_StructFileRep_hh


// Unit headers
#include <core/io/pdb/Field.fwd.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.fwd.hh>
//#include <core/io/pdb/file_data_fixup.fwd.hh>

// Package headers
#include <core/io/StructFileRepOptions.fwd.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/io/Remarks.hh>
#include <core/io/CrystInfo.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/io/izstream.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

// C++ headers
#include <iostream>
#include <map>
#include <string>


namespace core {
namespace io {

typedef std::string String;

/// @brief   A class that contains information for individual atoms.
/// @details Only fields that are present in the PDB file will be initialized;
/// others will have the default value.
/// This class basically reflects the structure of 'ATOM' lines in PDB file format.
class AtomInformation
{
public:
	/// @brief default constructor to initialize all values
	AtomInformation() :
		isHet( false ),
		serial( 0 ),
		name( "" ),
		altLoc( ' ' ),
		resName( "" ),
		chainID( ' ' ),
		resSeq( 0 ),
		iCode( ' ' ),
		x( 0.0 ), y( 0.0 ), z( 0.0 ),
		occupancy( 0.0 ),
		temperature( 0.0 ),
		segmentID( "" ),
		element( "" ),
		formalcharge(0),
		terCount( 0 )
	{}

	/// For now, all member names have the same names as fields in PDB standard.
	bool isHet;
	int serial;
	String name;
	char altLoc;
	String resName;
	char chainID;
	int resSeq;
	char iCode;
	double x, y, z;
	double occupancy;
	double temperature;
	String segmentID;
	String element;
	signed short int formalcharge;
	int terCount; //< number of TER or END cards encountered prior to this

	/// @brief List of lower-numbered atoms that this atom is connected to.
	utility::vector1< core::Size > connected_indices;

	/// @brief Debug printing, serializing to Tracer like object.
	friend std::ostream& operator <<(std::ostream &os, AtomInformation const & ai) {
		os << "<AtomInformation>{" << "isHet=" << ai.isHet << " serial=" << ai.serial << " name=" << ai.name << " resName=" << ai.resName
			<< " chainID=" << ai.chainID << " resSeq=" << ai.resSeq
			<< " x=" << ai.x << " y=" << ai.y << " z=" << ai.z
			<< " temperature=" << ai.temperature
			<< " occupancy=" << ai.occupancy
			<< " segmentID=" << ai.segmentID
			<< " element=" << ai.element
			<< "}";
		return os;
	}
};  // class AtomInformation

typedef utility::vector0< AtomInformation > ChainAtoms;  // A representation of all atoms in a chain.

/// @brief   Intermediate format for easy construction of core::conformation::Residue objects.
/// @details Subset of data from "ATOM" lines that is shared by all atoms in a residue.
class ResidueInformation
{
public:
	/// @brief default constructor to initialize all values
	ResidueInformation();
	/// @brief Initialize the resName/chainID/etc. data from an atom,
	/// but do not add the atom to the atoms_ array.
	ResidueInformation(AtomInformation const & ai);

	bool operator==(ResidueInformation const & that) const;
	bool operator!=(ResidueInformation const & that) const;

	String const & resName() const;
	char chainID()  const;
	int  resSeq()   const;
	char iCode()    const;
	int  terCount() const;
	std::string const & segmentID() const;

	void resName(  String const & setting );
	void chainID(  char setting );
	void resSeq(   int setting );
	void iCode(    char setting );
	void terCount( int setting );
	void set_xyz( std::string const &atomname, Vector const &vect );
	void set_temp( std::string const &atomname, core::Real const &val );
	void segmentID( std::string const & setting );

	utility::vector1< AtomInformation > const & atoms() const;
	void append_atom( AtomInformation const & atoms );
	void rename_atom( std::string const & orig_name, std::string const & new_name );
	void append_atoms( ResidueInformation const & source );

	/// @brief A convience accessor for the
	std::map< std::string, Vector > const & xyz() const;
	std::map< std::string, double > const & temps() const; //< map of names to B-factors;  redundant but used a lot in reader


	/// @brief Returns a short, printable designation for this residue.
	String resid() const;

private:
	/// For now, all member names have the same names as fields in PDB standard.
	String resName_;
	char chainID_;
	int resSeq_;
	char iCode_;
	int terCount_; //< number of TER or END cards encountered prior to this
	utility::vector1< AtomInformation > atoms_;
	std::map< std::string, Vector > xyz_; //< map of names to coords;  redundant but used a lot in reader
	std::map< std::string, double > temps_; //< map of names to B-factors;  redundant but used a lot in reader

	std::string segmentID_;


};  // class ResidueInformation


/// @brief  A structure for storing information from .pdb MODRES records.
/// @author Labonte <JWLabonte@jhu.edu>
struct ModifiedResidueInformation {
	std::string resName;
	char chainID;
	int seqNum;
	char iCode;
	std::string stdRes;
	std::string comment;
};  // struct ModifiedResidueInformation


/// @brief  A structure for storing information from .pdb LINK records.
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

/// @brief  A structure for storing information from .pdb SSBOND records.
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

/// @brief  A structure for storing information from .pdb CISPEP records.
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


std::ostream&
operator <<( std::ostream &os, StructFileRep const & sfr );

/// @brief StructFileRep class. Hold data created from PDB file.
/// @details This class is intended to be a big, dumb data container.  As such, it has
/// no methods aside from its constructor (which initializes data that need to be
/// initialized), and all data are public.
class StructFileRep : public utility::pointer::ReferenceCount
{
public:

	// Intermediate representation of data for ease of creating Pose object.
	typedef std::map< std::string, double > ResidueTemps;
	//typedef std::map< std::string, ResidueTemps > Temps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef std::map< std::string, ResidueCoords > Coords;
	typedef utility::vector1< std::string > Strings;

public:

	/// @brief Constructor.
	StructFileRep();

	/// @brief empty destructor in C++ file to reduce number of necessary includes.
	~StructFileRep();

	/// @brief Clone operator: create a copy and return an owning pointer to the copy.
	StructFileRepOP clone() const;


public: //All data are public because this is a big, dumb data container.

	// An instance of StructFileRep should not preserve any 'state', so its data and methods are public.
	/// @brief The name of the input file, if any.
	std::string const & filename() const;
	std::string & filename();

	/// @brief ????
	std::string const & modeltag() const;
	std::string & modeltag();

	// PDB Title Section //////////////////////////////////////////////////////
	/// @brief    HEADER, TITLE, EXPDTA, KEYWDS, and COMPND records
	/// @details "header" is a misnomer, as it actually stores HEADER, TITLE, EXPDTA, KEYWDS, and COMPND records.
	HeaderInformationCOP  header() const;
	HeaderInformationOP & header();

	// Data for OBSLTE, SPLT, CAVEAT, NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, and/or JRNL records should be declared
	// here if ever implemented.
	/// @brief PDB remarks.
	RemarksCOP  remarks() const;
	RemarksOP & remarks();


	// PDB Primary Structure Section //////////////////////////////////////////
	// Data for DBREF, DBREF1, DBREF2, and/or SEQADV records should be declared here if ever implemented.
	std::map< char, utility::vector1< std::string > > const & chain_sequences() const;  // for storing SEQRES records
	std::map< char, utility::vector1< std::string > > & chain_sequences();  // for storing SEQRES records

	// map for storing MODRES records:
	// key is 6-character resID of the modified residue
	std::map< std::string, ModifiedResidueInformation > const & modres_map() const;
	std::map< std::string, ModifiedResidueInformation > & modres_map();


	// PDB Heterogen Section //////////////////////////////////////////////////
	// Data for HET should be declared here if ever implemented, though from what I can tell, HET gives no extra
	// information.
	/// @brief list for storing HETNAM records:
	/// @details first value of the pair is hetID\n
	/// second value of the pair is the chemical name field
	utility::vector1< std::pair< std::string, std::string > > const & heterogen_names() const;
	utility::vector1< std::pair< std::string, std::string > > & heterogen_names();

	// map for storing HETSYN records:
	//  key is hetID
	//  value is the chemical synonym list field
	std::map< std::string, utility::vector1< std::string > > const & heterogen_synonyms() const;
	std::map< std::string, utility::vector1< std::string > > & heterogen_synonyms();

	// map for storing FORMUL records:
	//  key is hetID
	//  value is the chemical formula, including a potential asterisk character
	std::map< std::string, std::string > const & heterogen_formulae() const;
	std::map< std::string, std::string > & heterogen_formulae();

	/// @brief map for storing ResidueType base (non-variant) names; parsed from HETNAM records:
	/// @details key is 6-character resID
	std::map< std::string, std::string > const & residue_type_base_names() const;
	std::map< std::string, std::string > & residue_type_base_names();


	// PDB Secondary Structure Section ////////////////////////////////////////
	// Data for HELIX and/or SHEET records should be declared here if ever implemented.


	// PDB Connectivity Annotation Section ////////////////////////////////////
	/// @brief   Map for storing SSBOND records.
	/// @details key is 6-character resID of 1st residue in ssbond
	/// @note    (A vector is needed because to futureproof if we ever handle weird disorder
	/// situations.)
	std::map< std::string, utility::vector1< SSBondInformation > > const & ssbond_map() const;
	std::map< std::string, utility::vector1< SSBondInformation > >       & ssbond_map();

	/// @brief   Map for storing LINK records.
	/// @details Key is 6-character resID of 1st residue in link.
	/// @note    (A vector is needed because a single saccharide residue can have multiple branches.)
	std::map< std::string, utility::vector1< LinkInformation > > const & link_map() const;
	std::map< std::string, utility::vector1< LinkInformation > >       & link_map();

	/// @brief   Map for storing CISPEP records.
	/// @details Key is 6-character resID of 1st residue in the peptide bond.
	/// @note    (A vector is needed because a single saccharide residue can have multiple branches.)
	std::map< std::string, CisPeptideInformation > const & cispep_map() const;
	std::map< std::string, CisPeptideInformation > & cispep_map();


	// PDB Miscellaneous Features Section /////////////////////////////////////
	// Data for SITE records should be declared here if ever implemented.


	// PDB Crystallographic and Coordinate Transformation Section /////////////
	/// @brief Crystallographic information.
	CrystInfo const & crystinfo() const;  //fpd:  CRYST1 line
	CrystInfo       & crystinfo();        //fpd:  CRYST1 line

	// Data for MTRIX, ORIGX, and/or SCALE records should be declared here if ever implemented.


	// PDB Coordinate Section /////////////////////////////////////////////////
	/// @brief The actual atomic coordinates are stored here.
	utility::vector0< ChainAtoms > const & chains() const;
	utility::vector0< ChainAtoms >       & chains();

	// PDB Connectivity Section ///////////////////////////////////////////////
	// Data for CONECT records should be declared here for consistency.

	// Non-PDB Stuff //////////////////////////////////////////////////////////
	/// @brief The foldtree, represented as a string.

	/// @details Each file outputter must figure out how to write this out in its output format.
	std::string const & foldtree_string() const;
	std::string       & foldtree_string();

	/// @brief The PDB comments, represented as a map of string->string.
	/// @details Each file outputter must figure out how to write this out in its output format.
	std::map< std::string, std::string > const & pdb_comments() const;
	std::map< std::string, std::string >       & pdb_comments();

	/// @brief A catch-all place to store additional data for output.
	/// @details Each file outputter must figure out how to write this out in its output format.
	std::string const & additional_string_output() const;
	std::string       & additional_string_output();

	/// @brief Append more string data to the additional_string_output_ string in the SFR.
	/// @details Retains all string data already added.
	void append_to_additional_string_output( std::string const &input_string );

private:
	std::string filename_;
	std::string modeltag_;
	HeaderInformationOP header_;
	RemarksOP remarks_;
	std::map< char, utility::vector1< std::string > > chain_sequences_;  // for storing SEQRES records
	std::map< std::string, ModifiedResidueInformation > modres_map_;
	utility::vector1< std::pair< std::string, std::string > > heterogen_names_;
	std::map< std::string, utility::vector1< std::string > > heterogen_synonyms_;
	std::map< std::string, std::string > heterogen_formulae_;
	std::map< std::string, std::string > residue_type_base_names_;
	std::map< std::string, utility::vector1< SSBondInformation > > ssbond_map_;
	std::map< std::string, utility::vector1< LinkInformation > > link_map_;
	std::map< std::string, CisPeptideInformation > cispep_map_;
	CrystInfo crystinfo_;  //fpd:  CRYST1 line
	utility::vector0< ChainAtoms > chains_;
	std::string foldtree_string_;
	std::map < std::string, std::string > pdb_comments_;
	std::string additional_string_output_;

};  // class StructFileRep

} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_StructFileRep_HH
