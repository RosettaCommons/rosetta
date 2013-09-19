// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data.hh
/// @brief  Declarations for FileData and related classes.
/// @author Sergey Lyskov

// Note: AVOID ACCESSING THE OPTIONS SYSTEM DIRECTLY IN THIS FILE, ESPECIALLY FOR PDB INPUT!
// Doing so will mean the Resource Manager may not work properly.
// Instead, modify FileDataOptions to include the option.


#ifndef INCLUDED_core_io_pdb_file_data_hh
#define INCLUDED_core_io_pdb_file_data_hh


// Unit headers
#include <core/io/pdb/Field.fwd.hh>
#include <core/io/pdb/HeaderInformation.hh>
#include <core/io/pdb/file_data.fwd.hh>

// Package headers
#include <core/io/pdb/file_data_options.fwd.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/CrystInfo.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>


// C++ headers
#include <iostream>
#include <map>
#include <string>


namespace core {
namespace io {
namespace pdb {

typedef std::string String;

/// @brief   A class that contains information for individual atoms.
/// @details Only fields that are present in the PDB file will be initialized;
/// others will have the default value.
/// This class basically reflects the structure of 'ATOM' lines in PDB file format.
class AtomInformation
{
public:
	///@brief default constructor to initialize all values
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
	element( "" ),
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
	String element;
	int terCount; //< number of TER or END cards encountered prior to this

	/// @brief Debug printing, serializing to Tracer like object.
	friend std::ostream& operator <<(std::ostream &os, AtomInformation const & ai) {
		os << "<AtomInformation>{" << "isHet=" << ai.isHet << " serial=" << ai.serial << " name=" << ai.name << " resName=" << ai.resName
		<< " chainID=" << ai.chainID << " resSeq=" << ai.resSeq
		<< " x=" << ai.x << " y=" << ai.y << " z=" << ai.z
		<< " temperature=" << ai.temperature
		<< " occupancy=" << ai.occupancy
		<< " element=" << ai.element
		<< "}";
		return os;
	}
};  // class AtomInformation

typedef std::vector<AtomInformation> AtomChain;


/// @brief   Intermediate format for easy construction of core::conformation::Residue objects.
/// @details Subset of data from "ATOM" lines that is shared by all atoms in a residue.
class ResidueInformation
{
public:
	///@brief default constructor to initialize all values
	ResidueInformation();

	ResidueInformation(AtomInformation const & ai);

	bool operator==(ResidueInformation const & that) const;
	bool operator!=(ResidueInformation const & that) const;

	/// For now, all member names have the same names as fields in PDB standard.
	String resid;  // 6-character (partial) identifier used by reader
	String resName;
	char chainID;
	int resSeq;
	char iCode;
	int terCount; //< number of TER or END cards encountered prior to this
	utility::vector1< AtomInformation > atoms;
	std::map< std::string, Vector > xyz; //< map of names to coords;  redundant but used a lot in reader
	std::map< std::string, double > temps; //< map of names to B-factors;  redundant but used a lot in reader
};  // class ResidueInformation


/// @brief  A structure for storing information from PDB LINK records.
/// @author Labonte
struct LinkInformation {
	std::string name1;  // 1st atom name
	std::string resName1;
	char chainID1;
	int resSeq1;
	char iCode1;
	std::string resID1;  // a 6-character resID, as defined elsewhere in FileData (not from the PDB)
	std::string name2;  // 2nd atom name
	std::string resName2;
	char chainID2;
	int resSeq2;
	char iCode2;
	std::string resID2;
	core::Distance length;
};  // struct LinkInformation


/// @brief FileData class. Hold data created from PDB file.
class FileData
{
public:
	// Intermediate representation of data for ease of creating Pose object.
	typedef std::map< std::string, double > ResidueTemps;
	typedef std::map< std::string, ResidueTemps > Temps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef std::map< std::string, ResidueCoords > Coords;
	typedef utility::vector1< std::string > Strings;

public:
	FileData();

	/// @brief empty destructor in C++ file to reduce number of necessary includes.
	~FileData();

public:  // An instance of FileData should not preserve any 'state', so its data and methods are public.
	std::string filename;
	std::string modeltag;

	// PDB Title Section //////////////////////////////////////////////////////
	// "header" is a misnomer, as it actually stores HEADER, TITLE, EXPDTA, KEYWDS, and COMPND records.
	HeaderInformationOP header;
	// Data for OBSLTE, SPLT, CAVEAT, NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, and/or JRNL records should be declared
	// here if ever implemented.
	pose::RemarksOP remarks;


	// PDB Primary Structure Section //////////////////////////////////////////
	// Data for DBREF, SEQADV, SEQRES, and/or MODRES records should be declared here if ever implemented.


	// PDB Heterogen Section //////////////////////////////////////////////////
	// Data for HET and/or FORMULA records should be declared here if ever implemented.

	// list for storing HETNAM records:
	//  first value of the pair is hetID
	//  second value of the pair is the chemical name field
	utility::vector1<std::pair<std::string, std::string> > heterogen_names;

	// map for storing carbohydrate ResidueType base (non-variant) names; parsed from HETNAM records:
	// key is 6-character resID
	std::map<std::string, std::string> carbohydrate_residue_type_base_names;

	// Data for HETSYN records should be declared here if ever implemented.


	// PDB Secondary Structure Section ////////////////////////////////////////
	// Data for HELIX and/or SHEET records should be declared here if ever implemented.


	// PDB Connectivity Annotation Section ////////////////////////////////////
	// Data for SSBOND records should be declared here if ever implemented.

	// map for storing LINK records:
	// key is 6-character resID of 1st residue in link
	// (A vector is needed because a single saccharide residue can have multiple branches.)
	std::map<std::string, utility::vector1<LinkInformation> > link_map;

	// Data for CISPEP records should be declared here if ever implemented.


	// PDB Miscellaneous Features Section /////////////////////////////////////
	// Data for SITE records should be declared here if ever implemented.


	// PDB Crystallographic and Coordinate Transformation Section /////////////
	pose::CrystInfo crystinfo;  //fpd:  CRYST1 line
	// Data for MTRIX, ORIGX, and/or SCALE records should be declared here if ever implemented.


	// PDB Coordinate Section /////////////////////////////////////////////////
	std::vector< AtomChain > chains;


	// PDB Connectivity Section ///////////////////////////////////////////////
	// Data for CONECT records should be declared here for consistency.

public:
	void initialize_header_information();

	HeaderInformationOP	header_information() const;

	void store_header_record(Record & R);

	/// @brief The header records can span multiple lines while the
	/// pdb_dynamic_parser is done line-wise. Finalizing the header
	/// information ensures that all the information is fully processed.
 	void finalize_header_information();

	/// @brief Make sure to call finalize_header_information before
	/// calling this.
	void fill_header_records(std::vector< Record > & VR) const;


	/// @brief Store (non-standard) polymer linkages in a map.
	void store_link_record(Record & record);


	/// @brief Store heterogen name information in a map.
	void store_heterogen_names(std::string const & hetID, std::string & text);

	/// @brief Parse heterogen name data for a given carbohydrate and save the particular base (non-variant)
	/// ResidueType needed in a map.
	void parse_heterogen_name_for_carbohydrate_residues(std::string const & text);

	/// @brief Return the PDB resName, chainID, resSeq, and iCode for the given Rosetta sequence position.
	ResidueInformation get_residue_information(core::pose::Pose const & pose, core::uint const seqpos,
			bool use_PDB=true, bool renumber_chains=false) const;

	///@brief Append pdb information to FileData for a single residue.
	void append_residue(
			core::conformation::Residue const & rsd,
			core::Size & atom_index,
			core::pose::Pose const & pose,
			bool preserve_crystinfo = false
	);

	/// @brief Fill FileData object using information from given Pose object.
	void init_from_pose(core::pose::Pose const & pose);

	/// @brief Fill FileData object  using information from given Pose object and a set of options.
	void init_from_pose(core::pose::Pose const & pose, FileDataOptions const & options);

	/// @brief Fill FileData object using information from given Pose object,
	/// for a specified subset of residues
	void init_from_pose( core::pose::Pose const & pose, utility::vector1< core::Size > const & residue_indices );

	/// @brief Writes  <pose>  to a given stream in PDB file format
	static void dump_pdb(
		core::pose::Pose const & pose,
		std::ostream & out,
		std::string const & tag="",
		bool write_fold_tree = false
	);

	/// @brief Writes  <pose>  to a PDB file, returns false if an error occurs
	static bool dump_pdb(
		core::pose::Pose const & pose,
		std::string const & file_name,
		std::string const & tag="",
		bool write_fold_tree = false
	);

	static void dump_pdb(
		core::pose::Pose const & pose,
		std::ostream & out,
		utility::vector1< core::Size > const & residue_indices,
		std::string const & tag=""
	);


	bool update_atom_information_based_on_occupancy( AtomInformation & ai, FileDataOptions const & options,
			std::string const & resid ) const;

	///@brief randomize missing density
	void randomize_missing_coords( AtomInformation & ai ) const;

	/// @brief Debug printing
	friend std::ostream& operator <<(std::ostream &os, FileData const &);


	/// @brief Convert FileData into set of residues, sequences, coordinates etc.
	void create_working_data(
		utility::vector1< ResidueInformation > & rinfo
	);

	/// @brief Convert FileData into set of residues, sequences, coordinates etc. using a set of options
	void create_working_data(
		utility::vector1< ResidueInformation > & rinfo,
		FileDataOptions const & options
	);
};  // class FileData

/// @brief Adds data to the end of a pdb that are not a standard part of the pdb format.
void write_additional_pdb_data(
	std::ostream & out,
	pose::Pose const & pose,
	io::pdb::FileData const & fd,
	bool write_fold_tree = false
);

/// @brief Builds a pose into  <pose>, without repacking or optimizing
/// hydrogens; using the full-atom ResidueTypeSet
void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename
);

/// @brief Builds a pose into  <pose>, without repacking or optimizing
/// hydrogens; using the full-atom ResidueTypeSet and a set of options.
void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename,
	PDB_DReaderOptions const & pdr_options
);

void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
);

void build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	PDB_DReaderOptions const & pdr_options
);


void build_pose_as_is1(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	FileDataOptions const & options
);

bool is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & pdb_name,
	core::chemical::ResidueTypeCOPs const & rsd_type_list,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps);

bool is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & pdb_name,
	core::chemical::ResidueTypeCOPs const & rsd_type_list,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps,
	FileDataOptions const & options);


void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices
);

void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices,
	FileDataOptions const & options
);

void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices
);

void pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices,
	FileDataOptions const & options
);

} // namespace pdb
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_file_data_HH
