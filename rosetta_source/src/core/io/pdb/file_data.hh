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
///
/// @brief
/// @author Sergey Lyskov

#ifndef INCLUDED_core_io_pdb_file_data_hh
#define INCLUDED_core_io_pdb_file_data_hh


// Unit headers
#include <core/io/pdb/file_data.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/pose/Remarks.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>


#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/types.hh>

// C++ headers
#include <iostream>
#include <map>
#include <string>


namespace core {
namespace io {
namespace pdb {

typedef std::string String;

/// @brief Atom_information contain information for individual atom.
/// only fields that predent in file will be initialazed, other will have default value.
/// That class basicaly reflect structure of 'ATOM' line in PDB file format.
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
};

typedef std::vector<AtomInformation> AtomChain;


/// @brief Intermediate format for easy construction of core::conformation::Residue objects.
/// Subset of data from "ATOM" lines that is shared by all atoms in a residue.
class ResidueInformation
{
public:
	///@brief default constructor to initialize all values
	ResidueInformation() :
		resid( "" ),
		resName( "" ),
		chainID( ' ' ),
		resSeq( 0 ),
		iCode( ' ' ),
		terCount( 0 ),
		atoms(),
		xyz(),
		temps()
	{}

	ResidueInformation(AtomInformation const & ai) :
		resid( "" ),
		resName( ai.resName ),
		chainID( ai.chainID ),
		resSeq( ai.resSeq ),
		iCode( ai.iCode ),
		terCount( ai.terCount ),
		atoms(),
		xyz(),
		temps()
	{}

	bool operator==(ResidueInformation const & that) const
	{
		return resName == that.resName && chainID == that.chainID && resSeq == that.resSeq && iCode == that.iCode && terCount == that.terCount;
	}
	// Stupid language that allows separate definitions of == and !=
	bool operator!=(ResidueInformation const & that) const
	{ return !(*this == that); }

	/// For now, all member names have the same names as fields in PDB standard.
	String resid; //< 6-character (partial) identifier used by reader
	String resName;
	char chainID;
	int resSeq;
	char iCode;
	int terCount; //< number of TER or END cards encountered prior to this
	utility::vector1< AtomInformation > atoms;
	std::map< std::string, Vector > xyz; //< map of names to coords;  redundant but used a lot in reader
	std::map< std::string, double > temps; //< map of names to B-factors;  redundant but used a lot in reader
};

///
/// @brief FileData class. Hold data created from PDB file.
///
class FileData
{
public:
	FileData() : remarks(new pose::Remarks) {}

	/// @brief empty destructor in C++ file to reduce number of necessary includes.
	~FileData();

	 /// only one data member, that should not preserve any 'state' - so it is public.
	std::vector< AtomChain > chains;
	//std::vector< RemarkInfo > remarks;
	pose::RemarksOP remarks;
	std::string filename;
	std::string modeltag;

	/// @brief Fill FileData structure using information from given pose object.
	void init_from_pose(core::pose::Pose const & pose);

	/// @brief Fill FileData structure using information from given pose object,
	/// for a specified subset of residues
	void
	init_from_pose( core::pose::Pose const & pose, utility::vector1< core::Size > const & residue_indices );

	///@brief randomize missing density
	void randomize_missing_coords( AtomInformation & ai );


	/// @brief Debug printing
	friend std::ostream& operator <<(std::ostream &os, FileData const &);

	/// @brief Writes  <pose>  to a given stream in PDB file format
	static
	void dump_pdb(
		core::pose::Pose const & pose,
		std::ostream & out,
		std::string const & tag="",
		bool write_fold_tree = false
	);

	/// @brief Writes  <pose>  to a PDB file, returns false if an error occurs
	static
	bool dump_pdb(
		core::pose::Pose const & pose,
		std::string const & file_name,
		std::string const & tag="",
		bool write_fold_tree = false
	);

	static
	void
	dump_pdb(
		core::pose::Pose const & pose,
		std::ostream & out,
		utility::vector1< core::Size > const & residue_indices,
		std::string const & tag=""
	);

	// Intermediate representation of date for easy of creating Pose object.
	typedef std::map< std::string, double > ResidueTemps;
	typedef std::map< std::string, ResidueTemps > Temps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef std::map< std::string, ResidueCoords > Coords;
	typedef utility::vector1< std::string > Strings;

	/// @brief Convert FileData in to set of residues, sequences, coordinats etc.
	void create_working_data(
		utility::vector1< ResidueInformation > & rinfo
	);

	///@brief appends pdb information for a single residue
	void append_residue(
		core::conformation::Residue const & rsd,
		core::Size & atom_index,
		core::pose::Pose const & pose
	);

};

void
write_additional_pdb_data(
	std::ostream & out,
	pose::Pose const & pose,
	io::pdb::FileData const & fd,
	bool write_fold_tree = false
);

/// @brief Builds a pose into  <pose>, without repacking or optimizing
/// hydrogens; using the fullatom ResidueTypeSet
void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename
);

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
);

void
build_pose_as_is1(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing
);

bool
is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & pdb_name,
	core::chemical::ResidueTypeCAPs const & rsd_type_list,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps);


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices
);

void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices
);


} // namespace pdb
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pdb_file_data_HH
