// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/SheetDB.hh
/// @brief ALlows rosetta to read/manipulate a database of sheets
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_SheetDB_hh
#define INCLUDED_protocols_denovo_design_components_SheetDB_hh

// Unit headers

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>

// Package headers
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// C++ headers

// Utility Headers

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief this is the actual map used to store the sheets
typedef utility::vector1< core::pose::PoseCOP > SheetList;
typedef std::map< std::string, SheetList > MapByLengths;
typedef std::map< std::string, MapByLengths > MapByOrientations;
typedef std::map< core::Size, MapByOrientations > MapByStrands;

using namespace protocols::denovo_design::architects;

/// @brief class used to create/manipulate a database of sheets using a list of pdb files
class SheetDB : public utility::pointer::ReferenceCount {
public:
	/// @brief default constructor
	SheetDB();

	/// @brief query the database for a sheet and retrieve the list -- if not found, tries to read from disk
	SheetList
	sheet_list( Lengths const & lengths, StrandOrientations const & orientations, RegisterShifts const & shifts );

	/// @brief query the database for a sheet and retrieve the list -- if not found, tries to read from disk
	SheetList
	sheet_list_const( Lengths const & lengths, StrandOrientations const & orientations, RegisterShifts const & shifts ) const;

	/// @brief query the database, returns whether data exists
	bool
	has_data( Lengths const & lengths, StrandOrientations const & orientations, RegisterShifts const & shifts ) const;

	/// @brief set location of sheet database files
	inline void set_db_path( std::string const & db_path ) { db_path_ = db_path; }

	/// @brief add a sheet to the database
	void set_sheetlist(
		SheetList const & list,
		core::Size const nstrands,
		std::string const & in_orientations,
		std::string const & in_lengths );

	/// @brief adds a sheet list as it is presented
	void set_sheetlist_as_is(
		SheetList const & list,
		core::Size const nstrands,
		std::string const & orientations,
		std::string const & lengths );

	/// @brief add a sheet to the database
	void add_sheet(
		core::pose::PoseOP pose,
		core::Size const nstrands,
		std::string const & orientations,
		std::string const & lengths_shifts,
		bool const check_descriptor_against_pose );

	/// @brief adds sheet(s) based on a pose -- returns # of sheets added
	int add_sheets_from_pose( core::pose::Pose & pose );

	/// @brief adds sheet(s) based on a strand pairing set and pose -- returns # of sheets added
	int add_sheets_from_spairs( core::pose::Pose const & pose,
		protocols::fldsgn::topology::StrandPairingSet & spairs,
		protocols::fldsgn::topology::SS_Info2 const & ss_info );

	/// @brief queries the database for all strand lengths
	utility::vector1< core::Size > nstrand_values() const;

	/// @brief queries the database for all orientations with nstrands strands
	utility::vector1< std::string >
	orientations( core::Size const nstrands ) const;

	/// @brief queries the database for all lengths/shifts with nstrands strands and given orientations
	utility::vector1< std::string >
	lengths_shifts( core::Size const nstrands, std::string const & orientations ) const;

	static core::Size const MIN_STRAND_LEN = 4;
	static core::Size const SHEETINFO_REMARK = 994;
	static core::Size const SHEETSCORE_REMARK = 995;

	// mutators
public:
	inline void set_idealize( bool const idealize_val ) { idealize_ = idealize_val; }

protected:
	/// @brief reads sheetlist from the file database, and saves it to this database
	void load_sheetlist_from_file( core::Size const nstrands, std::string const & orientations, std::string const & strandinfo );

private:
	/// @brief query the database for a sheet and retrieve the list -- if not found, tries to read from disk
	SheetList
	sheet_list( core::Size const nstrands, std::string const & orientations, std::string const & strandinfo );

	/// @brief queries the database to see if the requested data exists
	bool
	has_data( core::Size const nstrands, std::string const & orient, std::string const & strandinfo ) const;

	/// @brief query the database for a sheet and retrieve the list
	SheetList
	sheet_list_const( core::Size const nstrands, std::string const & orientations, std::string const & strandinfo ) const;

	/// @brief query the database for a sheet and retrieve the list with all lengths and shifts included
	SheetList
	sheet_list( core::Size const nstrands, std::string const & orientations ) const;


private:
	/// @brief idealize poses before adding to db? (Default = true);
	bool idealize_;
	/// @brief the location of the database files
	std::string db_path_;
	/// @brief the guts of the database
	MapByStrands sheet_db_;
};

/// @brief counts number of strands based on number of jumps. Pose must be a disconnected sheet
core::Size num_strands( core::pose::Pose const & pose );

/// @brief gets a pair of strings, which refer to orientations, and lengths/shifts. Pose secstruct MUST be set
std::pair< std::string, std::string >
find_orientations_and_lengths( core::pose::Pose const & pose );

/// @brief creates a string key for lengths/shifts based on the given vectors for strand lengths and register shifts
std::string
make_strand_info_str( Lengths const & lengths, RegisterShifts const & shifts );

core::pose::PoseOP
reverse_chains( core::pose::Pose const & pose );

int
add_to_pose(
	core::pose::PoseOP newpose,
	core::pose::Pose const & pose,
	core::Size s_start,
	core::Size s_stop );

utility::vector1< core::pose::PoseOP >
extract_sheets_from_strandlist(
	utility::vector1< core::pose::PoseOP > const & strands,
	core::Size const nstrands );

utility::vector1< core::pose::PoseOP >
extract_sheets_from_pose(
	core::pose::Pose const & pose,
	protocols::fldsgn::topology::StrandPairings const & spairs,
	protocols::fldsgn::topology::SS_Info2 const & ss_info,
	bool const idealize );

std::string reverse_orientations( std::string const & orient );

utility::vector1< bool >
parse_orientations( std::string const & orientations );

std::pair< utility::vector1< core::Size >, utility::vector1< int > >
parse_lengths( std::string const & lengths );

std::string reverse_lengths( std::string const & orientations, std::string const & lengths );

std::string const &
choose_canonical( std::string const & l1, std::string const & l2 );

std::pair< std::string, std::string >
canonicalize( std::string const & orientations, std::string const & lengths );

} // components
} // denovo_design
} // protocols

#endif
