// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/pdb/PDBDataType.hh
/// @brief   Enumeration definition for PDBDataType.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_PDBDataType_HH
#define INCLUDED_core_io_pdb_PDBDataType_HH

namespace core {
namespace io {
namespace pdb {

/// @brief  Enumerators for all of the data types used in PDB files, according to:
/// http://wwpdb.org/documentation/file-format-content/format33/sect1.html#Field
enum PDBDataType {
	AChar = 1,  // alphabetic char (A-Z or a-z).
	Atom,  // string of size 4
	Character,
	Continuation,  // int
	Date,  // 9-char string (DD-MMM-YY)
	IDcode, // string of size 4
	Integer,
	//Token,  // string followed by colon and a space (variable length; used in Specification Lists)
	List,  // vector1< string > (separated by commas)
	LString,  // string (preserves all spacing, left-justified)
	Real_type,  // Real
	Record_name,  // string of size 6 (left-justfied)
	Residue_name,  // string of size 3
	SList,  // vector< string > (separated by semi-colons)
	Specification,  // string with Token followed by a value
	Specification_List,  // sequence of Specifications separated by semi-colons
	String_type,
	SymOP,  // int (4 to 6 digits, right-justified, nnnMMM
	N_PDB_DATA_TYPES = SymOP
};

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_PDBDataType_HH
