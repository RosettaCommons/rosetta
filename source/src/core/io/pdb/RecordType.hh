// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/io/pdb/RecordType.hh
/// @brief   Enumeration definition for RecordType.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_pdb_RecordType_HH
#define INCLUDED_core_io_pdb_RecordType_HH

namespace core {
namespace io {
namespace pdb {

/// @brief  Enumerators for all of the record types used in PDB files, according to:
/// http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
enum RecordType {
	// Title Section //////////////////////////////////////////////////////////
	HEADER = 1,
	OBSLTE,
	TITLE,
	SPLIT,
	CAVEAT,
	COMPND,
	SOURCE,
	KEYWDS,
	EXPDTA,
	NUMMD,
	MDLTYP,
	AUTHOR,
	REVDAT,
	SPRSDE,
	JRNL,
	REMARK,

	// Primary Structure Section //////////////////////////////////////////////
	DBREF,
	DBREF1,
	DBREF2,
	SEQADV,
	SEQRES,
	MODRES,

	// Heterogen Section //////////////////////////////////////////////////////
	HET,
	HETNAM,
	HETSYN,
	FORMUL,

	// Secondary Structure Section ////////////////////////////////////////////
	HELIX,
	SHEET,

	// Connectivity Annotation Section ////////////////////////////////////////
	SSBOND,
	LINK,
	CISPEP,

	// Miscellaneous Features Section /////////////////////////////////////////
	SITE,

	// Crystallographic & Coordinate Transformation Section ///////////////////
	CRYST1,
	ORIGX1,
	ORIGX2,
	ORIGX3,
	SCALE1,
	SCALE2,
	SCALE3,
	MTRIX1,
	MTRIX2,
	MTRIX3,

	// Coordinate Section /////////////////////////////////////////////////////
	MODEL,
	ATOM,
	ANISOU,
	TER,
	HETATM,
	ENDMDL,

	// Connectivity Section ///////////////////////////////////////////////////
	CONECT,

	// Bookkeeping Section ////////////////////////////////////////////////////
	MASTER,
	END,

	// Rosetta Type ///////////////////////////////////////////////////////////
	UNKNOW,

	N_RECORD_TYPES = UNKNOW
};

}  // namespace pdb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_RecordType_HH
