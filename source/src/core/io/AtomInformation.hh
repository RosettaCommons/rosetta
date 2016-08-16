// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/ResidueInformation.hh
/// @brief  Structure definition for ResidueInformation
/// @author Sergey Lyskov


#ifndef INCLUDED_core_io_AtomInformation_HH
#define INCLUDED_core_io_AtomInformation_HH

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace core {
namespace io {

// A structure that contains information for individual atoms.
/// @details Only fields that are present in the PDB file will be initialized;
/// others will have the default value.\n
/// This class basically reflects the structure of 'ATOM' lines in PDB file
/// format.
struct AtomInformation {
	/// @brief default constructor to initialize all values (except the
	/// connected indices)
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
		formalcharge( 0 ),
		terCount( 0 )
	{}

	// For now, all member names have the same names as fields in PDB standard.
	bool isHet;
	int serial;
	std::string name;
	char altLoc;
	std::string resName;
	char chainID;
	int resSeq;
	char iCode;
	double x, y, z;
	double occupancy;
	double temperature;
	std::string segmentID;
	std::string element;
	signed short int formalcharge;
	int terCount; //< number of TER or END cards encountered prior to this

	/// @brief List of lower-numbered atoms that this atom is connected to.
	utility::vector1< core::Size > connected_indices;
};  // struct AtomInformation

} // namespace io
} // namespace core

#endif  // INCLUDED_core_io_AtomInformation_HH
