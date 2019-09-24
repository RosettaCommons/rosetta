// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief declaration of implementation class for abstract class Residue
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_AtomICoor_fwd_hh
#define INCLUDED_core_chemical_AtomICoor_fwd_hh


// Unit headers

// Project headers
#include <core/types.hh>

// Utility headers

// C++ headers
#include <string>

namespace core {
namespace chemical {

class AtomICoor;

/// @brief ICoordAtomIDType
///
/// - INTERNAL: atoms which inherently belong to this ResidueType
/// - POLYMER_LOWER: atom at the polymer lower connection, such as backbone C in
/// the previous residue (N-term side)
/// - POLYMER_UPPER: atom at the polymer upper connection, such as backbone N in
/// the next residue (C-term side)
/// - CONNECT: atoms from a non-adjacent residue which connect to this residue
/// by non-polymer connection, such as disulfide
///
/// @details If you add anything to this enum, be sure to update the string_to_icoord_type() function.
enum class ICoordAtomIDType {
	INTERNAL = 1,
	POLYMER_LOWER,
	POLYMER_UPPER,
	CONNECT
};

/// @brief Convert a string designation into the corresponding ICoordAtomIDType enum.
ICoordAtomIDType
string_to_icoord_type( std::string const & );

/// @brief Get the connection number from a string representing an CONNECT type (e.g. 4 from `CONN4`)
Size
get_connection_number( std::string const & );

}
}

#endif
