// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/io/pdb/alt_codes_io.hh
/// @brief   Database input/output function declarations for alternative PDB 3-letter-code data.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_io_alt_codes_io_HH
#define INCLUDED_core_io_alt_codes_io_HH

// C++ headers
#include <map>
#include <string>

//Utility header
#include <utility/io/util.hh>


namespace core {
namespace io {

typedef std::map< std::string, std::tuple< std::string, std::string, utility::vector1< std::string> > > AltCodeMap;

/// @brief  Return a mapping of alternative PDB 3-letter codes to a paired set including the Rosetta 3-letter code and,
/// optionally, any HETNAM information for the residue desired.
AltCodeMap read_alternative_3_letter_codes_from_database_file( std::string const & filename );

}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_pdb_alt_codes_io_HH
