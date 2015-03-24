// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/ring_conformer_io.hh
/// @brief   Database input/output function declarations for ring-conformer-specific data.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_core_chemical_ring_conformer_io_HH
#define INCLUDED_core_chemical_ring_conformer_io_HH

// Unit header
#include <core/chemical/RingConformerSet.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <string>


namespace core {
namespace chemical {

/// @brief  Return a list of ring conformers, read from a database file.
utility::vector1< RingConformer > read_conformers_from_database_file_for_ring_size( std::string const & filename,
		core::Size ring_size );

}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_ring_conformer_io_HH
