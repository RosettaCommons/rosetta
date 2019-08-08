// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/RingConformerManager.cc
/// @brief   Method definitions for RingConformerManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerManager.hh>
#include <core/chemical/rings/ring_conformer_io.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/database/open.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rings.OptionKeys.gen.hh>

// C++ headers
#include <map>
#include <sstream>


namespace core {
namespace chemical {
namespace rings {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
utility::vector1< RingConformer > const &
RingConformerManager::conformers_for_ring_size_and_type(
	core::Size const ring_size,
	core::chemical::rings::RingSaturationType const type )
{
	return get_instance()->get_conformers_for_ring_size_and_type( ring_size, type );
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RingConformerManager::RingConformerManager() = default;

// Get the conformers requested, creating them if necessary.
// Called by the public static method conformers_for_ring_size().
utility::vector1< RingConformer > const &
RingConformerManager::get_conformers_for_ring_size_and_type(
	core::Size const ring_size,
	core::chemical::rings::RingSaturationType const type )
{
	using namespace std;
	using namespace utility;
	using namespace basic::options;

	// Only create sets one time, as needed, for each ring size and saturation type.
	if ( ! conformers_.count( make_pair( ring_size, type ) ) ) {
		stringstream filename( stringstream::out );
		filename << option[ OptionKeys::rings::ring_conformer_dbpath ]();
		filename << "/" << ring_size << "-membered_";
		// TODO: Switch to a switch or function when other saturation types are added.
		filename << ( ( type == ALIPHATIC ) ? "aliphatic" : "aromatic" ) << "_ring_conformers.data";
		vector1< RingConformer > conformers( read_conformers_from_database_file_for_ring_size(
			basic::database::full_name( filename.str() ), ring_size ) );
		conformers_.insert( make_pair( make_pair( ring_size, type ), conformers ) );
	}
	return conformers_[ make_pair( ring_size, type ) ];
}

}  // namespace rings
}  // namespace chemical
}  // namespace core
