// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/RingConformerManager.hh
/// @brief   Declarations and simple accessor/mutator definitions for RingConformerManager.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_rings_RingConformerManager_HH
#define INCLUDED_core_chemical_rings_RingConformerManager_HH

// Unit header
// No fwd.hh of RingConformer exists
namespace core { namespace chemical { namespace rings { struct RingConformer; } } }
//#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerManager.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>


namespace core {
namespace chemical {
namespace rings {

/// @details  This class is a singleton and manages RingConformer data that should only be read from the database one
/// time and shared among all RingConformerSets.
class RingConformerManager : public utility::SingletonBase< RingConformerManager > {
public:  // Declare friends ///////////////////////////////////////////////////
	friend class utility::SingletonBase< RingConformerManager >;


public:  // Static constant data access ///////////////////////////////////////
	/// @brief  Return a set of ring conformers for the requested ring size.
	static utility::vector1< RingConformer > const & conformers_for_ring_size( core::Size ring_size );


private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	RingConformerManager();

	// Get the conformers requested, creating them if necessary.
	// Called by the public static method conformers_for_ring_size().
	utility::vector1< RingConformer > const & get_conformers_for_ring_size( core::Size ring_size );


private:  // Private data /////////////////////////////////////////////////////
	std::map< core::uint, utility::vector1< RingConformer > > conformers_;
};

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_RingConformerManager_HH
