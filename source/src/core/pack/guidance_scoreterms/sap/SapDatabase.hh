// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sap/SapDatabase.hh
/// @brief Hdf5 version of rosetta database
/// @details
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_sap_SapDatabase_hh
#define INCLUDED_core_sap_SapDatabase_hh


// Unit Headers
#include <core/pack/guidance_scoreterms/sap/SapDatabase.fwd.hh>

// Package headers

// Project headers
#include <core/conformation/Residue.fwd.hh>

// Basic headers

// Utility headers
#include <core/types.hh>
#include <utility/SingletonBase.hh>

// C++ headers
#include <unordered_map>

// forward declaration for testing
class SapConstraintEnergyTests;

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

float const SAP_BLOCK_STORE_SCALE = 2.2;

float const MAX_ATOMIC_RADIUS = 2.5;
float const MAX_BLOCK_INTERACTION = MAX_ATOMIC_RADIUS*2 + 3;

Real const SAP_PROBE_SIZE = 1.1;


struct BlockParam {
	BlockParam() {}

	float max_sasa_score;
	uint8_t no_block;
	uint8_t full_block;    // This should be no more than 255 because of how we stored blocks

	BlockParam( float max, uint8_t no, uint8_t full )
	:
		max_sasa_score( max ),
		no_block( no ),
		full_block( full )
	{}
};


class SapDatabase : public utility::SingletonBase< SapDatabase > {
public:
	friend class utility::SingletonBase< SapDatabase >;
	friend class ::SapConstraintEnergyTests;

	utility::pointer::shared_ptr< std::unordered_map< std::string, BlockParam > const >
	atomtype_to_block_param() const;

	Real
	hydrophobic_weight( char aa ) const;

	Real
	max_sasa( char aa ) const;

	std::pair< char, std::string >
	get_name1_name3( core::conformation::Residue const & res, bool warn ) const;

	bool symm_debug() const;
	bool symm_debug_force_map() const;


private:

	SapDatabase();

	void
	load_block_data();

	void
	load_hydrophobic_data();

	void
	generate_max_sasa();

private:

	// Unlike the BlockParams in SapConstraintHelper, the max_sasa_score here is actually just the max atom sasa
	utility::pointer::shared_ptr< std::unordered_map< std::string, BlockParam > > atomtype_to_block_param_;

	std::unordered_map< char, Real > name1_to_hydrophobic_;
	std::unordered_map< char, Real > name1_to_max_sasa_;

	// This shouldn't be here but it makes the test suite a whole lot easier to write
	bool symm_debug_;
	bool symm_debug_force_map_;

};


} //sap
} //guidance_scoreterms
} //pack
} //core


#endif
