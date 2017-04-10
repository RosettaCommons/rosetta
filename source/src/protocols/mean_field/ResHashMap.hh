// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    ResHashMap.hh

/// @brief   Declarations and simple accessor/mutator definitions for ResHashMap.
/// @author  arubenstein

#ifndef INCLUDED_protocols_mean_field_ResHashMap_HH
#define INCLUDED_protocols_mean_field_ResHashMap_HH

// Unit header
#include <protocols/mean_field/ResHashMap.fwd.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers


// C++ headers
#include <iostream>
#include <functional>

// Boost headers
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

namespace protocols {
namespace mean_field {

/// @brief uses default boost::hash combine to hash rotamers based on rounded chi angles
struct ResHasher : std::unary_function< core::conformation::ResidueCOP, std::size_t > {
	// use std::size_t instead of core::Size just to be
	// consistent with boost::hash types

	/// @brief return hash value given ResidueCOP
	std::size_t operator()( core::conformation::ResidueCOP res ) const {

		std::size_t seed = 0;
		boost::hash_combine( seed, res->type().aa() );
		boost::hash_combine( seed, res->name() );

		utility::vector1 < core::Real > const chi = res->chi();

		for ( core::Size nchi = 1; nchi <= chi.size(); ++nchi ) {
			core::Size chi_int = static_cast< core::Size > ( round( chi[ nchi ] / core::Real( 5.0 ) ) );
			boost::hash_combine( seed, chi_int );
		}

		return seed;
	}
};

/// @brief checks whether two rotamer hashes match - if the rounded values of their chi angles match
struct ResPred : std::binary_function< core::conformation::ResidueCOP, core::conformation::ResidueCOP, bool > {

	/// @brief return bool depending on whether two Residues hashes match
	bool operator()( core::conformation::ResidueCOP first_res, core::conformation::ResidueCOP second_res ) const {

		utility::vector1 < core::Real > const fchi = first_res->chi();
		utility::vector1 < core::Real > const schi = second_res->chi();

		if ( fchi.size() != schi.size() || first_res->aa() != second_res->aa() || first_res->name() != second_res->name() ) {
			return false;
		}

		for ( core::Size nchi = 1; nchi <= fchi.size(); ++nchi ) {
			if ( static_cast< core::Size > ( round( fchi[ nchi ]/ core::Real( 5.0 ) ) ) != static_cast< core::Size > ( round( schi[ nchi ] / core::Real( 5.0 ) ) ) ) {
				return false;
			}
		}
		return true;
	}
};

/// @details  class ResHashMap holds a hashmap of rotamers
/// @remarks
class ResHashMap : public utility::pointer::ReferenceCount {
public:

	typedef core::SSize RotamerIndex;
	typedef boost::unordered_map < core::conformation::ResidueCOP, RotamerIndex, ResHasher, ResPred > ResHash;

	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Default constructor
	ResHashMap();

	/// @brief Destructor
	~ResHashMap();


	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief Generate string representation of ResHashMap for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that ResHashMap can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, ResHashMap const & object_to_output );

	/// @brief Attempts to insert a residue res with a RotamerIndex ind into the hash map if it is not already there
	RotamerIndex attempt_insert ( core::conformation::ResidueCOP res );

	/// @brief Attempts to get the rotamer index of a residue res
	RotamerIndex get_rot_ind ( core::conformation::ResidueCOP res ) const;

	// Accessors/Mutators

	/// @brief returns the last rotamer index that was assigned by the ResHashMap
	RotamerIndex last_ind_assigned () const
	{
		return last_ind_assigned_;
	}

private:
	// Private methods /////////////////////////////////////////////////////////

	// no copy constructor
	ResHashMap( ResHashMap const & object_to_copy );

	// no Assignment operator
	ResHashMap & operator=( ResHashMap const & object_to_copy );

	// Private data ////////////////////////////////////////////////////////////
	ResHash hash_;
	RotamerIndex last_ind_assigned_;

};  // class ResHashMap

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_ResHashMap_HH
