// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/rings/RingConformerSet.hh
/// @brief   Declarations and simple accessor/mutator definitions for RingConformerSet.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_rings_RingConformerSet_HH
#define INCLUDED_core_chemical_rings_RingConformerSet_HH

// Unit header
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerSet.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <iostream>


namespace core {
namespace chemical {
namespace rings {

/// @brief  Enumerators for the three Cremer-Pople "ring-puckering" parameters used to describe 4-, 5-, and 6-membered
/// ring conformers
enum CPParameter {
	q = 1,  // used to describe all ring sizes greater than 3
	PHI,    // used to describe all ring sizes greater than 4
	THETA   // used to describe all ring sizes greater than 5
};


class RingConformerSet : public utility::pointer::ReferenceCount {
public:
	// Standard methods //////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Standard constructor
	RingConformerSet(
		core::Size const ring_size,
		std::string const & lowest_conformer,
		utility::vector1< std::string > const & low_conformers );

	/// @brief  Copy constructor
	RingConformerSet( RingConformerSet const & object_to_copy );

	// Assignment operator
	RingConformerSet & operator=( RingConformerSet const & object_to_copy );

	// Destructor
	virtual ~RingConformerSet();


	// Standard Rosetta methods //////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Generate string representation of RingConformerSet for debugging purposes.
	virtual void show( std::ostream & output=std::cout ) const;


	// Accessors/Mutators ////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Return the ring size of the conformers in this set.
	core::Size ring_size() const { return ring_size_; }

	/// @brief  Return the size of the conformer set.
	core::Size size() const;

	/// @brief  Are the low-energy conformers known for this set?
	bool low_energy_conformers_are_known() const;


	/// @brief  Return a list of all nondegenerate conformers in the set.
	utility::vector1< RingConformer > const & get_all_nondegenerate_conformers() const;

	// AMW: cppcheck wants you to change to pass by reference; DO NOT
	// Why not? ~Labonte
	/// @brief  Return the conformer corresponding to the requested name.
	RingConformer const & get_ideal_conformer_by_name( std::string const & name ) const;

	// AMW: cppcheck wants you to change to pass by reference; DO NOT
	// Why not? ~Labonte
	/// @brief  Return the conformer that is the best fit for the provided Cremer-Pople parameters.
	RingConformer const & get_ideal_conformer_by_CP_parameters( utility::vector1< core::Real > const & parameters ) const;


	/// @brief  Return the conformer that is the best fit for the provided list of nu angles.
	RingConformer const & get_ideal_conformer_from_nus( utility::vector1< core::Angle > const & angles ) const;


	/// @brief  Return the conformer that is known from studies (if available) to be the lowest energy ring conformer.
	RingConformer const & get_lowest_energy_conformer() const;

	/// @brief  Return a random conformer from the set.
	RingConformer const & get_random_conformer() const;

	/// @brief  Return a random conformer from the subset of conformers that are local minima.
	RingConformer const & get_random_local_min_conformer() const;

	utility::vector1< RingConformer > get_local_min_conformers() const;

private:
	// Private methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Empty constructor
	RingConformerSet();

	// Initialize data members for the given ring size.
	void init(
		core::Size const ring_size,
		std::string const & lowest_conformer_in,
		utility::vector1< std::string > const & low_conformers );

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( RingConformerSet & object_to_copy_to, RingConformerSet const & object_to_copy_from );


	// Private data //////////////////////////////////////////////////////////////////////////////////////////////////
	core::Size ring_size_;  // almost always 5 or 6, but one could make a RingConformerSet for other sizes

	// Ring Conformer Subsets
	utility::vector1< RingConformer > nondegenerate_conformers_;
	utility::vector1< RingConformer > degenerate_conformers_;  // includes multiple copies of degenerate conformers
	RingConformer energy_minimum_conformer_;  // the global minimum conformer
	utility::vector1< RingConformer > energy_minima_conformers_;  // other relatively stable ring conformers
};  // class RingConformerSet


// Insertion operators (overloaded so that RingConformer and RingConformerSet can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, RingConformer const & object_to_output );

std::ostream & operator<<( std::ostream & output, RingConformerSet const & object_to_output );

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_RingConformerSet_HH
