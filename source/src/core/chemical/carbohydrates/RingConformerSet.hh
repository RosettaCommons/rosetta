// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/RingConformerSet.hh
/// @brief   Declarations and simple accessor/mutator definitions for RingConformerSet.
/// @author  Labonte

#ifndef INCLUDED_core_chemical_carbohydrates_RingConformerSet_HH
#define INCLUDED_core_chemical_carbohydrates_RingConformerSet_HH

// Unit header
#include <core/chemical/carbohydrates/RingConformerSet.fwd.hh>

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
namespace carbohydrates {

struct RingConformer : public utility::pointer::ReferenceCount {
	std::string specific_name;  // e.g., "1C4"
	std::string general_name;  // e.g., "chair"

	core::uint degeneracy; // E.g., 1C4 has a degeneracy of 3, since 3CO and 5C2 are equivalent.

	// a list of the torsion angles for the 1st n-2 nu angles, where n is the ring size.
	utility::vector1<core::Angle> ideal_angles;
};  // struct RingConformer


class RingConformerSet : public utility::pointer::ReferenceCount {
public:
	// Standard methods //////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Standard constructor
	RingConformerSet(core::uint ring_size);

	/// @brief  Copy constructor
	RingConformerSet(RingConformerSet const & object_to_copy);

	// Assignment operator
	RingConformerSet & operator=(RingConformerSet const & object_to_copy);

	// Destructor
	~RingConformerSet();


	// Standard Rosetta methods //////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Generate string representation of RingConformerSet for debugging purposes.
	virtual void show(std::ostream & output=std::cout) const;


	// Accessors/Mutators ////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Return a list of all nondegenerate conformers in the set.
	utility::vector1<RingConformerCOP> get_all_nondegenerate_conformers() const;

	/// @brief  Return the conformer that is the best fit for the provided list of nu angles.
	RingConformerCOP get_conformer_from_nus(utility::vector1<core::Angle> angles) const;

	/// @brief  Return the conformer that is known from studies (if available) to be the lowest energy ring conformer.
	RingConformerCOP get_lowest_energy_conformer() const;

	/// @brief  Return a random conformer from the set.
	RingConformerCOP get_random_conformer() const;

	/// @brief  Return a random conformer from the subset of conformers that are local minima.
	// TODO: better?: overload get_random_conformer and pass enum, such as "LOCAL_MIN"
	RingConformerCOP get_random_local_min_conformer() const;


private:
	// Private methods ///////////////////////////////////////////////////////////////////////////////////////////////
	// Empty constructor
	RingConformerSet();

	// Initialize data members for the given ring size.
	void init(core::uint ring_size);

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(RingConformerSet object_to_copy_to, RingConformerSet object_to_copy_from);


	// Static constant data access
	/// @brief A set of ring conformers for the requested ring size.
	static utility::vector1<RingConformer> const & conformers_for_ring_size(core::Size ring_size);


	// Private data //////////////////////////////////////////////////////////////////////////////////////////////////
	core::Size ring_size_;  // almost always 5 or 6, but one could make a RingConformerSet for other sizes
	// TODO: Should this data be const?  Can it be const?  The same applies to CarbohydrateInfo data.  I don't think
	// it can be const unless it is set in the initializer, but I am not using the initializer....

	// TODO: Make these map<uint, vector1<RingConformer> >s, with the vectors indexed by an enum?
	// Ring Conformer Subsets
	utility::vector1<RingConformerCOP> nondegenerate_conformers_;
	utility::vector1<RingConformerCOP> degenerate_conformers_;  // includes multiple copies of degenerate conformers
	utility::vector1<RingConformerCOP> energy_minima_conformers_;
	utility::vector1<RingConformerCOP> energy_maxima_conformers_;

};  // class RingConformerSet

// Insertion operator (overloaded so that RingConformerSet can be "printed" in PyRosetta).
std::ostream & operator<<(std::ostream & output, RingConformerSet const & object_to_output);


}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_RingConformerSet_HH
