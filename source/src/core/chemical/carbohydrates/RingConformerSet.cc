// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/RingConformerSet.cc
/// @brief   Method definitions for RingConformerSet.
/// @author  Labonte

// Unit header
#include <core/chemical/carbohydrates/RingConformerSet.hh>
#include <core/chemical/carbohydrates/database_io.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <iostream>

// Construct tracer.
static basic::Tracer TR("core.chemical.carbohydrates.RingConformerSet");


namespace core {
namespace chemical {
namespace carbohydrates {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////////////////////////////////////////////
// Standard constructor
/// @param    <ring_size>: an unsigned integer expressing the size of the saccharide ring
RingConformerSet::RingConformerSet(core::uint ring_size) : utility::pointer::ReferenceCount()
{
	init(ring_size);
}

// Copy constructor
RingConformerSet::RingConformerSet(RingConformerSet const & object_to_copy) :
		utility::pointer::ReferenceCount(object_to_copy)
{
	copy_data(*this, object_to_copy);
}

// Assignment operator
RingConformerSet &
RingConformerSet::operator=(RingConformerSet const & object_to_copy)
{
	// Abort self-assignment.
	if (this == &object_to_copy) {
		return *this;
	}

	copy_data(*this, object_to_copy);
	return *this;
}

// Destructor
RingConformerSet::~RingConformerSet() {}


// Standard Rosetta methods ///////////////////////////////////////////////////////////////////////////////////////////
// General methods
void
RingConformerSet::show(std::ostream & output) const
{
	using namespace std;

	output << endl;
}


// Accessors/Mutators
// Return the conformer that is the best fit for the provided list of nu angles.
RingConformerOP
RingConformerSet::get_conformer_from_nus(utility::vector1<core::Angle> angles)
{
	RingConformerOP conformer;

	return conformer;
}

// Return the conformer that is known from studies (if available) to be the lowest energy ring conformer.
RingConformerOP
RingConformerSet::get_lowest_energy_conformer()
{
	RingConformerOP conformer;

	return conformer;
}

// Return a random conformer from the set.
RingConformerOP
RingConformerSet::get_random_conformer()
{
	RingConformerOP conformer;

	return conformer;
}

// Return a random conformer from the subset of conformers that are local minima.
// TODO: better?: overload get_random_conformer and pass enum, such as "LOCAL_MIN"
RingConformerOP
RingConformerSet::get_random_local_min_conformer()
{
	RingConformerOP conformer;

	return conformer;
}


// Private methods ////////////////////////////////////////////////////////////////////////////////////////////////////
// Empty constructor
RingConformerSet::RingConformerSet() : utility::pointer::ReferenceCount()
{
	init(6);  // default to creating a conformer set for six-membered rings
}

// Initialize data members for the given ring size.
void
RingConformerSet::init(core::uint ring_size)
{
	using namespace utility;

	ring_size_ = ring_size;

	// TEMP
	RingConformer conformer;
	conformer.specific_name = "1C4";
	conformer.general_name = "chair";
	Angle angles[4];
	angles[0] = -60.0; angles[1] = 60.0; angles[2] = -60.0; angles[3] = 60.0;
	conformer.ideal_angles = vector1<Angle>(angles, angles + 4);
	conformers_.push_back(&conformer);
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingConformerSet::copy_data(
		RingConformerSet object_to_copy_to,
		RingConformerSet object_to_copy_from)
{
	object_to_copy_to.ring_size_ = object_to_copy_from.ring_size_;
	object_to_copy_to.conformers_ = object_to_copy_from.conformers_;
}


// Helper methods /////////////////////////////////////////////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that RingConformerSet can be "printed" in PyRosetta).
std::ostream &
operator<<(std::ostream & output, RingConformerSet const & object_to_output)
{
	object_to_output.show(output);
	return output;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
