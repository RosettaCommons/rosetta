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
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <iostream>
#include <sstream>


// Construct tracer.
static basic::Tracer TR("core.chemical.carbohydrates.RingConformerSet");

// Construct random-number generator.
static numeric::random::RandomGenerator RG(28);  // the 2nd perfect number


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

	output << "Possible ring conformers:" << endl;

	Size n_conformers = nondegenerate_conformers_.size();
	for (uint i = 1; i <= n_conformers; ++i) {
		output << "   " << nondegenerate_conformers_[i]->specific_name << endl;
	}

	output << endl;
}


// Accessors/Mutators
// Return a list of all nondegenerate conformers in the set.
utility::vector1<RingConformerCOP>
RingConformerSet::get_all_nondegenerate_conformers() const
{
	return nondegenerate_conformers_;
}

// Return the conformer that is the best fit for the provided list of nu angles.
RingConformerCOP
RingConformerSet::get_conformer_from_nus(utility::vector1<core::Angle> /*angles*/) const
{
	RingConformerOP conformer;

	return conformer;
}

// Return the conformer that is known from studies (if available) to be the lowest energy ring conformer.
RingConformerCOP
RingConformerSet::get_lowest_energy_conformer() const
{
	RingConformerOP conformer;

	return conformer;
}

// Return a random conformer from the set.
RingConformerCOP
RingConformerSet::get_random_conformer() const
{
	uint i = uint(RG.uniform() * degenerate_conformers_.size() + 1);

	return degenerate_conformers_[i];
}

// Return a random conformer from the subset of conformers that are local minima.
// TODO: better?: overload get_random_conformer and pass enum, such as "LOCAL_MIN"
RingConformerCOP
RingConformerSet::get_random_local_min_conformer() const
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

	Size n_conformers = conformers_for_ring_size(ring_size_).size();
	for (uint i = 1; i <= n_conformers; ++i) {
		RingConformerCOP conformer = &conformers_for_ring_size(ring_size_).at(i);
		nondegenerate_conformers_.push_back(conformer);
		for (uint copy_num = 1; copy_num <= conformer->degeneracy; ++copy_num) {
			degenerate_conformers_.push_back(conformer);
		}

		// Other subsets
		// TODO:
		// if (CarbohydrateInfo->vector_of_names_from_params_file.contains(conformer->specific_name))...
		// will need CAP back to owning CarbohydrateInfo or eventually ResidueType
		//energy_minima_conformers_ = vector1<RingConformerCOP>();  // TEMP
		//energy_maxima_conformers_ = vector1<RingConformerCOP>();  // TEMP
	}
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingConformerSet::copy_data(
		RingConformerSet object_to_copy_to,
		RingConformerSet object_to_copy_from)
{
	object_to_copy_to.ring_size_ = object_to_copy_from.ring_size_;
	object_to_copy_to.nondegenerate_conformers_ = object_to_copy_from.nondegenerate_conformers_;
	object_to_copy_to.degenerate_conformers_ = object_to_copy_from.degenerate_conformers_;
	//object_to_copy_to.energy_minima_conformers_ = object_to_copy_from.energy_minima_conformers_;
	//object_to_copy_to.energy_maxima_conformers_ = object_to_copy_from.energy_maxima_conformers_;
}


// Static constant data access
utility::vector1<RingConformer> const &
RingConformerSet::conformers_for_ring_size(core::Size ring_size)
{
	using namespace std;
	using namespace utility;

	static map<uint, vector1<RingConformer> > *CONFORMERS = NULL;

	// If statement ensures that the data is only created once, i.e., is constant.
	if (!CONFORMERS) {
		// A map of sets of ring conformers keyed by ring size
		CONFORMERS = new map<uint, vector1<RingConformer> >;
	}

	// Only create sets one time, as needed, for each ring size.
	if (!CONFORMERS->count(ring_size)) {
		stringstream filename(stringstream::out);
		filename << "chemical/carbohydrates/" << ring_size << "-membered_ring_conformers.data";
		vector1<RingConformer> conformers = read_conformers_from_database_file_for_ring_size(
				basic::options::option[basic::options::OptionKeys::in::path::database](1).name() +
				filename.str(), ring_size);
		CONFORMERS->insert(make_pair(ring_size, conformers));
	}

	return CONFORMERS->at(ring_size);
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
