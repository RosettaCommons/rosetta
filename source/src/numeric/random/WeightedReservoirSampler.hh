// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/random/WeightedReservoirSampler.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_NUMERIC_RANDOM_WEIGHTEDRESERVOIRSAMPLER_hh
#define INCLUDED_NUMERIC_RANDOM_WEIGHTEDRESERVOIRSAMPLER_hh

// C/C++ headers
#include <cmath>
#include <queue>

// External headers
#include <boost/utility.hpp>

// Utility headers
#include <utility/vector1.hh>

// Package headers
#include <numeric/random/random.hh>

namespace numeric {
namespace random {

/// @class A simple class for associating a real-valued weight with an item.
/// Overloads operator< for use in sorted containers.
template <typename T>
class WeightedReservoirItem {
public:
	/// @brief Associates the given item with a real-valued weight
	WeightedReservoirItem(T item, double weight) : item_(item), weight_(weight) {}

	/// @brief Returns the item
	T item() const {
		return item_;
	}

	/// @brief Returns the weight
	double weight() const {
		return weight_;
	}

	/// @brief Returns true if this item's weight is less than <o>'s weight
	bool operator<( WeightedReservoirItem<T> const & o ) const {
		return weight() < o.weight();
	}

private:
	/// @brief the item
	T item_;

	/// @brief the item's weight
	double weight_;
};

/// @class A library for weighted reservoir sampling. Retrieves M samples from a
/// population of N items in time O(N) and space O(M).
///
/// Requirements:
///   - T::operator< is defined
///   - Weights are strictly positive
template <typename T>
class WeightedReservoirSampler : boost::noncopyable {
public:
	/// @brief Constructs a new weighted reservoir sampler with the given capacity
	explicit WeightedReservoirSampler(int capacity)
	: capacity_(capacity), num_considered_(0) {
		assert(capacity > 0);
	}

	/// @brief Considers the given item for inclusion in the reservoir.
	/// Items with non-positive fitnesses are not considered.
	void consider_sample(T item, double fitness) {
		if ( fitness <= 0 ) {
			return;
		}

		++num_considered_;
		double weight = -std::log(uniform()) / fitness;

		if ( num_considered() <= capacity() ) {
			reservoir_.push(WeightedReservoirItem<T>(item, weight));
		} else {
			if ( weight < reservoir_.top().weight() ) {
				reservoir_.pop();
				reservoir_.push(WeightedReservoirItem<T>(item, weight));
			}
		}
	}

	/// @brief Populates <selected> with the contents of the reservoir.
	///
	/// TODO(cmiles) use a data structure that provides iterator access.
	/// As a result of calling this method, the reservoir is emptied and the
	/// sampler's state restored to its initial condition.
	void samples(utility::vector1<T>* selected) {
		assert(selected);
		while ( !reservoir_.empty() ) {
			const WeightedReservoirItem<T>& entry = reservoir_.top();
			selected->push_back(entry.item());
			reservoir_.pop();
		}
		reset();
	}

	/// @brief Restores the reservoir to its initial state
	void reset() {
		num_considered_ = 0;
		while ( !reservoir_.empty() ) {
			reservoir_.pop();
		}
	}

	/// @brief Returns the number of items considered by the sampler
	unsigned long num_considered() const {
		return num_considered_;
	}

	/// @brief Returns the capacity of the reservoir (i.e. number of samples)
	unsigned long capacity() const {
		return capacity_;
	}

private:
	/// @brief Capacity of the reservoir
	unsigned long capacity_;

	/// @brief Number of samples considered for inclusion
	unsigned long num_considered_;

	/// @brief Maintains the current set of samples selected from the population
	std::priority_queue<WeightedReservoirItem<T> > reservoir_;
};

}  // namespace random
}  // namespace numeric

#endif  // INCLUDED_NUMERIC_RANDOM_WEIGHTEDRESERVOIRSAMPLER_hh
