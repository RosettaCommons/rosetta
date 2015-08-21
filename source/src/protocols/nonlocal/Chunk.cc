// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/Chunk.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/Chunk.hh>

// External headers
#include <boost/scoped_ptr.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/math/distributions/normal.hpp>

// Utility headers
#include <numeric/random/DistributionSampler.hh>
#include <utility/pointer/ReferenceCount.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>

// Package headers
#include <protocols/nonlocal/Region.hh>

namespace protocols {
namespace nonlocal {

// @brief Auto-generated virtual destructor
Chunk::~Chunk() {}

/// @brief A small, positive value added to the standard deviation of each chunk.
/// The distribution requires its standard deviation to be strictly positive.
#define SALT 0.00000001

/// @brief A multiplicative factor applied to control the width of the normal
/// distribution. Affects the width of the distribution.
#define SD_MULTIPLIER 1.0

typedef boost::math::normal Normal;
using core::Size;
using core::kinematics::MoveMapCOP;
using utility::pointer::ReferenceCount;

// -- Construction and Assignment -- //
//
// Remember: scoped_ptr's aren't copyable. We need to allocate new memory for
// <sampler_> in the copy constructor and assignment operator.

Chunk::Chunk(const RegionOP& region, const MoveMapOP& movable)
: ReferenceCount(), region_(region), movable_(movable) {
	using namespace boost::accumulators;

	// compute normal distribution parameters
	accumulator_set<Size, stats<tag::mean, tag::variance> > acc;
	for ( Size i = start(); i <= stop(); ++i ) {
		acc(i);
	}

	double mu = mean(acc);
	double sigma = (std::sqrt(variance(acc)) * SD_MULTIPLIER) + SALT;
	Normal dist(mu, sigma);
	sampler_.reset(new numeric::random::DistributionSampler<Normal>(dist));
}

Chunk::Chunk(const Chunk& other)
: ReferenceCount(), region_(other.region_), movable_(other.movable_) {
	Normal dist = other.sampler_->distribution();
	sampler_.reset(new numeric::random::DistributionSampler<Normal>(dist));
}

Chunk& Chunk::operator=(const Chunk& other) {
	// check for self-assignment
	if ( this == &other ) {
		return *this;
	}

	region_ = other.region_;
	movable_ = other.movable_;
	Normal dist = other.sampler_->distribution();
	sampler_.reset(new numeric::random::DistributionSampler<Normal>(dist));
	return *this;
}

// -- Accessors -- //

Size Chunk::choose() const {
	assert(valid());
	while ( 1 ) {
		Size insert_pos = static_cast<Size>(sampler_->sample());

		// restrict samples to the closed interval [start(), stop()]
		if ( insert_pos < start() || insert_pos > stop() ) {
			continue;
		}

		if ( movable_->get_bb(insert_pos) ) {
			return insert_pos;
		}
	}
	return 0;
}

Size Chunk::start() const {
	return region_->start();
}

Size Chunk::stop() const {
	return region_->stop();
}

Size Chunk::length() const {
	return region_->length();
}

// -- Utility Functions -- //

bool Chunk::is_movable() const {
	for ( Size i = start(); i <= stop(); ++i ) {
		if ( movable_->get_bb(i) ) {
			return true;
		}
	}

	return false;
}

bool Chunk::valid() const {
	// Region's are agnostic as to their direction: start => stop, stop <= start.
	// For simplicity, our iteration methods assume a left-to-right orientation
	if ( start() >= stop() ) {
		return false;
	}

	// In order to avoid an infinite loop in choose(), at least one position on
	// the interval [start(), stop()] must be movable. Short-circuit evaluation.
	return is_movable();
}

}  // namespace nonlocal
}  // namespace protocols
