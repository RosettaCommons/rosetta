// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/uniform.hh
/// @brief  Uniform random number generator
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
///
/// @remarks
///  @li -


#ifndef INCLUDED_numeric_random_uniform_hh
#define INCLUDED_numeric_random_uniform_hh

#include <numeric/random/uniform.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <string>
#include <cstdlib> //required by GCC 4.3.2

#include <cassert>

namespace numeric {
namespace random {


/// @brief Uniform random number generator
class uniform_RG : public utility::pointer::ReferenceCount
{
public:
	uniform_RG() {}

	~uniform_RG() override = default;

	virtual void setSeed(int const seed) = 0;

	virtual void setSeed(std::string const & seed) = 0;  //< some generator take more than int to init

	virtual int getSeed() = 0;

	virtual double getRandom() = 0;

	virtual void saveState(std::ostream & out) = 0;

	virtual void restoreState(std::istream & in) = 0;

}; // uniform_RG


/// @brief Generator based on rand() < clib > function.
class standard_RG : public uniform_RG
{
public:

	inline standard_RG():
		seed_(0) // Clear, for consistent starting state
	{}

	inline ~standard_RG() override = default;

	inline void setSeed(int const seed) override { seed_ = seed; srand( seed_ ); }

	inline void setSeed(std::string const &) override { assert( false ); } // Not implemented yet!

	inline int getSeed() override { return seed_; }

	inline double getRandom() override { return (double)rand() / (double)RAND_MAX; }

	void saveState(std::ostream & /*out*/) override { assert( false ); } // Not implemented yet!

	void restoreState(std::istream & /*in*/) override { assert( false ); } // Not implemented yet!

private:
	int seed_;

}; // standard_RG


} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_uniform_HH
