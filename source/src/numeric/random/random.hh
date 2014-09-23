// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random/random.hh
/// @brief  Random number generator system
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
///
/// @remarks
///  @li -
///


#ifndef INCLUDED_numeric_random_random_hh
#define INCLUDED_numeric_random_random_hh


// Unit headers
#include <numeric/random/random.fwd.hh>

// Package headers
#include <numeric/random/uniform.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <string>
#include <vector>
#include <utility/vector1.hh> // need for template definition...


///  Currently supported RG types:
///  standard - build in C++ random generator
///  ran3 - old generator from previos version of rosetta

#if (defined min) && (defined WIN32)  // Workaround for MSVC and windows.h include which used #define min
	#undef min
#endif

#if (defined max) && (defined WIN32)  // Workaround for MSVC and windows.h include which used #define max
	#undef max
#endif


namespace numeric {
namespace random {

class RandomGenerator;

/// @brief Return the one-per-thread "singleton" random generator.
RandomGenerator & rg();

/// @brief Generate a random number between 0 and 1.  Threadsafe since each thread
/// uses its own random generator.
double uniform();

/// @brief Generate a random number pulled from a standard normal -- i.e. mean of
/// zero and standard deviation of 1.  Threadsafe since each thread uses its own
/// random generator
double gaussian();

/// @brief Return a number uniformly drawn from the inclusive range between low
/// and high.  Threadsafe since each thread uses its own random generator.
int random_range(int low, int high);


/// @brief Random number generator system
class RandomGenerator : public utility::pointer::ReferenceCount
{
private:
	// Private to force the instantiation of a RandomGenerator through the rg() function.
	RandomGenerator();

	// Do not make copies of RandomGenerators - pass by reference instead.
	RandomGenerator( RandomGenerator const & ); // Unimplemented private -- do not use.

	friend RandomGenerator & rg();

public:

	~RandomGenerator();

	/// Return from range [0, 1] (?) uniform random number
	///
	/// The implementation of random_range leads me to believe this is
	/// actually [0, 1), like most other random number generators.  -IWD
	double uniform();

	/// @brief Get Gaussian distribution random number
	double gaussian();

	/// @brief Returns a random int in the range specified by the arguments
	int random_range( int low, int high );

	/// @brief Return the seed used by this RNG.
	int get_seed() const;

	/// @brief Set the seed and the generator type synchronously.
	/// Currently the two supported generator types are "standard"
	/// and "mt19937" with the latter being the recommended form.
	void set_seed( std::string const & generator_type, int seed );

	/// @brief Return the seed used by this RNG.
	void set_seed( int seed );

	void saveState(std::ostream & out);

	void restoreState(std::istream & in);

  /// @brief return a random element from a utility::vector1.  What is
	/// this function doing inside the RandomGenerator class?
	template< class T >
	T const &
	random_element( utility::vector1< T > const & v )
	{
		assert( !v.empty() );
		return v[ RandomGenerator::random_range( 1, v.size() ) ];
	}

public: // implement the boost Uniform Random Generator concept

	typedef double result_type;

	// Maybe this should call uniform() instead?
	double operator()() { return gaussian(); }

	double min() const { return 0; }
	double max() const { return 1; }

private: // Fields

	//int seed_offset_; /// Our magic number goes there

	uniform_RG_OP generator_;

	/// data for Gaussian generation
	bool   gaussian_iset_;
	double gaussian_gset_;

}; // RandomGenerator

typedef utility::pointer::shared_ptr< RandomGenerator > RandomGeneratorOP;

} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_random_HH
