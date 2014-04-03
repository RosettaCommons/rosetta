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
#include <numeric/random/uniform.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <iostream>
#include <string>
#include <vector>
#include <utility/vector1.hh> // need for template definition...


///  Currently supported RG types:
///  standard - build in C++ random generator
///  ran3 - old generator from previos version of rosetta

namespace numeric {
namespace random {

class RandomGenerator;


// Single static random generator functions/object for convenience
// Use them only for benchmark unrelated tasks.
extern RandomGenerator RG;
double uniform(void);
double gaussian(void);
int random_range(int low, int high);

// Forward
class uniform_RG;

/// @brief Different type of initialization scheme:
/// _RND_NormalRun_ - all instances of RandomGenerator will be pointing to just one
///   generator. This is main production mode. We do not use multiple generators here
///   because putting many generators will effectively produce long-range correlations
///   between random numbers obtained from different places of program. (By
///   using N generators instead of one we effectively reducing dimensionality
///   of random generator by at least N. Now take in to account number of
///   nodes in computer cluster and its look like we can quickly run in to
///   trouble even with 600 dimensional Mersenne twister.)
///   Somewhat helpful presentation on this topic: "Don't Trust Parallel Monte Carlo!" by
///   Peter Hellekalek, available online at: http://random.mat.sbg.ac.at/~peter/pads98.ps
///
/// _RND_TestRun_ - each instance of RandomGenerator will have its own generator with its
///   own seed. We using this approach because we trying to increase stability of
///   unit test / performance tests / scientific test.
///   WARNING: This mode design specifically for testing purpose and SHOULD NOT be used
///   in production environment (see comment above).
enum RND_RunType { _RND_NormalRun_, _RND_TestRun_ };

/// @brief Random number generator system
class RandomGenerator
{
private:
	// RandomGenerators must be initialized with a magicNumber.
	RandomGenerator(); // Unimplemented private -- do not use.
	// Do not make copies of RandomGenerators - pass by reference instead.
	RandomGenerator( RandomGenerator const & ); // Unimplemented private -- do not use.
public:
	RandomGenerator(int const magicNumber);
	~RandomGenerator();

	/// Return from range [0, 1] (?) uniform random number
	///
	/// The implementation of random_range leads me to believe this is
	/// actually [0, 1), like most other random number generators.  -IWD
	inline double uniform() { return generator->getRandom(); }

	/// @brief Get Gaussian distribution random number
	double gaussian();

	/// @brief Returns a random int in the range specified by the arguments
	int random_range( int low, int high );

	/// @brief Return the seed used by this RNG.
	int get_seed() { return generator->getSeed(); }

	/// @brief Return the seed used by this RNG.
	void set_seed(int seed) { generator->setSeed(seed); }

	void saveState(std::ostream & out);

	void restoreState(std::istream & in);

		// Static functions
	/// init all rundom number generators in program, must be called after main()
	///  start executing
	static void	initializeRandomGenerators(
		int const start_seed,
		RND_RunType run_type,
		std::string const & type = ""
	);

	///@brief Saves the state of all random number generators to given stream.
	static void saveAllStates(std::ostream & out);

	///@brief Restores the state of all random number generators from given stream.
	static void restoreAllStates(std::istream & in);

  /// return a random element from a utility::vector1
	template< class T >
	T const &
	random_element( utility::vector1< T > const & v )
	{
		assert( !v.empty() );
		return v[ RandomGenerator::random_range( 1, v.size() ) ];
	}

public: // implement the boost Uniform Random Generator concept

	typedef double result_type;

	double operator()() { return gaussian(); }

	double min() const { return 0; }
	double max() const { return 1; }

private: // Fields
	int seed_offset; /// Our magic number goes there

	utility::pointer::owning_ptr<uniform_RG> generator;

	/// flags for Gaussian generation
	bool gaussian_iset;
	double gaussian_gset;


	static std::vector<RandomGenerator*> &allGenerators();

	static utility::pointer::owning_ptr<uniform_RG> createIntNumberGenerator(
		std::string const & type
	);
}; // RandomGenerator


} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_random_HH
