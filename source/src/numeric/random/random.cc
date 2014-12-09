// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/random/random.cc
/// @brief  Random number generator system
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

// Unit headers
#include <numeric/random/random.hh>

// Package headers
#include <numeric/random/uniform.hh>
#include <numeric/random/mt19937.hh>

// Utility headers
#include <utility/exit.hh>

// C++ 11 Compatibility
#include <utility/thread/backwards_thread_local.hh>

// C++ headers
#include <iostream>
#include <cmath>

// C++ 11 support
#include <utility/thread/backwards_thread_local.hh>


namespace numeric {
namespace random {

#ifdef CXX11
// Use a "thread_local" declaration to ensure that the each thread has its own
// random_generator pointer.
thread_local RandomGeneratorOP random_generator( 0 );
#else
static RandomGeneratorOP random_generator( 0 );
#endif

RandomGenerator & rg()
{
	// If you want thread saftety mechanisms but don't want to use C++11 (foldit?) then
	// you'll want to edit this logic.
	if ( random_generator == 0 ) {
		random_generator = RandomGeneratorOP( new RandomGenerator );
	}
	return *random_generator;
}

double uniform() { return rg().uniform(); }
double gaussian() { return rg().gaussian(); }
int random_range(int low, int high) { return rg().random_range(low, high); }


using namespace std;

/// uniform_RG factory
uniform_RG_OP
createRG(string const & type )
{
	if( type == "standard" ) return uniform_RG_OP( new standard_RG() );
	//if( type == "ran3" ) return new ran3_RG();
	if( type == "mt19937" ) return uniform_RG_OP( new mt19937_RG() );

	utility_exit_with_message("Unknown random number generator type: " + type);
	return 0;
}

RandomGenerator::RandomGenerator() :
	gaussian_iset_( true ),
	gaussian_gset_( 0 )
{}

RandomGenerator::~RandomGenerator()
{}


double
RandomGenerator::uniform() {
	return generator_->getRandom();
}

/// @brief
/// SL: this is function is ported from old Rosetta++.
///
/// Returns a gaussian random number (normally distributed deviate with
/// zero mean and unit variance) using ran3 as a source of uniform deviates.
/// Always call with the same idum
///
/// @references  Numerical Recipes, section 7.2, a.k.a. "GASDEV"
///
/// @author JJG 4/01
double RandomGenerator::gaussian()
{
	double v1, v2, rsq, fac;

	double rgaussian; // Return value
	if ( gaussian_iset_ ) {
		do {
			v1 = 2.0f * uniform() - 1.0f;
			v2 = 2.0f * uniform() - 1.0f;
			rsq = ( v1 * v1 ) + ( v2 * v2 );
		} while ( rsq >= 1.0 || rsq == 0.0 );
		fac = std::sqrt(-(2.0*std::log(rsq)/rsq));
		gaussian_gset_ = v1*fac;
		rgaussian = v2*fac;
		gaussian_iset_ = false;
	} else {
		rgaussian = gaussian_gset_;
		gaussian_iset_ = true;
	}
//	std::cout << "NR Gaussian: " << gaussian << std::endl;
	return rgaussian;
}


/// @brief Returns a random int in the range specified by the arguments,
/// with both enpoints being included in the possible output.
///
/// @author XA
int RandomGenerator::random_range(int low, int high)
{
	if ( low > high ) {
		int temp;
		temp = low;
		low = high;
		high = temp;
	}

	int const range( high - low + 1 );
	return ( static_cast< int >( range * numeric::random::uniform() ) + low );
}

int RandomGenerator::get_seed() const
{
	return generator_->getSeed();
}

void
RandomGenerator::set_seed(
	std::string const & generator_type,
	int seed
)
{
	generator_ = createRG( generator_type );
	generator_->setSeed( seed );
}


void RandomGenerator::set_seed( int seed ) {
	generator_->setSeed(seed);
}

void RandomGenerator::saveState(std::ostream & out)
{
	//out << " " << seed_offset;
	out << " " << gaussian_iset_;
	// Cast raw bits to int for perfect precision:
	out << " " << *((uint64_t*) &gaussian_gset_);
	generator_->saveState(out);
	out << "\n";
}

void RandomGenerator::restoreState(std::istream & in)
{
	in >> gaussian_iset_;
	// Cast raw bits from int for perfect precision:
	uint64_t gset_bits = 0;
	in >> gset_bits;
	gaussian_gset_ = *((double*) &gset_bits);
	generator_->restoreState(in);
}

} // namespace random
} // namespace numeric
