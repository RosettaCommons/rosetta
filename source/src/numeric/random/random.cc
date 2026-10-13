// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

// Numeric headers
#include <numeric/constants.hh>

// C++ headers
#include <iostream>
#include <cmath>

// C++ 11 support


namespace numeric {
namespace random {

// Use a "thread_local" declaration to ensure that the each thread has its own
// random_generator pointer.
static THREAD_LOCAL RandomGeneratorOP random_generator( nullptr );

RandomGenerator & rg()
{
	// If you want thread saftety mechanisms but don't want to use C++11 (foldit?) then
	// you'll want to edit this logic.
	if ( random_generator == nullptr ) {
		random_generator = RandomGeneratorOP( new RandomGenerator );
	}
	return *random_generator;
}

double uniform() { return rg().uniform(); }
double gaussian() { return rg().gaussian(); }
int random_range(int low, int high) { return rg().random_range(low, high); }

/// @brief Generate a random deviate from the Von Mises distribution centered at 0
/// with concentration parameter kappa, using the Best-Fisher (1979) algorithm.
/// @details When kappa is 0, returns a uniform random angle in (-pi, pi].
/// For large kappa, the distribution approximates a Gaussian with sigma = 1/sqrt(kappa).
/// @author Andy Watkins
double von_mises( double kappa ) {
	double const pi = numeric::constants::d::pi;

	// Special case: kappa == 0 gives uniform distribution on the circle
	if ( kappa < 1e-12 ) {
		return ( uniform() * 2.0 - 1.0 ) * pi;
	}

	// Best-Fisher algorithm (Best & Fisher, 1979, Applied Statistics 28:152-157)
	double const tau = 1.0 + std::sqrt( 1.0 + 4.0 * kappa * kappa );
	double const rho = ( tau - std::sqrt( 2.0 * tau ) ) / ( 2.0 * kappa );
	double const r = ( 1.0 + rho * rho ) / ( 2.0 * rho );

	while ( true ) {
		double const u1 = uniform();
		double const z = std::cos( pi * u1 );
		double const f = ( 1.0 + r * z ) / ( r + z );
		double const c = kappa * ( r - f );

		double const u2 = uniform();

		if ( c * ( 2.0 - c ) > u2 || std::log( c / u2 ) + 1.0 >= c ) {
			double const u3 = uniform();
			double const theta = ( u3 > 0.5 ) ? std::acos( f ) : -std::acos( f );
			return theta;
		}
	}
}


using namespace std;

/// uniform_RG factory
uniform_RG_OP
createRG(string const & type )
{
	if ( type == "standard" ) return utility::pointer::make_shared< standard_RG >();
	//if( type == "ran3" ) return new ran3_RG();
	if ( type == "mt19937" ) return utility::pointer::make_shared< mt19937_RG >();

	utility_exit_with_message("Unknown random number generator type: " + type);
	return nullptr;
}

RandomGenerator::RandomGenerator() :
	gaussian_iset_( true ),
	gaussian_gset_( 0 )
{}

RandomGenerator::~RandomGenerator() = default;

/// @details The RNG has not been intiliaed if the generator_ object points at null.
bool RandomGenerator::initialized() const { return generator_ != nullptr; }

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
	double rgaussian; // Return value
	if ( gaussian_iset_ ) {
		double v1, v2, rsq, fac;
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
	// std::cout << "NR Gaussian: " << gaussian << std::endl;
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

int RandomGenerator::random_range2(int low, int high){
	if ( low == high ) {
		return low;
	} else {
		return random_range( low, high);
	}

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
