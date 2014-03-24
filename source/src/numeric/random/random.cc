// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/random/random.cc
/// @brief  Random number generator system
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)
///
/// @remarks
///  @li -
///


// Unit headers
#include <numeric/random/random.hh>

// Package headers
#include <numeric/random/uniform.hh>
//#include <numeric/random/ran3.hh>
#include <numeric/random/mt19937.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <iostream>
#include <cmath>


namespace numeric {
namespace random {

RandomGenerator RG(0);

double uniform(void) { return RG.uniform(); }
double gaussian(void) { return RG.gaussian(); }
int random_range(int low, int high) { return RG.random_range(low, high); }


using namespace std;

vector<RandomGenerator*> &RandomGenerator::allGenerators()
{
	static vector<RandomGenerator*> * allGen = new vector<RandomGenerator*>;
	return *allGen;
}


/// uniform_RG factory
uniform_RG * createRG(string const & type )
{
	if( type == "standard" ) return new standard_RG();
	//if( type == "ran3" ) return new ran3_RG();
	if( type == "mt19937" ) return new mt19937_RG();

	utility_exit_with_message("Unknown random number generator type: " + type);
	return 0;
}


void RandomGenerator::initializeRandomGenerators(int const _start_seed, RND_RunType run_type,
												 string const & RGtype)
{
	uint32_t start_seed(_start_seed);

	if( run_type == _RND_TestRun_ ) {
		for(Size i=0; i<allGenerators().size(); i++) {
			uint32_t offset(allGenerators()[i]->seed_offset);
			offset += start_seed;

			allGenerators()[i]->generator = createRG(RGtype);

			//allGenerators()[i]->generator->setSeed(allGenerators()[i]->seed_offset+start_seed);
			///^ ToDo: fix possible int overflow here

			allGenerators()[i]->generator->setSeed(offset);
		}
	}
	else { // Normal run, we use only one generator here then
		utility::pointer::owning_ptr<uniform_RG> RG = createRG(RGtype);
		for(Size i=0; i<allGenerators().size(); i++) allGenerators()[i]->generator = RG;
		RG->setSeed(start_seed);
	}
}

void RandomGenerator::saveState(std::ostream & out)
{
	out << " " << seed_offset;
	out << " " << gaussian_iset;
	// Cast raw bits to int for perfect precision:
	out << " " << *((uint64_t*) &gaussian_gset);
	generator->saveState(out);
	out << "\n";
}

void RandomGenerator::restoreState(std::istream & in)
{
	in >> seed_offset;
	in >> gaussian_iset;
	// Cast raw bits from int for perfect precision:
	uint64_t gset_bits = 0;
	in >> gset_bits;
	gaussian_gset = *((double*) &gset_bits);
	generator->restoreState(in);
}

void RandomGenerator::saveAllStates(std::ostream & out)
{
	for(Size i=0; i<allGenerators().size(); i++) {
		allGenerators()[i]->saveState(out);
	}
}

void RandomGenerator::restoreAllStates(std::istream & in)
{
	for(Size i=0; i<allGenerators().size(); i++) {
		allGenerators()[i]->restoreState(in);
	}
}

RandomGenerator::RandomGenerator(int const magicNumber) :
	seed_offset(magicNumber),
	gaussian_iset( true ),
	gaussian_gset( 0 )
{
	for ( Size i = 0; i < allGenerators().size(); i++ ) { //< checking for magic number dublicate
		if ( allGenerators()[i]->seed_offset == magicNumber ) {
			cerr << "Duplicate magic number in random object initialization! Number:"
				 << magicNumber << endl;
			utility_exit();
		}
	}
	allGenerators().push_back(this);
}

RandomGenerator::~RandomGenerator()
{
	for(Size i = 0; i < allGenerators().size(); i++) {
		if( allGenerators()[i] == this ) {
			allGenerators().erase(allGenerators().begin() + i);
			break;
		}
	}
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
	if ( gaussian_iset ) {
		do {
			v1 = 2.0f * uniform() - 1.0f;
			v2 = 2.0f * uniform() - 1.0f;
			rsq = ( v1 * v1 ) + ( v2 * v2 );
		} while ( rsq >= 1.0 || rsq == 0.0 );
		fac = std::sqrt(-(2.0*std::log(rsq)/rsq));
		gaussian_gset = v1*fac;
		rgaussian = v2*fac;
		gaussian_iset = false;
	} else {
		rgaussian = gaussian_gset;
		gaussian_iset = true;
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


} // namespace random
} // namespace numeric
