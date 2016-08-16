// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/GA_Minimizer.cc
/// @brief  Minimizer based on Genetic Algorithm
/// @author Sergey Lyskov


#include <core/optimization/GA_Minimizer.hh>

#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <algorithm>

#include <core/optimization/Multifunc.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {

static THREAD_LOCAL basic::Tracer TR( "core.optimization.GA_Minimizer" );


using core::Size;

/// Return true with given probability
static bool yes_no_random(Real probability)
{
	if ( numeric::random::rg().uniform() < probability ) return true;
	else return false;
}


// starting position, and solution is returned here
Real GA_Minimizer::run( Multivec & v, int max_time)
{
	allowed_time_ = max_time;

	best_.v.resize( v.size() );
	for ( Size i=1; i<=v.size(); i++ ) {
		best_.v[i] = v[i];
	}
	best_.r = func_(best_.v);

	TR << "Initial score (time=0): " << best_.r << /*" V=" << best_.v << */ std::endl;


	int ctime = allowed_time_;

	TR << "Randomizing..." << std::endl;
	EItem r = randomize(best_, ctime);
	for ( Size i=1; i<=v.size(); i++ ) {
		v[i] = r.v[i];
	}

	return func_(v);
}

EItem GA_Minimizer::randomize(const EItem& sit, int &time)
{
	std::vector<EItem> pop;
	for ( int i=0; i<20; i++ ) {
		EItem t = sit;  t.tag = 'r';
		for ( Size j=1; j<=t.v.size(); j++ ) {
			// Original: Random::DoubleUniformRandom()*2. - 1.;
			t.v[j] = sit.v[j] + numeric::random::rg().uniform()*2. - 1.;
		}

		t.r = func_(t.v);

		pop.push_back(t);
		time--;
	}
	EItem t = sit;  t.tag = 'O';  t.r = func_(t.v);
	pop.push_back(t);

	best_ = pop[0];

	loop(pop, time);  return best_;
}

EItem GA_Minimizer::loop(std::vector<EItem> & pop, int &time)
{
	EItem shift;  for ( Size i=1; i<=best_.v.size(); i++ ) shift.v.push_back(0);
	int mres = 10;

	for ( ; time>0; ) {
		if ( mres < -5 ) {
			TR << "Shifting...(" << allowed_time_-time<< ") -> " << pop[0].r << std::endl;
			std::vector<EItem> npop;

			if ( add_original_ ) {
				EItem t = pop[0];  fill(t.v.begin(), t.v.end(), 1.e-10);
				npop.push_back( t );
			}

			for ( int i=0; i<10; i++ ) {
				EItem t = pop[0];  t.tag = 'r';
				for ( Size j=1; j<=t.v.size(); j++ ) {
					t.v[j] = numeric::random::rg().uniform()*2. - 1.; //Random::DoubleUniformRandom()*2. - 1.;
				}

				t.r = func_(t.v);

				npop.push_back(t);
			}

			for ( Size i=1; i<=shift.v.size(); i++ ) shift.v[i] += pop[0].v[i];

			mres = 7;

			pop = npop;
			time -= pop.size();
		} else {
			step(pop, time, mres, shift);
			if ( pop[0].r < min_error_ ) {
				TR << " Best < MaxError... - stoping, time=" << time << "." << std::endl;
				for ( Size i=1; i<=shift.v.size(); i++ ) pop[0].v[i] += shift.v[i];
				//StopTime = time;
				return pop[0];
			}
		}
	}

	//StopTime = time;
	for ( Size i=1; i<=shift.v.size(); i++ ) pop[0].v[i] += shift.v[i];
	return pop[0];
}

void GA_Minimizer::step(std::vector<EItem> &pop, int &c_time, int &mres, EItem &shift)
{
	int i_ev = pop.size();  double best = pop[0].r;

	for ( int i=pop.size()-1; i>=0; i-- ) {  // for each item in pop do: mutation/crossover/...
		// Plain Mutation
		if ( yes_no_random(1.) ) {
			c_time--;
			int last_i = pop.size();
			pop.push_back(pop[i]);
			mutate(pop[last_i]);
		}

		// Cross over mutation
		if ( yes_no_random(1.) ) {
			c_time--;
			int last_i = pop.size();
			pop.push_back(pop[i]);

			int i2 = numeric::random::rg().random_range(0, i_ev-1);  //Random::RangeRandom(i_ev);
			if ( i2 != i ) cross_over(pop[last_i], pop[i], pop[i2]);
		}
	}

	EItem tmp;  // evaluate new items
	tmp.v.resize( shift.v.size() );
	for ( Size i=i_ev; i<pop.size(); i++ ) {
		for ( Size j=1; j<=pop[i].v.size(); j++ ) {
			tmp.v[j] = pop[i].v[j] + shift.v[j];
		}

		pop[i].r = func_(tmp.v);
	}

	// sort pop by perfomarce, and take best 10 items.
	sort(pop.begin(), pop.end(), EItem::sort_R_function);
	if ( pop.size() > 10 ) pop.erase(pop.begin()+10, pop.end());

	// see if best value improved
	if ( best > pop[0].r + minimize_tolerance_ ) {
		TR << pop[0].tag;
		mres++;
	} else {
		mres--;
	}

	if ( best_.r > pop[0].r ) {        // global best value;
		best_ = pop[0];
		for ( Size j=1; j<=best_.v.size(); j++ ) {
			best_.v[j] = best_.v[j] + shift.v[j];
		}

		TR << " Time=" << allowed_time_ - c_time << " Score=" << best_.r /*<< " V=" << best_.v */<< std::endl;
	}
}


void GA_Minimizer::mutate(EItem &V)
{
	for ( Size i=1; i<=V.v.size(); i++ ) {
		if ( yes_no_random( mutation_probability_ ) ) {
			double r = numeric::random::rg().gaussian() + .7; //Random::NormalRandom(1, .7);
			V.v[i] *= r * 1.0;
			//V.v[i] += numeric::random::rg().gaussian();
		}
	}
	V.tag='m';
}

void GA_Minimizer::cross_over(EItem &V, EItem &A, EItem &B)
{
	for ( Size i=1; i<=V.v.size(); i++ ) {
		if ( yes_no_random(.5) ) V.v[i] = A.v[i];
		else V.v[i] = B.v[i];
	}
	mutate(V) ;
	V.tag='c';
}

} // namespace optimization
} // namespace core
