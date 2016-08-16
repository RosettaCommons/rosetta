// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/ParticleSwarmMinimizer.cc
///
/// @brief
/// @author Ian W. Davis


#include <core/optimization/ParticleSwarmMinimizer.hh>

#include <utility/exit.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>

#include <algorithm>

#include <utility/vector1.hh>


namespace core {
namespace optimization {

using namespace ObjexxFCL::format;

/// @brief stream output operator for Particle types
std::ostream &
operator<< ( std::ostream & os, Particle const & p ) {
	os << " best fitness: " << ObjexxFCL::format::F( 9,5,-1.0 * p.fitness_pbest() )
		<< ", current fitness: " << ObjexxFCL::format::F( 9,6,-1.0 * p.fitness_ ) << ", current dofs: [";
	for ( core::Size i=1; i <= p.p_.size(); ++i ) { os << ObjexxFCL::format::F( 8,4,p.p_[i] ) << ", "; }
	os << " ]";
	return os;
}


// Used to sort particles from best to worst
bool cmp_particles(ParticleOP a, ParticleOP b)
{
	return a->fitness_pbest() > b->fitness_pbest();
}


ParticleSwarmMinimizer::ParticleSwarmMinimizer(Multivec p_min, Multivec p_max):
	utility::pointer::ReferenceCount(),
	size_(p_min.size()),
	C_inertia_start_(0.9),
	C_inertia_end_(0.4),
	C_pbest_(2.0),
	C_lbest_(2.0),
	C_gbest_(0.0),
	first_nbr_(-2),
	last_nbr_(2),
	p_min_(p_min),
	p_max_(p_max),
	p_range_(),
	v_max_()
{
	runtime_assert(p_min_.size() == p_max_.size());
	p_range_.resize(size_, 0.0);
	v_max_.resize(size_, 0.0);
	for ( Size i = 1; i <= size_; ++i ) {
		runtime_assert( p_min_[i] < p_max_[i] );
		p_range_[i] = p_max_[i] - p_min_[i];
		v_max_[i] = 0.1 * p_range_[i];
	}
}


ParticleSwarmMinimizer::~ParticleSwarmMinimizer() {}


ParticleOPs ParticleSwarmMinimizer::run(Size num_cycles, Multifunc & f_fitness, Size num_part /*= 50*/)
{
	ParticleOPs particles;
	for ( Size i = 1; i <= num_part; ++i ) {
		ParticleOP p( new Particle(size_) );
		for ( Size j = 1; j <= size_; ++j ) {
			p->p_[j] = p_min_[j] + numeric::random::rg().uniform()*p_range_[j];
		}
		// debugging output
		/*std::cout << "PSM: created new particle: dofs: [ ";
		for ( core::Size k=1; k <= size_; ++k ) {
		std::cout << F(8,4,p->p_[k]) << ", ";
		}
		std::cout << " ]" << std::endl;*/
		particles.push_back(p);
	}
	run(num_cycles, f_fitness, particles);
	return particles;
}


ParticleOPs ParticleSwarmMinimizer::run(Size num_cycles, Multifunc & f_fitness, Size num_part, Multivec init_values )
{
	ParticleOPs particles;
	for ( Size i = 1; i <= num_part; ++i ) {
		ParticleOP p( new Particle(size_) );
		for ( Size j = 1; j <= size_; ++j ) {
			p->p_[j] = init_values[j] + numeric::random::rg().uniform() - numeric::random::rg().uniform(); // want to go up *and* down by a little bit
			// init values should never be outside min/max range
			if ( p->p_[j] < p_min_[j] ) { p->p_[j] = numeric::random::rg().uniform(); }
			else if ( p->p_[j] > p_max_[j] ) { p->p_[j] = p_max_[j]; }
		}
		// debugging output
		/*std::cout << "PSM: created custom init particle: dofs: [ ";
		for ( core::Size k=1; k <= size_; ++k ) {
		std::cout << F(8,4,p->p_[k]) << ", ";
		}
		std::cout << " ]" << std::endl;*/
		particles.push_back(p);
	}
	run(num_cycles, f_fitness, particles);
	return particles;
}


void ParticleSwarmMinimizer::run(Size num_cycles, Multifunc & f_fitness, ParticleOPs & particles)
{
	Size const N = particles.size();
	//runtime_assert( int(N) >= last_nbr_ - first_nbr_ + 1 );
	// Ensure particle vector sizes are consistent and long enough
	for ( Size i = 1; i <= N; ++i ) {
		particles[i]->ensure_size( size_ );
	}
	for ( Size cycle = 1; cycle <= num_cycles; ++cycle ) {
		// linear ramp on inertial weight
		Real const frac_done = Real(cycle) / Real(num_cycles > 1 ? num_cycles-1 : 1);
		Real const C_inertia = (1.0-frac_done)*C_inertia_start_ + frac_done*C_inertia_end_;
		// score everyone and update p(ersonal)best
		score_all_particles(f_fitness,particles);
		//
		// TODO: apply local optimization to best particle in the swarm
		//
		// Determine l(ocal)best and g(lobal)best
		ParticleOPs lbests;
		for ( Size i = 1; i <= N; ++i ) {
			ParticleOP const & p = particles[i];
			lbests.push_back(p); // start by assuming each particle is best among its neighbors
			for ( int jj = int(i)+first_nbr_; jj <= int(i)+last_nbr_; ++jj ) {
				// wrap index around:
				int j = jj;
				if ( j < 1 ) j += N;
				else if ( j > int(N) ) j -= N;
				if ( lbests[i]->fitness_ < particles[j]->fitness_ ) {
					lbests[i] = particles[j];
				}
			}
		}
		ParticleOP gbest = particles[1];
		for ( Size i = 1; i <= N; ++i ) {
			if ( gbest->fitness_ < particles[i]->fitness_ ) {
				// debugging output
				/*std::cout << "PSM: New global best: fitness: " << -1 * particles[i]->fitness_ << ", dofs: [ ";
				for ( core::Size k=1; k <= size_; ++k ) {
				std::cout << F(8,4,particles[i]->p_[k]) << ", ";
				}
				std::cout << " ]" << std::endl;*/
				gbest = particles[i];
			}
		}

		// update velocity
		for ( Size j = 1; j <= N; ++j ) {
			ParticleOP const & p = particles[j];
			ParticleOP const & lbest = lbests[j];
			for ( Size i = 1; i <= size_; ++i ) {
				Real const pi = p->p_[i];

				Real vi = ( C_inertia*p->v_[i]
					+ numeric::random::rg().uniform()*C_pbest_*(p->pbest()[i] - pi)
					+ numeric::random::rg().uniform()*C_lbest_*(lbest->p_[i] - pi)
					+ numeric::random::rg().uniform()*C_gbest_*(gbest->p_[i] - pi) );

				// sometimes particles react too quickly, or move too fast to their local/global best and end up
				// getting stuck at 0.0.
				// slow down how fast the particles move, but not by imposing a speed limit but instead by
				// throttling them when they decide on their new speed.
				/* Real vi = ( C_inertia*p->v_[i]
				+ numeric::random::rg().uniform()*C_pbest_*(p->pbest()[i] - pi)
				+ numeric::random::rg().uniform()*C_lbest_*(lbest->p_[i] - pi)
				+ numeric::random::rg().uniform()*C_gbest_*(gbest->p_[i] - pi) ) * ( ( p_range_[i] ) / num_cycles ); */

				Real const vmax = v_max_[i];
				if ( vi > vmax ) {
					//std::cout << "PSM: particle " << I(2,j) << " DOF " << i << " velocity " << vi << " reset to " << vmax << ". v_: " << p->v_[i]
					// << ", p_: " << p->p_[i] << ", p_best_: " << p->pbest()[i] << ", lbest: " << lbest->p_[i] << std::endl;
					vi = vmax;
				} else if ( vi < -vmax ) {
					//std::cout << "PSM: particle " << I(2,j) << " DOF " << i << " velocity " << vi << " reset to " << -vmax << ". v_: " << p->v_[i]
					// << ", p_: " << p->p_[i] << ", p_best_: " << p->pbest()[i] << ", lbest: " << lbest->p_[i] << std::endl;
					vi = -vmax;
				}
				p->v_[i] = vi;
			}
		}
		// update positions
		for ( Size j = 1; j <= N; ++j ) {
			ParticleOP const & p = particles[j];
			for ( Size i = 1; i <= size_; ++i ) {
				Real ppi = p->p_[i] + p->v_[i];
				if ( ppi < p_min_[i] ) {
					//std::cout << "PSM: particle " << I(2,j) << " DOF " << i << " reset to minimum. p_: " << p->p_[i]
					// << ", v_: " << p->v_[i] << ", p_best_: " << p->pbest()[i] << ", lbest: " << lbests[j]->p_[i] << ", ppi: " << ppi << std::endl;
					ppi = p_min_[i];
				} else if ( ppi > p_max_[i] ) {
					//std::cout << "PSM: particle " << I(2,j) << " DOF " << i << " reset to maximum. p_: " << p->p_[i]
					// << ", v_: " << p->v_[i] << ", p_best_: " << p->pbest()[i] << ", lbest: " << lbests[j]->p_[i] << ", ppi: " << ppi << std::endl;
					ppi = p_max_[i];
				}
				p->p_[i] = ppi;
			}
		}
	}

	// score everyone and update p(ersonal)best one last time
	score_all_particles(f_fitness,particles);
	// sort by fitness, fitest ones first
	std::sort(particles.begin(), particles.end(), cmp_particles);
}


void ParticleSwarmMinimizer::score_all_particles(Multifunc & f_fitness, ParticleOPs & particles) {
	Size const N = particles.size();
	for ( Size i = 1; i <= N; ++i ) {
		particles[i]->score(f_fitness);
	}
}


/// @brief helper function for displaying current particle information; calls the output operator on each particle
void ParticleSwarmMinimizer::print_particles( ParticleOPs & ps, std::string header ) {
	for ( core::Size i=1; i <= ps.size(); ++i ) {
		ParticleOP p = ps[i];
		std::cout << header << ", particle: " << I(2,i);
		std::cout << *p << std::endl;
	}
}


} // namespace optimization
} // namespace core
