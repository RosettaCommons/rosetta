// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/ParticleSwarmMinimizer.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_core_optimization_ParticleSwarmMinimizer_hh
#define INCLUDED_core_optimization_ParticleSwarmMinimizer_hh

#include <core/optimization/ParticleSwarmMinimizer.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/optimization/Multifunc.hh>
#include <core/optimization/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {


/// @brief Simple data container for PSO algorithm.
class Particle : public utility::pointer::ReferenceCount
{
public:

	Particle(Size size):
		p_(size, 0.0),
		fitness_(0.0),
		v_(size, 0.0),
		best_valid_( false ),
		pbest_(size, 0.0),
		fitness_pbest_(0.0)
	{}

	Particle(Multivec const & p_in):
		p_(p_in),
		fitness_(0.0),
		v_(p_in.size(), 0.0),
		best_valid_( false ),
		pbest_(p_in),
		fitness_pbest_(0.0)
	{}

	virtual ~Particle() {}

	Real score(Multifunc & f)
	{
		// Reverse the sign: for normal multifunc, lower is better
		// For historical reasons, the code maximizes the "fitness"
		fitness_ = -f(p_);
		if ( ! best_valid_ || fitness_pbest_ < fitness_ ) {
			best_valid_ = true;
			pbest_ = p_; // make a copy
			fitness_pbest_ = fitness_;
		}
		return fitness_;
	}

	Real set_score(Real & new_score)
	{
		// Reverse the sign: for normal multifunc, lower is better
		// For historical reasons, the code maximizes the "fitness"
		fitness_ = -new_score;
		if ( ! best_valid_ || fitness_pbest_ < fitness_ ) {
			best_valid_ = true;
			pbest_ = p_; // make a copy
			fitness_pbest_ = fitness_;
		}
		return fitness_;
	}

	/// @brief Make sure that all arrays are large enough -- prevents index-out-of-bound errors.
	void ensure_size(Size minsize)
	{
		if ( p_.size() < minsize ) p_.resize(minsize);
		Size const s = p_.size();
		if ( v_.size() < s ) v_.resize(s);
		if ( pbest_.size() < s ) pbest_.resize(s);
	}

	/// @brief This is why data should be private: you get to ensure it's valid when you read it.
	Multivec const &
	pbest() const {
		debug_assert( best_valid_ );
		return pbest_;
	}

	Real
	fitness_pbest() const {
		debug_assert( best_valid_ );
		return fitness_pbest_;
	}

public:
	Multivec p_;         //< current position (list of reals)
	Real fitness_;       //< current fitness, higher is better (real)
	Multivec v_;         //< current velocity (list of reals)

private:
	bool best_valid_;    //< Does the fitness_pbest_ value represent a real value of the multifunc?
	Multivec pbest_;     //< copy of the "p" with the highest fitness seen by this particle
	Real fitness_pbest_; //< fitness value for pbest

}; // Particle

std::ostream & operator << ( std::ostream & os, Particle const & p );


/// @brief Particle Swarm Optimization engine.
///
/// @details Algorithm details based heavily on
///
///     Chen, Liu, Huang, Hwang, Ho (2006).
///     "SODOCK:  Swarm Optimization for Highly Flexible Protein-Ligand Docking"
///     J Comput Chem 28: 612-623, 2007
///
/// Also on
///     http://en.wikipedia.org/wiki/Particle_swarm_optimization
///     http://www.swarmintelligence.org/
///
/// One can imagine writing another version that distributed the work via MPI...
class ParticleSwarmMinimizer : public utility::pointer::ReferenceCount
{
public:

	ParticleSwarmMinimizer(Multivec p_min, Multivec p_max);
	virtual ~ParticleSwarmMinimizer();

	ParticleOPs run(Size num_cycles, Multifunc & f_fitness, Size num_part = 50);
	ParticleOPs run(Size num_cycles, Multifunc & f_fitness, Size num_part, Multivec init_values );
	void run(Size num_cycles, Multifunc & f_fitness, ParticleOPs & particles);
	void print_particles( ParticleOPs & particles, std::string header );

protected:

	virtual void score_all_particles(Multifunc & f_fitness, ParticleOPs & particles);

private:

	Size size_; //< number of degrees of freedom
	Real C_inertia_start_;
	Real C_inertia_end_;
	Real C_pbest_; //< the "cognitive" parameter
	Real C_lbest_; //< the neighborhood "social" parameter
	Real C_gbest_; //< the global "social" parameter
	int first_nbr_;
	int last_nbr_;
	Multivec p_min_;
	Multivec p_max_;
	Multivec p_range_;
	Multivec v_max_;

}; // ParticleSwarmMinimizer


} // namespace optimization
} // namespace core

#endif // INCLUDED_core_optimization_ParticleSwarmMinimizer_HH
