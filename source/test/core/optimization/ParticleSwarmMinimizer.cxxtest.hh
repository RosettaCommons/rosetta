// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/io/izstream.cxxtest.hh
/// @brief  zipstream unit test suite
/// @author Ian Davis

// Package headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.optimization.ParticleSwarmMinimizer.cxxtest");

// Accessory class
class SimpleMultifunc : public core::optimization::Multifunc {
public:

	/// @brief Stupid object function:  all DOFs should sum to 1.0
	virtual
	core::Real
	operator ()( core::optimization::Multivec const & phipsi ) const
	{
		core::Real sum = 0.0;
		for(core::Size i = 1; i <= phipsi.size(); ++i) sum += phipsi[i];
		return std::abs( 1.0 - sum );
	}

	virtual
	void
	dfunc( core::optimization::Multivec const & /*phipsi*/, core::optimization::Multivec & /*dE_dphipsi*/ ) const
	{
		utility_exit_with_message("Doesn't support dfunc!");
	}
};


class ParticleSwarmMinimizerTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	/// @brief Test particle swarm optimizer with a very simple fitness function
	void test_simple_multifunc() {
		using namespace core::optimization;
		SimpleMultifunc f;
		core::Size const num_dofs = 5;
		Multivec p_min(num_dofs, 0.), p_max(num_dofs, 1.);
		ParticleSwarmMinimizer pso(p_min, p_max);
		ParticleOPs particles = pso.run(5000, f);
		for(core::Size i = 1; i <= particles.size(); ++i) {
			ParticleOP p = particles[i];
			// Remember that fitness = -f() ...
			TR << "Particle " << i << ": fitness=" << p->fitness_pbest() << " [";
			for(core::Size j = 1; j <= p->pbest().size(); ++j) {
				TR << " " << p->pbest()[j];
			}
			TR << " ]" << std::endl;
		}
		TS_ASSERT( f(particles[1]->p_) < 1e-6 );
	}

};

