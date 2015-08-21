// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/md/thermostat.hh
/// @brief  Thermostat for MD simulation
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_md_thermostat_hh
#define INCLUDED_protocols_md_thermostat_hh

#include <protocols/md/MDConstants.hh>
#include <core/types.hh>
#include <core/optimization/types.hh>

#include <cmath>
#include <iostream>

namespace protocols {
namespace md {

using namespace core;
using namespace core::optimization;
using namespace protocols::md;

class Thermostat {

public:

	Thermostat( Real const &temperature0, Size const ndof ){
		temp0_ = temperature0*Boltzmann;
		tau_ = 0.015;
		nstep_per_update_ = 10;
		ndof_ = ndof;
	}

	Thermostat( Real const &temperature0, Real const &tau, Size const ndof ){
		temp0_ = temperature0*Boltzmann;
		tau_ = tau;
		nstep_per_update_ = 10;
		ndof_ = ndof;
	}

	~Thermostat(){}

	void
	rescale( Multivec &vel, Real const& dt, Multivec const &mass ){
		Real curr_temperature = get_temperature( vel, mass );

		Real delta  = (temp0_/Boltzmann)/curr_temperature;
		Real lambda = std::sqrt(1.0 + dt/tau_*(delta-1.0));

		for ( Size i_dof = 1; i_dof <= vel.size(); ++i_dof ) {
			vel[i_dof] *= lambda;
		}
	}

	Real
	get_temperature( Multivec const &vel, Multivec const &mass)
	{
		Real temperature = 0.0;

		for ( Size i_dof = 1; i_dof<=vel.size(); ++i_dof ) {
			int i_atm = (i_dof+2)/3;
			// pass Virtual atoms
			if ( mass[i_atm] < 1e-3 ) continue;

			temperature += mass[i_atm]*vel[i_dof]*vel[i_dof];
		}
		temperature /= ndof_*Boltzmann;
		return temperature;
	}

	Size
	nstep_per_update(){ return nstep_per_update_; }

private:
	Real temp0_;
	Real tau_;
	Size nstep_per_update_;
	Size ndof_;
};
}
}

#endif
