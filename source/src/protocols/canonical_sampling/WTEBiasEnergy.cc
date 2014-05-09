// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/WTEBiasEnergyMover.cc
/// @brief WTEBiasEnergy methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/WTEBiasEnergy.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

//static basic::Tracer tr( "protocols.canonical_sampling.BiasEnergy" );
//static numeric::random::RandomGenerator RG(2592747);

namespace protocols {
namespace canonical_sampling {

using namespace core;

WTEBiasEnergy::WTEBiasEnergy( core::Size stride, core::Real omega, core::Real gamma ) :
	Parent( stride,omega,gamma)
{}

WTEBiasEnergy::WTEBiasEnergy() {}


Real WTEBiasEnergy::extract_collective_var( core::pose::Pose const& pose ) const {
	return pose.energies().total_energy();
}

} //moves
} //protocols

