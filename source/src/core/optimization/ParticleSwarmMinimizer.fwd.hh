// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/ParticleSwarmMinimizer.fwd.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_core_optimization_ParticleSwarmMinimizer_fwd_hh
#define INCLUDED_core_optimization_ParticleSwarmMinimizer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.fwd.hh>


namespace core {
namespace optimization {


class Particle; // fwd declaration
typedef utility::pointer::shared_ptr< Particle > ParticleOP;
typedef utility::pointer::shared_ptr< Particle const > ParticleCOP;
typedef utility::vector1< ParticleOP > ParticleOPs;

class ParticleSwarmMinimizer; // fwd declaration
typedef utility::pointer::shared_ptr< ParticleSwarmMinimizer > ParticleSwarmMinimizerOP;
typedef utility::pointer::shared_ptr< ParticleSwarmMinimizer const > ParticleSwarmMinimizerCOP;


} // namespace optimization
} // namespace core

#endif // INCLUDED_core_optimization_ParticleSwarmMinimizer_FWD_HH
