// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FitnessFunction.hh
/// @brief determines the fitness of Entity for GeneticAlgorithm
/// @author ashworth

#ifndef INCLUDED_protocols_genetic_algorithm_FitnessFunction_hh
#define INCLUDED_protocols_genetic_algorithm_FitnessFunction_hh

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/genetic_algorithm/Entity.fwd.hh>

#include <core/types.hh>

namespace protocols {
namespace genetic_algorithm {

class FitnessFunction : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::shared_ptr< FitnessFunction > OP;
	typedef utility::pointer::shared_ptr< FitnessFunction const > COP;
	~FitnessFunction() override= default;
	virtual core::Real evaluate( Entity & entity ) = 0;
};

} // namespace genetic_algorithm
} // namespace protocols

#endif
