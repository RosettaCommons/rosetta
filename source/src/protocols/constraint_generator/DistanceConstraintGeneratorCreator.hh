// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/DistanceConstraintGeneratorCreator.hh
/// @brief Generates AtomPair constraints to enforce a given distance between two residue subsets
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_constraint_generator_DistanceConstraintGeneratorCreator_hh
#define INCLUDED_protocols_constraint_generator_DistanceConstraintGeneratorCreator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>

namespace protocols {
namespace constraint_generator {

class DistanceConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	virtual protocols::constraint_generator::ConstraintGeneratorOP
	create_constraint_generator() const;

	virtual std::string
	keyname() const;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols/constraint_generator_DistanceConstraintGenerator_fwd_hh
