// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/CoordinateConstraintGeneratorCreator.hh
/// @brief Generates coodinate constraints
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_CoordinateConstraintGeneratorCreator_hh
#define INCLUDED_protocols_constraint_generator_CoordinateConstraintGeneratorCreator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>

namespace protocols {
namespace constraint_generator {

class CoordinateConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	virtual protocols::constraint_generator::ConstraintGeneratorOP
	create_constraint_generator() const;

	virtual std::string
	keyname() const;

};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols/constraint_generator_CoordinateConstraintGenerator_fwd_hh
