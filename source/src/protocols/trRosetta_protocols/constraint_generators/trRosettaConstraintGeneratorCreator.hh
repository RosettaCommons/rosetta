// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGeneratorCreator.hh
/// @brief A module that runs a trRosetta neural network on an input multiple sequence alignment and
/// uses the output to apply distance and/or angle constraints to a pose for subsequent structure
/// prediction or refinement.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_protocols_constraint_generators_trRosettaConstraintGeneratorCreator_hh
#define INCLUDED_protocols_trRosetta_protocols_constraint_generators_trRosettaConstraintGeneratorCreator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>

namespace protocols {
namespace trRosetta_protocols {
namespace constraint_generators {

class trRosettaConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	protocols::constraint_generator::ConstraintGeneratorOP
	create_constraint_generator() const override;

	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //constraint_generators
} //trRosetta_protocols
} //protocols

#endif //INCLUDED_protocols/trRosetta_protocols/constraint_generators_trRosettaConstraintGenerator_fwd_hh
