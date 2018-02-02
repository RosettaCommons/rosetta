// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/SegmentedAtomPairConstraintGeneratorCreator.hh
/// @brief Given a set of non-continuous selected segments, generates differently scored atom pair constraints
///       for the resides in each segment and between segments.
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#ifndef INCLUDED_protocols_fold_from_loops_constraint_generator_SegmentedAtomPairConstraintGeneratorCreator_hh
#define INCLUDED_protocols_fold_from_loops_constraint_generator_SegmentedAtomPairConstraintGeneratorCreator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace fold_from_loops {
namespace constraint_generator {

class SegmentedAtomPairConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	protocols::constraint_generator::ConstraintGeneratorOP
	create_constraint_generator() const override;

	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
} //protocols
} //fold_from_loops

#endif //INCLUDED_protocols_fold_from_loops_SegmentedAtomPairConstraintGeneratorCreator_hh
