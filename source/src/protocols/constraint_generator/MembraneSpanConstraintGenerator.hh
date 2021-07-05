// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MembraneSpanConstraintGenerator.hh
/// @brief Generates constraints for membrane spans to stay in the membrane
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_constraint_generator_MembraneSpanConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_MembraneSpanConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/MembraneSpanConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace constraint_generator {

///@brief Generates constraints for membrane spans to stay in the membrane
class MembraneSpanConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	MembraneSpanConstraintGenerator();

	~MembraneSpanConstraintGenerator() override;

	static std::string
	class_name();

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

	/// @brief Provide citations to the passed CitationCollectionList.
	/// This allows the constraint generator to provide citations for itself
	/// and for any modules that it invokes.
	/// @details Cites Jonathan Weinstein.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

protected:
	void
	parse_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;

public:
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_MembraneSpanConstraintGenerator_hh

