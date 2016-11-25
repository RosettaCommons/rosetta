// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/TerminiConstraintGenerator.hh
/// @brief Generates distance constraints between the upper and lower termini
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_TerminiConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_TerminiConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/TerminiConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace constraint_generator {

///@brief Generates distance constraints between the upper and lower termini
class TerminiConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	TerminiConstraintGenerator();

	~TerminiConstraintGenerator() override;

	static std::string
	class_name() { return "TerminiConstraintGenerator"; }

	ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	set_min_distance( core::Real const dist );

	void
	set_max_distance( core::Real const dist );

	void
	set_sd( core::Real const sd );

	void
	set_weight( core::Real const weight );

	core::Real
	weight() const;

	core::Real
	min_distance() const;

	core::Real
	max_distance() const;

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

private:
	core::Real min_distance_;
	core::Real max_distance_;
	core::Real sd_;
	core::Real weight_;

};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_TerminiConstraintGenerator_hh

