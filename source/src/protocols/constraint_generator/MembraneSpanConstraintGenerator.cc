// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MembraneSpanConstraintGenerator.cc
/// @brief Generates constraints for membrane spans to stay in the membrane
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

// Unit headers
#include <protocols/constraint_generator/MembraneSpanConstraintGenerator.hh>
#include <protocols/constraint_generator/MembraneSpanConstraintGeneratorCreator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
// Core headers
#include <core/scoring/constraints/MembraneSpanConstraint.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Boost headers

static basic::Tracer TR( "protocols.constraint_generator.MembraneSpanConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
MembraneSpanConstraintGeneratorCreator::create_constraint_generator() const
{
	return utility::pointer::make_shared< MembraneSpanConstraintGenerator >();
}

std::string
MembraneSpanConstraintGeneratorCreator::keyname() const
{
	return MembraneSpanConstraintGenerator::class_name();
}

MembraneSpanConstraintGenerator::MembraneSpanConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( MembraneSpanConstraintGenerator::class_name() )
{
}

MembraneSpanConstraintGenerator::~MembraneSpanConstraintGenerator() = default;

protocols::constraint_generator::ConstraintGeneratorOP
MembraneSpanConstraintGenerator::clone() const
{
	return utility::pointer::make_shared< MembraneSpanConstraintGenerator >( *this );
}

std::string
MembraneSpanConstraintGenerator::class_name()
{
	return "MembraneSpanConstraintGenerator";
}

void
MembraneSpanConstraintGenerator::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
{
}

core::scoring::constraints::ConstraintCOPs
MembraneSpanConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	//using core::scoring::constraints::ConstraintOP;
	//ConstraintOP amb_cst( new core::scoring::constraints::MembraneSpanConstraint( pose ) );
	//
	core::scoring::constraints::ConstraintCOPs csts;
	csts.push_back( utility::pointer::make_shared< core::scoring::constraints::MembraneSpanConstraint >( pose ) );
	return csts;
}

/// @brief Provide citations to the passed CitationCollectionList.
/// This allows the constraint generator to provide citations for itself
/// and for any modules that it invokes.
/// @details Cites Jonathan Weinstein.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
MembraneSpanConstraintGenerator::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	using namespace basic::citation_manager;
	UnpublishedModuleInfoOP citation(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		class_name(), CitedModuleType::ConstraintGenerator,
		"Jonathan Weinstein", "Department of Biomolecular Science, Weizmann Institute of Science", "jonathan.weinstein@weizmann.ac.il",
		"Created the MembraneSpanConstraintGenerator."
		)
	);
	citations.add(citation);
}


void
MembraneSpanConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	MembraneSpanConstraintGenerator::provide_xml_schema( xsd );
}

void
MembraneSpanConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;

	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"Forms constraints on each span to stay in the membrane",
		attlist );
}


} //protocols
} //constraint_generator

