// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/SheetRemodelConstraintGenerator.cc
/// @brief Remodel constraint generator for adding sheet constraints
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit header
#include <protocols/fldsgn/SheetRemodelConstraintGenerator.hh>
#include <protocols/fldsgn/SheetRemodelConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/fldsgn/SheetConstraintGenerator.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fldsgn.SheetRemodelConstraintGenerator" );

namespace protocols {
namespace fldsgn {

SheetRemodelConstraintGenerator::SheetRemodelConstraintGenerator():
	protocols::forge::remodel::RemodelConstraintGenerator(),
	cg_( new SheetConstraintGenerator )
{
}

SheetRemodelConstraintGenerator::~SheetRemodelConstraintGenerator() = default;

// XRW TEMP std::string
// XRW TEMP SheetRemodelConstraintGenerator::get_name() const
// XRW TEMP {
// XRW TEMP  return SheetRemodelConstraintGenerator::mover_name();
// XRW TEMP }

protocols::moves::MoverOP
SheetRemodelConstraintGenerator::clone() const
{
	return protocols::moves::MoverOP( new SheetRemodelConstraintGenerator( *this ) );
}

void
SheetRemodelConstraintGenerator::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	cg_->parse_my_tag( tag, data );
}

void
SheetRemodelConstraintGenerator::generate_remodel_constraints( core::pose::Pose const & pose )
{
	if ( !cg_ ) {
		std::stringstream msg;
		msg << mover_name() << "::generate_remodel_constraints(): "
			<< "constraint generator not set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	add_constraints( cg_->apply( pose ) );
}

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SheetRemodelConstraintGeneratorCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new SheetRemodelConstraintGenerator );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SheetRemodelConstraintGeneratorCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SheetRemodelConstraintGenerator::mover_name();
// XRW TEMP }

std::string SheetRemodelConstraintGenerator::get_name() const {
	return mover_name();
}

std::string SheetRemodelConstraintGenerator::mover_name() {
	return "SheetCstGenerator";
}

void SheetRemodelConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	SheetConstraintGenerator::attributes_for_sheet_constraint_generator( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SheetRemodelConstraintGeneratorCreator::keyname() const {
	return SheetRemodelConstraintGenerator::mover_name();
}

protocols::moves::MoverOP
SheetRemodelConstraintGeneratorCreator::create_mover() const {
	return protocols::moves::MoverOP( new SheetRemodelConstraintGenerator );
}

void SheetRemodelConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SheetRemodelConstraintGenerator::provide_xml_schema( xsd );
}


} //protocols
} //fldsgn

