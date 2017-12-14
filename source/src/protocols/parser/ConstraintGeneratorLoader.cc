// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/parser/ConstraintGeneratorLoader.hh>
#include <protocols/parser/ConstraintGeneratorLoaderCreator.hh>

// Project headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.ConstraintGeneratorLoader" );

ConstraintGeneratorLoader::ConstraintGeneratorLoader() = default;
ConstraintGeneratorLoader::~ConstraintGeneratorLoader() = default;

void ConstraintGeneratorLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	utility::vector0< utility::tag::TagCOP > const & selector_tags( tag->getTags() );
	for ( utility::tag::TagCOP ii_tag: selector_tags ) {
		constraint_generator::ConstraintGeneratorOP cst_gen = constraint_generator::ConstraintGeneratorFactory::get_instance()->new_constraint_generator(
			ii_tag->getName(),
			ii_tag,
			datamap
		);

		// if "name" is specified, add it to the data map under that name. Otherwise use the type name.
		bool const data_add_status = datamap.add( "ConstraintGenerators" ,
			ii_tag->getOption( "name", ii_tag->getName() ),
			cst_gen );
		if ( !data_add_status ) {
			utility_exit_with_message( "ConstraintGenerator '" + ii_tag->getName() + "' already exists in the basic::datacache::DataMap. Please rename." );
		}
	}
	TR.flush();
}

std::string
ConstraintGeneratorLoader::loader_name() { return "CONSTRAINT_GENERATORS"; }

std::string
ConstraintGeneratorLoader::cst_gen_loader_ct_namer( std::string const & element_name )
{
	return "constraint_generator_loader_" + element_name + "_type";
}

void ConstraintGeneratorLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	constraint_generator::ConstraintGeneratorFactory::get_instance()->define_constraint_generator_xml_schema_group( xsd );

	XMLSchemaSimpleSubelementList cst_gen_loader_subelements;
	cst_gen_loader_subelements.add_group_subelement( & constraint_generator::ConstraintGeneratorFactory::constraint_generator_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator cst_gen_ct;
	cst_gen_ct.element_name( loader_name() ).complex_type_naming_func( & cst_gen_loader_ct_namer )
		.description( "ConstraintGenerators may be defined as subelements of the " + loader_name() + " element, and then will be placed into the DataMap"
		" for later retrieval by Movers and Filters or anything else that might use a ConstraintGenerator. All immediate subelements should have the 'name' attribute"
		" as that is how they will be identified in the DataMap. NOTE: The BluePrintBDR only takes a SheetRemodelConstraintGenerator, which is actually"
		" a mover and is declared in the MOVERS block as \"SheetCstGenerator\". Passing an actualy constraint generator to this mover will fail."
		" Constraint generators are currently only used by the AddConstraints mover."
		)
		.set_subelements_repeatable( cst_gen_loader_subelements )
		.write_complex_type_to_schema( xsd );

}

DataLoaderOP
ConstraintGeneratorLoaderCreator::create_loader() const { return DataLoaderOP( new ConstraintGeneratorLoader ); }

std::string
ConstraintGeneratorLoaderCreator::keyname() const { return ConstraintGeneratorLoader::loader_name(); }

ConstraintGeneratorLoaderCreator::DerivedNameFunction
ConstraintGeneratorLoaderCreator::schema_ct_naming_function() const
{
	return & ConstraintGeneratorLoader::cst_gen_loader_ct_namer;
}

void
ConstraintGeneratorLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConstraintGeneratorLoader::provide_xml_schema( xsd );
}



} //namespace parser
} //namespace protocols
