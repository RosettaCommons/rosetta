// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandDockingLoaders.cc
/// @brief  Implementation of the InterfaceBuilderLoader and MoveMapBuilderLoader classes
/// @author Gordon Lemmon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) -- moved here from DockDesignParser.cc

// Unit Headers
#include <protocols/ligand_docking/LigandDockingLoaders.hh>
#include <protocols/ligand_docking/LigandDockingLoaderCreators.hh>

// Project Headers
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>


#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.LigandDockingLoaders" );

InterfaceBuilderLoader::InterfaceBuilderLoader() {}
InterfaceBuilderLoader::~InterfaceBuilderLoader() = default;

void InterfaceBuilderLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	for ( TagCOP interface_builder_tag : tag->getTags() ) {
		runtime_assert( interface_builder_tag->getName() == InterfaceBuilder::element_name() );
		std::string const name( interface_builder_tag->getOption< std::string >( "name" ));

		if ( data.has("interface_builders", name) ) {
			TR << "WARNING WARNING movemap_builder of name \"" << name
				<< "\" already exists. Skipping\n" << interface_builder_tag << std::endl;
			continue;
		}
		///// Add this movemap to the data map
		protocols::ligand_docking::InterfaceBuilderOP interface_builder( new protocols::ligand_docking::InterfaceBuilder() );
		interface_builder->parse_my_tag( interface_builder_tag, data );
		data.add( "interface_builders" , name, interface_builder);
	}
	TR.flush();
}

std::string InterfaceBuilderLoader::loader_name() { return "INTERFACE_BUILDERS"; }

std::string
InterfaceBuilderLoader::interface_builder_ct_namer( std::string const & element_name )
{
	return "interface_builder_" + element_name + "_type";
}

void InterfaceBuilderLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	InterfaceBuilder::provide_xml_schema( xsd );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( InterfaceBuilder::element_name(), & interface_builder_ct_namer );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.complex_type_naming_func( & interface_builder_ct_namer )
		.set_subelements_repeatable( subelements )
		.description( "Define a set of InterfaceBuilders and load them into the DataMap for later use by MoveMapBuilders" )
		.write_complex_type_to_schema( xsd );
}


parser::DataLoaderOP
InterfaceBuilderLoaderCreator::create_loader() const { return parser::DataLoaderOP( new InterfaceBuilderLoader ); }

std::string
InterfaceBuilderLoaderCreator::keyname() const { return InterfaceBuilderLoader::loader_name(); }

InterfaceBuilderLoaderCreator::DerivedNameFunction
InterfaceBuilderLoaderCreator::schema_ct_naming_function() const
{
	return & InterfaceBuilderLoader::interface_builder_ct_namer;
}

void
InterfaceBuilderLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterfaceBuilderLoader::provide_xml_schema( xsd );
}

MoveMapBuilderLoader::MoveMapBuilderLoader() {}
MoveMapBuilderLoader::~MoveMapBuilderLoader() = default;

void MoveMapBuilderLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	for ( TagCOP movemap_builder_tag : tag->getTags() ) {
		std::string const name( movemap_builder_tag->getOption< std::string> ( "name" ) );

		if ( data.has("movemap_builders", name) ) {
			TR << "WARNING WARNING movemap_builder of name \"" << name
				<< "\" already exists. Skipping\n" << movemap_builder_tag << std::endl;
			continue;
		}
		///// Add this movemap to the data map
		protocols::ligand_docking::MoveMapBuilderOP movemap_builder( new protocols::ligand_docking::MoveMapBuilder() );
		movemap_builder->parse_my_tag( movemap_builder_tag, data );
		data.add( "movemap_builders" , name, movemap_builder);
	}
	TR.flush();
}

std::string MoveMapBuilderLoader::loader_name() { return "MOVEMAP_BUILDERS"; }
std::string MoveMapBuilderLoader::movemap_builder_ct_namer( std::string const & element_name )
{
	return "move_map_builder_" + element_name + "_type";
}
void MoveMapBuilderLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	MoveMapBuilder::provide_xml_schema( xsd );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( MoveMapBuilder::element_name(), & movemap_builder_ct_namer );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.complex_type_naming_func( & movemap_builder_ct_namer )
		.set_subelements_repeatable( subelements )
		.description( "Define a set of MoveMapBuilders and load them into the DataMap for later use by other Movers and Filters" )
		.write_complex_type_to_schema( xsd );


}

parser::DataLoaderOP
MoveMapBuilderLoaderCreator::create_loader() const { return parser::DataLoaderOP( new MoveMapBuilderLoader ); }

std::string
MoveMapBuilderLoaderCreator::keyname() const { return MoveMapBuilderLoader::loader_name(); }

MoveMapBuilderLoaderCreator::DerivedNameFunction
MoveMapBuilderLoaderCreator::schema_ct_naming_function() const
{
	return & MoveMapBuilderLoader::movemap_builder_ct_namer;
}

void MoveMapBuilderLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MoveMapBuilderLoader::provide_xml_schema( xsd );
}


LigandAreaLoader::LigandAreaLoader() {}
LigandAreaLoader::~LigandAreaLoader() = default;

void LigandAreaLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	for ( TagCOP ligand_area_tag : tag->getTags() ) {
		runtime_assert( ligand_area_tag->getName() == LigandArea::element_name() );
		std::string const name( ligand_area_tag->getOption< std::string >( "name" ) );

		if ( data.has("ligand_areas", name) ) {
			TR << "WARNING WARNING ligand_area of name \"" << name
				<< "\" already exists. Skipping\n" << ligand_area_tag << std::endl;
			continue;
		}
		///// Add this movemap to the data map
		protocols::ligand_docking::LigandAreaOP ligand_area( new protocols::ligand_docking::LigandArea() );
		ligand_area->parse_my_tag( ligand_area_tag );
		data.add( "ligand_areas" , name, ligand_area);
	}
	TR.flush();
}

std::string LigandAreaLoader::loader_name() { return "LIGAND_AREAS"; }
std::string LigandAreaLoader::ligand_area_ct_namer( std::string const & element_name )
{
	return "ligand_area_builder_" + element_name + "_type";
}
void LigandAreaLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	LigandArea::provide_xml_schema( xsd );
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( LigandArea::element_name(), & ligand_area_ct_namer );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.complex_type_naming_func( & ligand_area_ct_namer )
		.set_subelements_repeatable( subelements )
		.description( "Define a set of LigandAreas and load them into the DataMap for later use by MoveMapBuilders and other Movers and Filters" )
		.write_complex_type_to_schema( xsd );

}


parser::DataLoaderOP
LigandAreaLoaderCreator::create_loader() const { return parser::DataLoaderOP( new LigandAreaLoader ); }

std::string
LigandAreaLoaderCreator::keyname() const { return LigandAreaLoader::loader_name(); }

LigandAreaLoaderCreator::DerivedNameFunction
LigandAreaLoaderCreator::schema_ct_naming_function() const
{
	return & LigandAreaLoader::ligand_area_ct_namer;
}

void LigandAreaLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LigandAreaLoader::provide_xml_schema( xsd );
}


} //namespace jd2
} //namespace protocols
