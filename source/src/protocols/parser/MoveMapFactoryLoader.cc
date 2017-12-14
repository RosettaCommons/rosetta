// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/MoveMapFactoryLoader.cc
/// @brief  Implementation of the XML parser's DataLoader class for MoveMaps
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/parser/MoveMapFactoryLoader.hh>
#include <protocols/parser/MoveMapFactoryLoaderCreator.hh>


// Project Headers
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/movemap/util.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.MoveMapFactoryLoader" );

MoveMapFactoryLoader::MoveMapFactoryLoader() = default;
MoveMapFactoryLoader::~MoveMapFactoryLoader() = default;

void MoveMapFactoryLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace core::select::movemap;

	for ( utility::tag::TagCOP const & subtag : tag->getTags() ) {
		if ( ! subtag->hasOption("name") ) {
			utility_exit_with_message( "Can't create unnamed MoveMap" );
		}
		std::string const & name( subtag->getOption<std::string>("name") );
		MoveMapFactoryOP new_mmf( new MoveMapFactory );
		runtime_assert( new_mmf != nullptr );
		new_mmf->parse_my_tag( subtag, data );
		data.add( mmf_cat_in_datamap(), name, new_mmf );
		TR << "Defined MoveMap named \"" << name << "\"" << std::endl;
	}
	TR.flush();

}

std::string MoveMapFactoryLoader::mmf_cat_in_datamap() { return core::select::movemap::movemap_factory_category(); }

std::string MoveMapFactoryLoader::loader_name() { return "MOVE_MAP_FACTORIES"; }

std::string MoveMapFactoryLoader::mmf_loader_ct_namer( std::string const & element_name )
{
	return "movemap_factory_loader_" + element_name + "_type";
}

void MoveMapFactoryLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using core::select::movemap::MoveMapFactory;

	MoveMapFactory::provide_xml_schema( xsd );

	XMLSchemaSimpleSubelementList mmf_subelements;
	mmf_subelements.add_already_defined_subelement(
		MoveMapFactory::element_name(),
		& MoveMapFactory::mmf_ct_namer
	);

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.description( "Define MoveMapFactories that can be used to generate MoveMaps just in time, tailored to the conformation of the Pose" )
		.complex_type_naming_func( mmf_loader_ct_namer )
		.set_subelements_repeatable( mmf_subelements )
		.write_complex_type_to_schema( xsd );
}


DataLoaderOP
MoveMapFactoryLoaderCreator::create_loader() const { return DataLoaderOP( new MoveMapFactoryLoader ); }

std::string
MoveMapFactoryLoaderCreator::keyname() const { return MoveMapFactoryLoader::loader_name(); }

MoveMapFactoryLoaderCreator::DerivedNameFunction
MoveMapFactoryLoaderCreator::schema_ct_naming_function() const
{
	return & MoveMapFactoryLoader::mmf_loader_ct_namer;
}

void
MoveMapFactoryLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MoveMapFactoryLoader::provide_xml_schema( xsd );
}


} //namespace parser
} //namespace protocols
