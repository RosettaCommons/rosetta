// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/parser/PackerPaletteLoader.cc
/// @brief  XML parsing for PackerPalettes.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/parser/PackerPaletteLoader.hh>
#include <protocols/parser/PackerPaletteLoaderCreator.hh>

// Project headers
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/palette/PackerPaletteFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace parser {

static basic::Tracer TR( "protocols.parser.PackerPaletteLoader" );

PackerPaletteLoader::PackerPaletteLoader() {}
PackerPaletteLoader::~PackerPaletteLoader() {}

void PackerPaletteLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	using namespace utility::tag;
	using core::pack::palette::PackerPaletteOP;
	typedef utility::vector0< TagCOP > TagCOPs;

	TagCOPs const & selector_tags( tag->getTags() );
	for ( core::Size ii = 0; ii < selector_tags.size(); ++ii ) {
		TagCOP ii_tag = selector_tags[ ii ];
		PackerPaletteOP selector = core::pack::palette::PackerPaletteFactory::get_instance()->newPackerPalette(
			ii_tag->getName(),
			datamap,
			ii_tag
		);

		// if "name" is specified, add it to the data map under that name. Otherwise use the type name.
		bool const data_add_status = datamap.add( "packer_palettes" ,
			ii_tag->getOption( "name", ii_tag->getName() ),
			selector );
		if ( !data_add_status ) {
			utility_exit_with_message( "PackerPalette '" + ii_tag->getName() + "' already exists in the basic::datacache::DataMap. Please rename." );
		}
	}
	TR.flush();
}

std::string PackerPaletteLoaderCreator::packerpalette_loader_ct_namer( std::string const & element_name )
{
	return "packer_palette_loader_" + element_name + "_type";
}

DataLoaderOP
PackerPaletteLoaderCreator::create_loader() const { return DataLoaderOP( new PackerPaletteLoader ); }

std::string
PackerPaletteLoaderCreator::keyname() const { return PackerPaletteLoader::loader_name(); }

DataLoaderCreator::DerivedNameFunction
PackerPaletteLoaderCreator::schema_ct_naming_function() const {
	return & PackerPaletteLoaderCreator::packerpalette_loader_ct_namer;
}


std::string
PackerPaletteLoader::loader_name() {
	return "PACKER_PALETTES";
}

void
PackerPaletteLoader::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {

	using namespace utility::tag;

	core::pack::palette::PackerPaletteFactory::get_instance()->define_packer_palette_xml_schema( xsd );

	XMLSchemaSimpleSubelementList pp_subelements;
	pp_subelements.add_group_subelement( & core::pack::palette::PackerPaletteFactory::packer_palette_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.description( "Define PackerPalettes, which establish the set of residue types (the \"palette\") that one might use in design.  PackerPalettes define the total set of residue types available, while TaskOperations prune the types allowed at a particular position or during a particular design step." )
		.complex_type_naming_func( PackerPaletteLoaderCreator::packerpalette_loader_ct_namer )
		.set_subelements_repeatable( pp_subelements )
		.write_complex_type_to_schema( xsd );

}

void PackerPaletteLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const { PackerPaletteLoader::provide_xml_schema( xsd ); }

} //namespace parser
} //namespace protocols
