// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryLoader.cc
/// @brief  Implementation of the XML parser's ChemistryLoader
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/chemistries/ChemistryLoader.hh>
#include <protocols/chemistries/ChemistryLoaderCreator.hh>

// Project headers
#include <protocols/chemistries/Chemistry.hh>
#include <protocols/chemistries/ChemistryFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace chemistries {

static basic::Tracer TR( "protocols.chemistries.ChemistryLoader" );

ChemistryLoader::ChemistryLoader() {}
ChemistryLoader::~ChemistryLoader() {}

void ChemistryLoader::load_data(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	using namespace utility::tag;
	using namespace protocols::chemistries;
	typedef utility::vector0< TagCOP > TagCOPs;

	TagCOPs const & chemistry_tags( tag->getTags() );
	for ( core::Size ii = 0; ii < chemistry_tags.size(); ++ii ) {
		TagCOP ii_tag = chemistry_tags[ ii ];
		ChemistryOP chemistry = ChemistryFactory::get_instance()->new_chemistry(
			ii_tag,
			datamap
		);

		// if "name" is specified, add it to the data map under that name. Otherwise use the type name.
		bool const data_add_status = datamap.add( "chemistry" ,
			ii_tag->getOption( "name", ii_tag->getName() ),
			chemistry );
		if ( !data_add_status ) {
			utility_exit_with_message( "Chemistry '" + ii_tag->getName() + "' already exists in the basic::datacache::DataMap. Please rename." );
		}
	}
}

std::string
ChemistryLoader::loader_name() { return "CHEMISTRY"; }

std::string
ChemistryLoader::chemistry_loader_ct_namer( std::string const & element_name )
{
	return "chemistry_loader_" + element_name + "_type";
}

void ChemistryLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace protocols::chemistries;

	ChemistryFactory::get_instance()->define_chemistry_xml_schema( xsd );

	XMLSchemaSimpleSubelementList chem_loader_subelements;
	chem_loader_subelements.add_group_subelement( & ChemistryFactory::chemistry_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator chem_ct;
	chem_ct.element_name( loader_name() ).complex_type_naming_func( & chemistry_loader_ct_namer )
		.description( "Chemistry objects may be defined as subelements of the " + loader_name() + " element, and then will be placed into the DataMap"
		" for later retrieval by Movers and Filters or anything else that might use a Chemistry. All immediate subelements should have the 'name' attribute"
		" as that is how they will be identified in the DataMap. " )
		.set_subelements_repeatable( chem_loader_subelements )
		.write_complex_type_to_schema( xsd );

}

parser::DataLoaderOP
ChemistryLoaderCreator::create_loader() const { return parser::DataLoaderOP( new ChemistryLoader ); }

std::string
ChemistryLoaderCreator::keyname() const { return ChemistryLoader::loader_name(); }

ChemistryLoaderCreator::DerivedNameFunction
ChemistryLoaderCreator::schema_ct_naming_function() const
{
	return & ChemistryLoader::chemistry_loader_ct_namer;
}

void
ChemistryLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ChemistryLoader::provide_xml_schema( xsd );
}

} //namespace chemistries
} //namespace protocols
