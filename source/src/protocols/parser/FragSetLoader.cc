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
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/parser/FragSetLoader.hh>
#include <protocols/parser/StandardLoaderCreators.hh>

// for fragments
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <protocols/parser/FragmentReader.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>


static basic::Tracer TR( "protocols.jd2.parser.FragSetLoader" );

namespace protocols {
namespace parser {

FragSetLoader::FragSetLoader() {}
FragSetLoader::~FragSetLoader() {}

void FragSetLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	using protocols::parser::FragmentReader;
	using protocols::parser::FragmentReaderOP;
	typedef std::map< std::string, FragmentReaderOP > FragmentReaderMap;

	FragmentReaderMap frag_readers_map;

	std::string const frag_reader_element_name = FragmentReader::xml_element_name();
	std::string const frag_reader_block_name = "FRAGMENTS";

	// This code says that technically, the FRAGMENTS element isn't required, but it doesn't
	// seem like the fragments that would appear beneath it could possibly be read if it's absent?
	if ( tag->hasTag( frag_reader_block_name ) ) {
		for ( TagCOP subtag : tag->getTag( frag_reader_block_name )->getTags() ) {
			runtime_assert( subtag->getName() == frag_reader_element_name );
			//std::string const name ( subtag->getName() ); // this name is used when fragsets are defined later.
			std::string const name( subtag->getOption< std::string >( "name", "" ));
			runtime_assert( !name.empty() );
			FragmentReaderOP frop( new FragmentReader( subtag ) );
			frag_readers_map[ name ] = frop;
		}
	} else {
		TR.Fatal << "No tag of FRAGMENTS" << std::endl;
		utility_exit_with_message("FRAGMENTS tag is missing for FragSetLoader");
	}

	for ( TagCOP subtag : tag->getTags() ) {
		std::string const name ( subtag->getName() );
		if ( name == "FRAGMENTS" ) continue;

		runtime_assert( name == "FragSet" );

		// The name for the FragmentSet that is about to be defined; requireda
		std::string const set_name( subtag->getOption<std::string>( "name" ) );

		// Comma-separated list of frag readers that contribute to the fragment set
		// Not a very good name for this attribute -- should be something like "frag readers"
		std::string const frag_name ( subtag->getOption<std::string>( "frag_name", "" ) );

		// The file to write the resulting fragment set to -- not required
		std::string const output ( subtag->getOption<std::string>( "output", "" ) );
		runtime_assert( frag_name != "" );

		core::fragment::FragSetOP fragset( new core::fragment::OrderedFragSet );
		utility::vector1< std::string > fnames( utility::string_split( frag_name, ',' ) );
		for ( std::string fname : fnames ) {
			std::map< std::string, FragmentReaderOP >::const_iterator itr;
			itr = frag_readers_map.find( fname );
			if ( itr != frag_readers_map.end() ) {
				FragmentReaderOP frop ( frag_readers_map[ fname ] );
				frop->apply( fragset );
			} else {
				TR.Fatal << "frag_name " << fname << " does not exist." << std::endl;
				utility_exit_with_message("Fragment of given name is not present.");
			}
		}
		runtime_assert( fragset->nr_frames() != 0 );
		data.add( "fragsets", set_name,  fragset );
		// output flagments to fyile
		if ( !output.empty() ) {
			core::fragment::FragmentIO().write_data( output, *fragset );
		}
	}

}

std::string
FragSetLoader::loader_name() {
	return "FRAGSETS";
}

std::string FragSetLoader::frag_set_loader_ct_namer( std::string const & ct_name )
{
	return "frag_set_loader_" + ct_name + "_type";
}

void FragSetLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	FragmentReader::provide_xml_schema( xsd );

	std::string const fragreader_block_element_name = "FRAGMENTS";
	XMLSchemaSimpleSubelementList fragreader_subelements;
	fragreader_subelements.add_already_defined_subelement( FragmentReader::xml_element_name(), & FragSetLoader::frag_set_loader_ct_namer );
	XMLSchemaComplexTypeGenerator fragreader_block_ct;
	fragreader_block_ct.element_name( fragreader_block_element_name )
		.complex_type_naming_func( & FragSetLoader::frag_set_loader_ct_namer )
		.set_subelements_repeatable( fragreader_subelements )
		.description( "The " + fragreader_block_element_name + " element is used to define a set of FragmentReaders, each of which can provide "
		"their fragments to a fragment reader; fragment readers must be composed from one or more FragmentReaders." )
		.write_complex_type_to_schema( xsd );

	AttributeList fragset_attributes;
	fragset_attributes + optional_name_attribute( "The name for the fragment set being defined. This name will be the identifier"
		" used to add the fragment set to the datamap so that other Movers/Filters can access it using this name" )
		+ XMLSchemaAttribute( "frag_name", xs_string, "The comma-separated list of names of the previously defined"
		" FragReader objects that this FragmetSet will be composed from" )
		+ XMLSchemaAttribute( "output", xs_string, "The file to which the FragmentSet should be written, if desired" );

	XMLSchemaSimpleSubelementList fragreader_block_subelement;
	fragreader_block_subelement.add_already_defined_subelement( fragreader_block_element_name,
		& FragSetLoader::frag_set_loader_ct_namer );

	XMLSchemaSimpleSubelementList fragset_subelements;
	fragset_subelements.add_simple_subelement( "FragSet", fragset_attributes, "The FragSet element will"
		" be added to the datamap and can be accessed thereafter by other Movers and Filters");

	XMLSchemaComplexTypeGenerator fragset_ct;
	fragset_ct.element_name( loader_name() )
		.complex_type_naming_func( & frag_set_loader_ct_namer )
		.description( "The FRAGMENTS element can be used to load FragmentSets into the datamap; it consists of two parts:"
		" first, a list of FragmentReaders that load fragments from various sources (the vall, a fragment file, or from a particular"
		" PDB file), and second a list of FragmentSets that are composed from the FragmentReaders that have already been defined." )
		.add_ordered_subelement_set_as_optional( fragreader_block_subelement )
		.add_ordered_subelement_set_as_repeatable( fragset_subelements )
		.write_complex_type_to_schema( xsd );
}

DataLoaderOP
FragSetLoaderCreator::create_loader() const { return DataLoaderOP( new FragSetLoader ); }

std::string
FragSetLoaderCreator::keyname() const {
	return FragSetLoader::loader_name();
}

FragSetLoaderCreator::DerivedNameFunction
FragSetLoaderCreator::schema_ct_naming_function() const
{
	return & FragSetLoader::frag_set_loader_ct_namer;
}

void
FragSetLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FragSetLoader::provide_xml_schema( xsd );
}


} //namespace parser
} //namespace protocols
