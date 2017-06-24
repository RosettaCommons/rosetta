// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsDefinerLoader.cc
/// @brief  Implementation the LoopsDefinerLoader class which implements the DataLoader interface
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefinerLoader.hh>
#include <protocols/loops/loops_definers/LoopsDefinerCreator.hh>
#include <protocols/loops/loops_definers/LoopsDefinerLoaderCreator.hh>

// Project Headers
#include <protocols/loops/loops_definers/LoopsDefiner.hh>
#include <protocols/loops/loops_definers/LoopsDefinerFactory.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.hh>


// Utility headers
#include <utility/tag/Tag.hh>

// Boost Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using std::string;
using std::endl;
using core::pose::Pose;
using utility::tag::TagCOP;
using basic::datacache::DataMap;
using utility::vector0;
using protocols::loops::loops_definers::LoopsDefinerOP;
using protocols::loops::loops_definers::LoopsDefinerFactory;

namespace protocols {
namespace loops {
namespace loops_definers {

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loops_definers.LoopsDefinerLoader" );

LoopsDefinerLoader::LoopsDefinerLoader() {}
LoopsDefinerLoader::~LoopsDefinerLoader() {}

void LoopsDefinerLoader::load_data(
	Pose const & pose,
	TagCOP const tag,
	basic::datacache::DataMap & data
) const
{
	for ( TagCOP subtag : tag->getTags() ) {
		string const & type( subtag->getName() );
		if ( ! subtag->hasOption("name") ) {
			utility_exit_with_message( "Can't create unnamed Loops definition (type: " + type + ")" );
		}
		string const & name( subtag->getOption<string>("name") );
		if ( data.has( "loops_definers", name ) ) {
			TR.Error << "LoopsDefiner of name \"" << name
				<< "\" (with type " << type << ") already exists. \n" << subtag << endl;
			utility_exit_with_message("Duplicate definition of LoopsDefiner with name " + name);
		}
		LoopsDefinerOP loops_definer( LoopsDefinerFactory::get_instance()->create_loops_definer( type ) );
		loops_definer->parse_my_tag(subtag, data, pose);
		data.add("loops_definers", name, loops_definer );
		TR << "Created LoopsDefiner named \"" << name << "\" of type " << type << endl;
	}
	TR.flush();
}

std::string LoopsDefinerLoader::loader_name() { return "LOOP_DEFINITIONS"; }
std::string LoopsDefinerLoader::loop_def_loader_ct_namer( std::string const & element_name )
{
	return "loop_def_loader_" + element_name + "_type";
}

void LoopsDefinerLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	LoopsDefinerFactory::get_instance()->define_loop_definer_xml_schema( xsd );

	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & LoopsDefinerFactory::loop_definer_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( loader_name() )
		.complex_type_naming_func( & loop_def_loader_ct_namer )
		.set_subelements_repeatable( subelements )
		.description( "The " + loader_name() + " element allows you to define loops from various sources (e.g. from files, from databases, or directly within the XML file); loops will be named and placed in the DataMap for other Movers and Filters to retrieve." )
		.write_complex_type_to_schema( xsd );
}


parser::DataLoaderOP
LoopsDefinerLoaderCreator::create_loader() const { return parser::DataLoaderOP( new LoopsDefinerLoader ); }

string
LoopsDefinerLoaderCreator::keyname() const { return LoopsDefinerLoader::loader_name(); }

LoopsDefinerLoaderCreator::DerivedNameFunction
LoopsDefinerLoaderCreator::schema_ct_naming_function() const
{
	return & LoopsDefinerLoader::loop_def_loader_ct_namer;
}

void
LoopsDefinerLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopsDefinerLoader::provide_xml_schema( xsd );
}

} //namespace
} //namespace
} //namespace
