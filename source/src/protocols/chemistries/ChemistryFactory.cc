// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/ChemistryFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary Chemistry objects
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/chemistries/ChemistryFactory.hh>

// Package headers
#include <protocols/chemistries/Chemistry.hh>
#include <protocols/chemistries/ChemistryCreator.hh>
#include <protocols/chemistries/util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace chemistries {

static basic::Tracer TR("protocols.chemistries.ChemistryFactory");

void
ChemistryFactory::factory_register( ChemistryCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string err_msg = "Factory Name Conflict: Two or more ChemistryCreators registered with the name " + creator->keyname();
		if ( throw_on_double_registration_ ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
		} else {
			utility_exit_with_message(  err_msg );
		}
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool ChemistryFactory::has_type( std::string const & chemistry_type ) const
{
	return creator_map_.find( chemistry_type ) != creator_map_.end();
}

ChemistryOP ChemistryFactory::new_chemistry(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	std::string const & chemistry_name = tag->getName();

	if ( ! has_type( chemistry_name ) ) {

		std::string err_msg =  "No ChemistryCreator with the name '" + chemistry_name + "' has been registered with the ChemistryFactory";
		TR.Error << err_msg << std::endl;
		TR.Error << "Known Chemistries:";
		// TODO: Probably should alphabetize this at some point.
		for ( std::map< std::string, ChemistryCreatorOP >::const_iterator itr( creator_map_.begin() ); itr != creator_map_.end(); ++itr ) {
			TR.Error << " " << itr->first << ", ";
		}
		TR.Error << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	std::map< std::string, ChemistryCreatorOP >::const_iterator iter = creator_map_.find( chemistry_name );
	ChemistryOP new_chemistry = iter->second->create_chemistry();
	new_chemistry->parse_my_tag( tag, datamap );
	return new_chemistry;
}

void ChemistryFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}

void ChemistryFactory::define_chemistry_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			creator_map_,
			chemistry_xml_schema_group_name(),
			& complex_type_name_for_chemistry,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Could not generate an XML Schema for Chemistry from ChemistryFactory; offending class"
			" must call protocols::chemistries::complex_type_name_for_chemistry when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::string ChemistryFactory::chemistry_xml_schema_group_name()
{
	return "chemistry";
}

ChemistryFactory::ChemistryFactory() :
	throw_on_double_registration_( false )
{}


} //namespace chemistries
} //namespace protocols
