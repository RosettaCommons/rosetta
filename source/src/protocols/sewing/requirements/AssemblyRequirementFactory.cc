// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   AssemblyRequirementFactory.cc
/// @brief  Factory for creating Requirements objects
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>

// Package Headers
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <protocols/sewing/movers/AssemblyMover.hh>
// Project Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

//Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>



namespace protocols {
namespace sewing  {
namespace requirements {

static basic::Tracer TR("protocols.sewing.requirements.AssemblyRequirementFactory");
/*
#if defined MULTI_THREADED && defined CXX11
std::atomic< AssemblyRequirementFactory * > AssemblyRequirementFactory::instance_( 0 );
#else
AssemblyRequirementFactory * AssemblyRequirementFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex AssemblyRequirementFactory::singleton_mutex_;

std::mutex & AssemblyRequirementFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif
/// @brief static function to get the instance of ( pointer to) this singleton class
AssemblyRequirementFactory * AssemblyRequirementFactory::get_instance()
{
boost::function< AssemblyRequirementFactory * () > creator = boost::bind( &AssemblyRequirementFactory::create_singleton_instance );
utility::thread::safely_create_singleton( creator, instance_ );
return instance_;
}

AssemblyRequirementFactory *
AssemblyRequirementFactory::create_singleton_instance()
{
return new AssemblyRequirementFactory;
}
*/
/// @details Private constructor insures correctness of singleton.
AssemblyRequirementFactory::AssemblyRequirementFactory() {}

AssemblyRequirementFactory::~AssemblyRequirementFactory() = default;

void
AssemblyRequirementFactory::factory_register(
	AssemblyRequirementCreatorCOP creator
) {
	requirement_types_[ creator->keyname() ] = creator;
}

AssemblyRequirementOP
AssemblyRequirementFactory::get_requirement(
	std::string const & type_name
) {
	TR.Trace << "Generating requirement of type " << type_name << std::endl;
	AssemblyRequirementCreatorMap::const_iterator iter = requirement_types_.find( type_name );
	if ( iter != requirement_types_.end() ) {
		return iter->second->create_requirement();
	} else {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized AssemblyRequirement "
			<< "'" << type_name << "'." << std::endl
			<< "check spelling or "
			<< "register a new AssemblyRequirement in the AssemblyRequirementFactory" << std::endl
			<< "known AssemblyRequirement types are:" << std::endl;

		BOOST_FOREACH ( const AssemblyRequirementCreatorMap::value_type& type, requirement_types_ ) {
			error_msg << "\t" << type.first << std::endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

std::string
AssemblyRequirementFactory::assembly_requirement_ct_namer( std::string tag_name ){
	return "assembly_requirement_" + tag_name + "_complex_type";
}


std::string
AssemblyRequirementFactory::assembly_requirement_group_name(){
	return "assembly_requirement";
}

void
AssemblyRequirementFactory::define_assembly_requirement_subtag( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//Define the AssemblyRequirement group
	define_xml_schema_group( requirement_types_, assembly_requirement_group_name() , & AssemblyRequirementFactory::assembly_requirement_ct_namer, xsd );

	//define the subelement list (only one type)
	XMLSchemaSimpleSubelementList subelements;
	subelements
		//.add_already_defined_subelement( "AssemblyRequirement", & assembly_scorer_ct_namer );
		.add_group_subelement( & assembly_requirement_group_name );
	//define the complex type (no attributes, just repeatable subelements
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & movers::AssemblyMover::assembly_mover_subtag_ct_namer )
		.element_name( "AssemblyRequirements" )
		.set_subelements_repeatable( subelements )
		.description( "Subtags of this tag define the set of requirements that will be used when evaluating SEWING assemblies" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}




//utility::vector1<std::string> AssemblyRequirementFactory::get_all_features_names()
//{
// utility::vector1<std::string> collection;
// AssemblyRequirementCreatorMap::const_iterator iter = types_.begin();
// while ( iter != types_.end() ) {
//  collection.push_back(iter->first);
//  iter++;
// }
// return collection;

//}

} // namespace
} // namespace
} // namespace
