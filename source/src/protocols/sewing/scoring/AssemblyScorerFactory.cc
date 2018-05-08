// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   RequirementFactory.cc
/// @brief  Factory for creating Requirements objects
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>

// Package Headers
#include <protocols/sewing/scoring/AssemblyScorerCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>
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
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.AssemblyScorerFactory");
/*
#if defined MULTI_THREADED && defined CXX11
std::atomic< AssemblyScorerFactory * > AssemblyScorerFactory::instance_( 0 );
#else
AssemblyScorerFactory * AssemblyScorerFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex AssemblyScorerFactory::singleton_mutex_;

std::mutex & AssemblyScorerFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif
*/

/*
/// @brief static function to get the instance of ( pointer to) this singleton class
AssemblyScorerFactory * AssemblyScorerFactory::get_instance()
{
boost::function< AssemblyScorerFactory * () > creator = boost::bind( &AssemblyScorerFactory::create_singleton_instance );
utility::thread::safely_create_singleton( creator, instance_ );
return instance_;
}

AssemblyScorerFactory *
AssemblyScorerFactory::create_singleton_instance()
{
return new AssemblyScorerFactory;
}
*/

/// @details Private constructor insures correctness of singleton.
AssemblyScorerFactory::AssemblyScorerFactory() {}

AssemblyScorerFactory::~AssemblyScorerFactory() = default;

void
AssemblyScorerFactory::factory_register(
	AssemblyScorerCreatorCOP creator
) {
	assembly_scorer_types_[ creator->keyname() ] = creator;
}

AssemblyScorerOP
AssemblyScorerFactory::get_assembly_scorer(
	std::string const & type_name
) {
	TR.Trace << "Generating assembly scorer of type " << type_name << std::endl;
	AssemblyScorerCreatorMap::const_iterator iter = assembly_scorer_types_.find( type_name );
	if ( iter != assembly_scorer_types_.end() ) {
		return iter->second->create_assembly_scorer();
	} else {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized AssemblyScorer "
			<< "'" << type_name << "'." << std::endl
			<< "check spelling or "
			<< "register a new AssemblyScorer in the AssemblyScorerFactory" << std::endl
			<< "known AssemblyScorer types are:" << std::endl;

		BOOST_FOREACH ( const AssemblyScorerCreatorMap::value_type& type, assembly_scorer_types_ ) {
			error_msg << "\t" << type.first << std::endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

std::string
AssemblyScorerFactory::assembly_scorer_group_name(){
	return "assembly_scorer";
}


void
AssemblyScorerFactory::define_assembly_scorer_subtag( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//Define the AssemblyScorer group
	define_xml_schema_group( assembly_scorer_types_, assembly_scorer_group_name() , & AssemblyScorerFactory::assembly_scorer_ct_namer, xsd );

	//define the subelement list (only one type)
	XMLSchemaSimpleSubelementList subelements;
	subelements
		//.add_already_defined_subelement( "AssemblyScorer", & assembly_scorer_ct_namer );
		.add_group_subelement( & AssemblyScorerFactory::assembly_scorer_group_name );
	//define the complex type (no attributes, just repeatable subelements
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & movers::AssemblyMover::assembly_mover_subtag_ct_namer )
		.element_name( "AssemblyScorers" )
		.set_subelements_repeatable( subelements )
		.description( "The subtags of this tag define the AssemblyScoreFunction that will be used to evaluate assemblies" )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}
std::string
AssemblyScorerFactory::assembly_scorer_ct_namer( std::string tag_name ){
	return "assembly_scorer_" + tag_name + "_complex_type";
}
//utility::vector1<std::string> AssemblyScorerFactory::get_all_features_names()
//{
// utility::vector1<std::string> collection;
// AssemblyScorerCreatorMap::const_iterator iter = types_.begin();
// while ( iter != types_.end() ) {
//  collection.push_back(iter->first);
//  iter++;
// }
// return collection;

//}

} // namespace
} // namespace
} // namespace
