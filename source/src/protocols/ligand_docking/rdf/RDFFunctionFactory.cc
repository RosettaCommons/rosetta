// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/rdf/RDFFunctionFactory.cc
/// @brief  Factory for creating RDFBase objects
/// @author Sam DeLuca

// Unit Headers
#include <protocols/ligand_docking/rdf/RDFFunctionFactory.hh>
#include <protocols/ligand_docking/rdf/RDFFunctionCreator.hh>
#include <protocols/ligand_docking/rdf/RDFBase.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector0.hh>

// Boost Headers
#include <boost/foreach.hpp>

// C++ Headers
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>





namespace protocols {
namespace ligand_docking {
namespace rdf {
	
using std::endl;
using std::string;
using std::stringstream;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::tag::TagCOP;

static basic::Tracer tr("protocols.ligand_docking.rdf.RDFFunctionFactory");

RDFFunctionFactory * RDFFunctionFactory::instance_( 0 );

/// @details Private constructor insures correctness of singleton.
RDFFunctionFactory::RDFFunctionFactory() {}

RDFFunctionFactory::RDFFunctionFactory(
	const RDFFunctionFactory &
) {}

RDFFunctionFactory::~RDFFunctionFactory() {}


RDFFunctionFactory *
RDFFunctionFactory::get_instance()
{
	if ( instance_ == 0 ) {
		instance_ = new RDFFunctionFactory;
	}
	return instance_;
}


void
RDFFunctionFactory::factory_register(
	RDFFunctionCreatorCOP creator
) {
	types_[ creator->type_name() ] = creator;
}


RDFBaseOP
RDFFunctionFactory::get_rdf_function(std::string const & type_name)
{
	tr.Trace << "generate RDF function of type " << type_name << std::endl;
	RDFFunctionCreatorMap::const_iterator iter = types_.find( type_name );
	if (iter != types_.end()) {
		return iter->second->create_rdf_function();
	} else {
		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized RDF Function "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new RDF Function in the RDFFunctionFactory" << endl
			<< "known RDF Function types are:" << endl;

		BOOST_FOREACH(const RDFFunctionCreatorMap::value_type& type, types_){
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}
utility::vector1<std::string> RDFFunctionFactory::get_all_function_names()
{
	utility::vector1<std::string> collection;
	RDFFunctionCreatorMap::const_iterator iter = types_.begin();
	while ( iter != types_.end() ) {
		collection.push_back(iter->first);
		iter++;
	}
	return collection;

}
RDFBaseOP
RDFFunctionFactory::get_rdf_function(
	TagCOP const tag,
	basic::datacache::DataMap & data
) {
	assert(tag->getName() == "RDF");

	string type_name;
	if(!tag->hasOption("name")){
		utility_exit_with_message("'RDF' tags require a name field");
	} else {
		type_name = tag->getOption<string>("name");
	}

	RDFBaseOP rdf_function(get_rdf_function(type_name));

	rdf_function->parse_my_tag(tag, data);
	return rdf_function;
}
	
}
} // namespace
} // namespace
