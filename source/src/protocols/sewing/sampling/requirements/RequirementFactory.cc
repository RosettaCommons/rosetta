// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   RequirementFactory.cc
/// @brief  Factory for creating Requirements objects
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/requirements/RequirementFactory.hh>
#include <protocols/sewing/sampling/requirements/GlobalRequirement.hh>
#include <protocols/sewing/sampling/requirements/IntraSegmentRequirement.hh>
#include <protocols/sewing/sampling/requirements/RequirementCreator.hh>

// Package Headers
#include <basic/Tracer.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers
#include <sstream>

//Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>



namespace protocols {
namespace sewing  {
namespace sampling {
namespace requirements {

static basic::Tracer tr("protocols.sewing.sampling.requirements.RequirementFactory");

#if defined MULTI_THREADED && defined CXX11
std::atomic< RequirementFactory * > RequirementFactory::instance_( 0 );
#else
RequirementFactory * RequirementFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex RequirementFactory::singleton_mutex_;

std::mutex & RequirementFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
RequirementFactory * RequirementFactory::get_instance()
{
	boost::function< RequirementFactory * () > creator = boost::bind( &RequirementFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

RequirementFactory *
RequirementFactory::create_singleton_instance()
{
	return new RequirementFactory;
}

/// @details Private constructor insures correctness of singleton.
RequirementFactory::RequirementFactory() {}

RequirementFactory::RequirementFactory(
	const RequirementFactory &
) {}

RequirementFactory::~RequirementFactory() {}

void
RequirementFactory::factory_register(
	GlobalRequirementCreatorCOP creator
) {
	global_types_[ creator->type_name() ] = creator;
}

void
RequirementFactory::factory_register(
	IntraSegmentRequirementCreatorCOP creator
) {
	intra_segment_types_[ creator->type_name() ] = creator;
}


GlobalRequirementOP
RequirementFactory::get_global_requirement(
	std::string const & type_name
) {
	tr.Trace << "Generating global requirement of type " << type_name << std::endl;
	GlobalRequirementCreatorMap::const_iterator iter = global_types_.find( type_name );
	if ( iter != global_types_.end() ) {
		return iter->second->create_requirement();
	} else {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized GlobalRequirement "
			<< "'" << type_name << "'." << std::endl
			<< "check spelling or "
			<< "register a new Requirement in the RequirementFactory" << std::endl
			<< "known Requirement types are:" << std::endl;

		BOOST_FOREACH ( const GlobalRequirementCreatorMap::value_type& type, global_types_ ) {
			error_msg << "\t" << type.first << std::endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}


IntraSegmentRequirementOP
RequirementFactory::get_intra_segment_requirement(
	std::string const & type_name
) {
	tr.Trace << "Generating intra-segment requirement of type " << type_name << std::endl;
	IntraSegmentRequirementCreatorMap::const_iterator iter = intra_segment_types_.find( type_name );
	if ( iter != intra_segment_types_.end() ) {
		return iter->second->create_requirement();
	} else {
		std::stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized IntraSegmentRequirement "
			<< "'" << type_name << "'." << std::endl
			<< "check spelling or "
			<< "register a new Requirement in the RequirementFactory" << std::endl
			<< "known Requirement types are:" << std::endl;

		BOOST_FOREACH ( const IntraSegmentRequirementCreatorMap::value_type& type, intra_segment_types_ ) {
			error_msg << "\t" << type.first << std::endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

//utility::vector1<std::string> RequirementFactory::get_all_features_names()
//{
// utility::vector1<std::string> collection;
// RequirementCreatorMap::const_iterator iter = types_.begin();
// while ( iter != types_.end() ) {
//  collection.push_back(iter->first);
//  iter++;
// }
// return collection;
//}

} // namespace
} // namespace
} // namespace
} // namespace
