// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loops_definers/LoopsDefinerFactory.cc
/// @brief  Factory for creating LoopsDefiner objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefinerFactory.hh>
#include <protocols/loops/loops_definers/LoopsDefinerCreator.hh>
#include <protocols/loops/loops_definers/LoopsDefiner.hh>

// Package Headers
#include <basic/Tracer.hh>

// Project Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>

// C++ Headers
#include <string>
#include <sstream>

// Singleton instance and mutex static data members
namespace utility {

using protocols::loops::loops_definers::LoopsDefinerFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< LoopsDefinerFactory >::singleton_mutex_{};
template <> std::atomic< LoopsDefinerFactory * > utility::SingletonBase< LoopsDefinerFactory >::instance_( 0 );
#else
template <> LoopsDefinerFactory * utility::SingletonBase< LoopsDefinerFactory >::instance_( 0 );
#endif

}

namespace protocols {
namespace loops {
namespace loops_definers {

using std::endl;
using std::string;
using std::stringstream;
using utility::vector1;


static thread_local basic::Tracer tr( "protocols.loops.loops_definers.LoopsDefinerFactory" );

LoopsDefinerFactory *
LoopsDefinerFactory::create_singleton_instance()
{
	return new LoopsDefinerFactory;
}

/// @details Private constructor insures correctness of singleton.
LoopsDefinerFactory::LoopsDefinerFactory() {}

LoopsDefinerFactory::~LoopsDefinerFactory() {}

void
LoopsDefinerFactory::factory_register(
	LoopsDefinerCreatorOP creator
) {
	types_[ creator->type_name() ] = creator;
}

bool
LoopsDefinerFactory::has_type (
	string const & type_name
) const {
	LoopsDefinerCreatorMap::const_iterator iter = types_.find( type_name );
	return iter != types_.end();
}

LoopsDefinerOP
LoopsDefinerFactory::create_loops_definer(
	std::string const & type_name
) {

	tr.Trace << "generate LoopsDefiner of type " << type_name << std::endl;
	LoopsDefinerCreatorMap::const_iterator iter = types_.find( type_name );
	if (iter != types_.end()) {
		return iter->second->create_loops_definer();
	} else {
		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized LoopsDefiner "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new LoopsDefiner with the LoopsDefinerFactory" << endl
			<< "known LoopsDefiner types are:" << endl;

		BOOST_FOREACH(const LoopsDefinerCreatorMap::value_type& type, types_){
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

vector1< string >
LoopsDefinerFactory::get_all_loops_definer_names(
) const {
	vector1< string > collection;
	LoopsDefinerCreatorMap::const_iterator iter = types_.begin(), end = types_.end();
	while ( iter != end ) {
		collection.push_back(iter->first);
		++iter;
	}
	return collection;
}


} // namespace
} // namespace
} // namespace
