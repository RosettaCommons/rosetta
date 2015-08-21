// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RotamerRecoveryFactory.cc
/// @brief  Factory for creating RotamerRecoverys objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RRProtocol.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryCreator.hh>

// Package Headers
#include <basic/Tracer.hh>

// Project Headers

// C++ Headers
#include <string>
#include <sstream>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>

// Singleton instance and mutex static data members
namespace utility {

using protocols::rotamer_recovery::RotamerRecoveryFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< RotamerRecoveryFactory >::singleton_mutex_{};
template <> std::atomic< RotamerRecoveryFactory * > utility::SingletonBase< RotamerRecoveryFactory >::instance_( 0 );
#else
template <> RotamerRecoveryFactory * utility::SingletonBase< RotamerRecoveryFactory >::instance_( 0 );
#endif

}

namespace protocols {
namespace rotamer_recovery {

using std::endl;
using std::string;
using std::stringstream;
using core::pose::Pose;
using core::scoring::ScoreFunctionCOP;

static thread_local basic::Tracer tr( "protocols.rotamer_recovery.RotamerRecoveryFactory" );

RotamerRecoveryFactory *
RotamerRecoveryFactory::create_singleton_instance()
{
	return new RotamerRecoveryFactory;
}

/// @details Private constructor insures correctness of singleton.
RotamerRecoveryFactory::RotamerRecoveryFactory() {}

RotamerRecoveryFactory::~RotamerRecoveryFactory() {}

/*
20140630: Additional factory_register() methods disabled because of ambiguouity below.

src/utility/factory/WidgetRegistrator.hh:37:28: error: call to member function 'factory_register' is ambiguous
FACTORY::get_instance()->factory_register( CREATOROP( new CREATOR ) );
~~~~~~~~~~~~~~~~~~~~~~~~~^~~~~~~~~~~~~~~~
src/protocols/rotamer_recovery/RotamerRecoveryFactory.hh:122:33: note: in instantiation of member function
'utility::factory::WidgetRegistrator<protocols::rotamer_recovery::RotamerRecoveryFactory, protocols::rotamer_recovery::RRReporterHumanCreator>::WidgetRegistrator' requested here
RotamerRecoveryRegistrator() : parent() {}
^
src/protocols/init/init.RotamerRecoveryRegistrators.ihh:32:61: note: in instantiation of member function
'protocols::rotamer_recovery::RotamerRecoveryRegistrator<protocols::rotamer_recovery::RRReporterHumanCreator>::RotamerRecoveryRegistrator' requested here
static RotamerRecoveryRegistrator< RRReporterHumanCreator > RRReporterHumanCreator_registrator;
^
src/protocols/rotamer_recovery/RotamerRecoveryFactory.hh:66:7: note: candidate function
void factory_register( RRProtocolCreatorCOP creator );
^
src/protocols/rotamer_recovery/RotamerRecoveryFactory.hh:67:7: note: candidate function
void factory_register( RRComparerCreatorCOP creator );
^
src/protocols/rotamer_recovery/RotamerRecoveryFactory.hh:68:7: note: candidate function
void factory_register( RRReporterCreatorCOP creator );
*/

void
RotamerRecoveryFactory::factory_register( utility::pointer::ReferenceCountOP creator ) {

	{
		RRProtocolCreatorCOP p( utility::pointer::dynamic_pointer_cast< RRProtocolCreator const >( creator ) );
		if ( p ) {
			protocol_types_[ p->type_name() ] = p;
		}
	}

	{
		RRComparerCreatorCOP p( utility::pointer::dynamic_pointer_cast< RRComparerCreator const >( creator ) );
		if ( p ) {
			comparer_types_[ p->type_name() ] = p;
		}
	}

	{
		RRReporterCreatorCOP p( utility::pointer::dynamic_pointer_cast< RRReporterCreator const >( creator ) );
		if ( p ) {
			reporter_types_[ p->type_name() ] = p;
		}
	}

}

/*
void
RotamerRecoveryFactory::factory_register(
RRProtocolCreatorCOP creator
) {
protocol_types_[ creator->type_name() ] = creator;
}

void
RotamerRecoveryFactory::factory_register(
RRComparerCreatorCOP creator
) {
comparer_types_[ creator->type_name() ] = creator;
}

void
RotamerRecoveryFactory::factory_register(
RRReporterCreatorCOP creator
) {
reporter_types_[ creator->type_name() ] = creator;
}
*/

RRProtocolOP
RotamerRecoveryFactory::get_rotamer_recovery_protocol(
	string const & type_name
) {
	tr.Trace << "get rotamer recovery protocol of type '" << type_name << "'" << endl;
	RRProtocolCreatorMap::const_iterator iter = protocol_types_.find( type_name );
	if ( iter != protocol_types_.end() ) {
		return iter->second->create_protocol();
	} else {

		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized rotamer recovery protocol "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new RRProtocol with the RotamerRecoveryFactory" << endl
			<< "known RRProtocol types are:" << endl;

		BOOST_FOREACH ( const RRProtocolCreatorMap::value_type& type, protocol_types_ ) {
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

RRComparerOP
RotamerRecoveryFactory::get_rotamer_recovery_comparer(
	string const & type_name
) {
	tr.Trace << "get rotamer recovery comparer of type '" << type_name << "'" << endl;
	RRComparerCreatorMap::const_iterator iter = comparer_types_.find( type_name );
	if ( iter != comparer_types_.end() ) {
		return iter->second->create_comparer();
	} else {

		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized rotamer recovery comparer "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new RRComparer with the RotamerRecoveryFactory" << endl
			<< "known RRComparer types are:" << endl;

		BOOST_FOREACH ( const RRComparerCreatorMap::value_type& type, comparer_types_ ) {
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

RRReporterOP
RotamerRecoveryFactory::get_rotamer_recovery_reporter(
	string const & type_name
) {
	tr.Trace << "get rotamer recovery reporter of type '" << type_name << "'" << endl;
	RRReporterCreatorMap::const_iterator iter = reporter_types_.find( type_name );
	if ( iter != reporter_types_.end() ) {
		return iter->second->create_reporter();
	} else {

		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized rotamer recovery reporter "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new RRReporter with the RotamerRecoveryFactory" << endl
			<< "known RRReporter types are:" << endl;

		BOOST_FOREACH ( const RRReporterCreatorMap::value_type& type, reporter_types_ ) {
			error_msg << "\t" << type.first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
	return 0;
}

RotamerRecoveryOP
RotamerRecoveryFactory::get_rotamer_recovery(
	string const & protocol_name,
	string const & comparer_name,
	string const & reporter_name
) {
	return RotamerRecoveryOP( new RotamerRecovery(
		get_rotamer_recovery_protocol(protocol_name),
		get_rotamer_recovery_comparer(comparer_name),
		get_rotamer_recovery_reporter(reporter_name)) );
}

} // namespace
} // namespace
