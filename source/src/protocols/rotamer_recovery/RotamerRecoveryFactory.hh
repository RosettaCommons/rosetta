// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RotamerRecoveryFactory.hh
/// @brief Factory for creating RotamerRecovery objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_rotamer_recovery_RotamerRecoveryFactory_hh
#define INCLUDED_protocols_rotamer_recovery_RotamerRecoveryFactory_hh

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.fwd.hh>

// Project Headers
#include <protocols/rotamer_recovery/RotamerRecovery.fwd.hh>
#include <protocols/rotamer_recovery/RRProtocol.fwd.hh>
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryCreator.fwd.hh>

// Platform Headers
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace rotamer_recovery {

/// Create Rotamer_Recovery Reporters
class RotamerRecoveryFactory {
private:
	// Private constructor to make it singleton managed
	RotamerRecoveryFactory();
	RotamerRecoveryFactory(const RotamerRecoveryFactory & src); // unimplemented

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static RotamerRecoveryFactory * create_singleton_instance();

	RotamerRecoveryFactory const &
	operator=( RotamerRecoveryFactory const & ); // unimplemented

public:

	// Warning this is not called because of the singleton pattern
	virtual ~RotamerRecoveryFactory();

	static RotamerRecoveryFactory * get_instance();

	void factory_register( RRProtocolCreatorCOP creator );
	void factory_register( RRComparerCreatorCOP creator );
	void factory_register( RRReporterCreatorCOP creator );
	RRProtocolOP get_rotamer_recovery_protocol( std::string const & type_name );
	RRComparerOP get_rotamer_recovery_comparer( std::string const & type_name );
	RRReporterOP get_rotamer_recovery_reporter( std::string const & type_name );

	RotamerRecoveryOP
	get_rotamer_recovery(
		std::string const & protocol,
		std::string const & comparer,
		std::string const & reporter);

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:

#if defined MULTI_THREADED && defined CXX11
	static std::atomic< RotamerRecoveryFactory * > instance_;
#else
	static RotamerRecoveryFactory * instance_;
#endif

	typedef std::map< std::string, protocols::rotamer_recovery::RRProtocolCreatorCOP > RRProtocolCreatorMap;
	RRProtocolCreatorMap protocol_types_;

	typedef std::map< std::string, protocols::rotamer_recovery::RRComparerCreatorCOP > RRComparerCreatorMap;
	RRComparerCreatorMap comparer_types_;

	typedef std::map< std::string, protocols::rotamer_recovery::RRReporterCreatorCOP > RRReporterCreatorMap;
	RRReporterCreatorMap reporter_types_;

};


/// @brief This templated class will register an instance of an
/// RotamerRecoveryCreator (class T) with the
/// RotamerRecoveryFactory.  It will ensure that no
/// RotamerRecoveryCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class RotamerRecoveryRegistrator : public utility::factory::WidgetRegistrator< RotamerRecoveryFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< RotamerRecoveryFactory, T > parent;
public:
	RotamerRecoveryRegistrator() : parent() {}
};


} // namespace
} // namespace

#endif
