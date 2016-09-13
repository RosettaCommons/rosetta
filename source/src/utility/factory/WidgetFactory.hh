// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/factory/WidgetFactory.hh
/// @brief  WidgetFactory base class for load-time registration of WidgetCreators to the WidgetFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_factory_WidgetFactory_hh
#define INCLUDED_utility_factory_WidgetFactory_hh


// Unit headers
#include <utility/factory/WidgetFactory.fwd.hh>
#include <utility/exit.hh>

// Package headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <utility/assert.hh>
#include <map>
#include <string>


namespace utility {
namespace factory {


/// @brief Factory base class holds a map between strings and owning pointers of
/// the creator classes.  This should only be used as a base class for a WidgetFactory
/// which expects to map between strings and Creators.  WidgetRegistrators
/// may be used reguardless of how the Factory method maps to its Creators.
template< class Creator >
class WidgetFactory
{
public:
	typedef typename utility::pointer::shared_ptr< Creator > CreatorOP;
	typedef typename std::map< std::string, CreatorOP >      CreatorMap;

public:

	WidgetFactory() = default;
	virtual ~WidgetFactory() {}

	void factory_register( CreatorOP creator ) {
		/// This would be an appropriate place to lock a mutex.
		typename CreatorMap::const_iterator iter = creators_.find( creator->widget_name() );
		if ( iter != creators_.end() ) {
			utility_exit_with_message( "factory_register (" + factory_name() +
				") failed as Widget with the name" + creator->widget_name() + " has already been registered" );
		}
		creators_[ creator->widget_name() ] = creator;
		/// This would be an appropriate place to unlock a mutex.
	}

	virtual std::string factory_name() const = 0;

protected:
	CreatorMap const & creators() const {
		return creators_;
	}

private:

	CreatorMap creators_;
};

} // namespace factory
} // namespace utility


#endif // INCLUDED_utility_factory_WidgetFactory_HH
