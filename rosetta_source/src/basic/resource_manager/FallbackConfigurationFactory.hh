// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfigurationFactory.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_factory_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_factory_HH

//unit headers
#include <basic/resource_manager/FallbackConfigurationFactory.fwd.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfigurationCreator.fwd.hh>

// Utility headers
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers
#include <map>
#include <string>

namespace basic {
namespace resource_manager {

class FallbackConfigurationFactory {
private:
	FallbackConfigurationFactory(); // singleton, private constructor

public:

	FallbackConfigurationOP
	create_fallback_configuration( std::string const & resource_description	) const;

	static FallbackConfigurationFactory * get_instance();

	void
	factory_register( FallbackConfigurationCreatorOP creator );
	
	bool
	has_fallback_for_resource( std::string const & desc ) const;

private:
	static FallbackConfigurationFactory * instance_;

	typedef std::map< std::string, FallbackConfigurationCreatorOP > FallbackConfigurationCreatorsMap;
	FallbackConfigurationCreatorsMap creators_map_;

};

template < class T >
class FallbackConfigurationRegistrator : public utility::factory::WidgetRegistrator< FallbackConfigurationFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< FallbackConfigurationFactory, T > parent;
public:
	FallbackConfigurationRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_fallback_configuration_factory_HH
