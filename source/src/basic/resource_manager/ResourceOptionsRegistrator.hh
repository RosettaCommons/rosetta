// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceLoader.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceOptionsRegistrator_hh
#define INCLUDED_basic_resource_manager_ResourceOptionsRegistrator_hh

//unit headers

//package headers
#include <basic/resource_manager/ResourceOptionsFactory.hh>

//utility headers
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers

namespace basic {
namespace resource_manager {

/// @brief The %ResourceOptionsRegistrator class is responsible for creating an instance of the (templated)
/// ResourceOptionsCreator class and giving it to the ResourceOptionsFactory in its construtor.  Instances
/// of this class placed in the init.cc files (e.g. core/init/init.cc) ensure that the ResourceOptionsFactory
/// is fully populated with the ResourceOptionsCreators by the time that the call to devel::init() completes.
template < class T >
class ResourceOptionsRegistrator : public utility::factory::WidgetRegistrator< ResourceOptionsFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResourceOptionsFactory, T > parent;
public:
	ResourceOptionsRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceOptionsRegistrator_hh
