// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLocatorRegistrator.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_managerResourceLocatorRegistrator_hh
#define INCLUDED_basic_resource_managerResourceLocatorRegistrator_hh

//unit headers

//project headers
#include <basic/resource_manager/ResourceLocatorFactory.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers

namespace basic {
namespace resource_manager {

/// @brief The %ResourceLocatorRegistrator gives an instance of a ResourceLocatorCreator to the
/// ResourceLocatorFactory in its constructor, calling the ResourceLocatorFactory's factory_register()
/// method.  This call is actually accomplished by the WidgetRegistrator parent class.  A single
/// (templated) instance of this class for each ResourceLocatorCreator should be placed in the
/// appropriate init.cc file (i.e. ResourceLocatorCreators that live in the protocols library
/// should be put in src/protocols/init/init.ResourceLocatorRegistrators.ihh).
template < class T >
class ResourceLocatorRegistrator : public utility::factory::WidgetRegistrator< ResourceLocatorFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResourceLocatorFactory, T > parent;
public:
	ResourceLocatorRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_managerResourceLocatorRegistrator_hh
