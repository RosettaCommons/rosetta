// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceLoaderRegistrator.hh
/// @brief  Declaration of class ResourceLoaderRegistrator
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
///         Brian Weitzner (brian.weitzner@gmail.com)
///         Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_resource_manager_ResourceLoaderRegistrator_hh
#define INCLUDED_basic_resource_manager_ResourceLoaderRegistrator_hh

//unit headers

//project headers
#include <basic/resource_manager/ResourceLoaderFactory.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers

namespace basic {
namespace resource_manager {

/// @brief The %ResourceLoaderRegistrator class is a simple templated registration class
/// that will, at construction, create a ResourceLoader and register it with the
/// ResouceLoaderFactory
///
/// Instances of the %ResourceLoaderRegistrator class should be added to the "init.cc" files
/// of the libraries they belong to so that they will be constructed along the init(argv,argc)
/// calling pathway at the very beginning of program execution.
template < class T >
class ResourceLoaderRegistrator : public utility::factory::WidgetRegistrator< ResourceLoaderFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResourceLoaderFactory, T > parent;
public:
	ResourceLoaderRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_manager_ResourceLoaderRegistrator_hh
