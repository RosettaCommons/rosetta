// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceManagerCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceManagerCreator_HH
#define INCLUDED_basic_resource_manager_ResourceManagerCreator_HH

//unit headers
#include <basic/resource_manager/ResourceManagerCreator.fwd.hh>

// package headers
#include <basic/resource_manager/ResourceManager.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <string>

namespace basic {
namespace resource_manager {

/// @brief Derived classes will be used by the ResourceManagerFactory to decide which
/// of the various ResourceManagers should be instantiated.  The ResourceManager
/// is a singleton, but, different ResourceManagers can be instantiated in different
/// contexts.
class ResourceManagerCreator : public utility::pointer::ReferenceCount {
public:
	friend class ResourceManagerFactory;


	~ResourceManagerCreator() override;

private:
	/// @brief Returns a raw pointer, not an owning pointer, to the ResourceManager.
	/// The ResourceManager is a singleton class, and will be held in a raw pointer by the
	/// base class.  Each derived ResourceManagerCreator should be a friend of the
	/// corresponding ResourceManager.  This function is private, meaning it should
	/// only be called by its friend class, the ResourceManagerFactory.
	virtual
	ResourceManager *
	create_resource_manager() const = 0;

	/// @brief Give the name of the ResourceManager that this Creator instantiates
	/// to the ResourceManagerFactory.  The ResourceManagerFactory reads the option
	/// system and instantiates the appropriate ResourceManager by string lookup:
	/// this is a hokey way of "decoupling" the logic for ResourceManager instantiation
	/// that must happen high up in the library hierarchy: the ResourceManagerFactory
	/// does not, e.g., need to #include the JD2ResourceManagerCreator from protocols::jd2.
	virtual
	std::string manager_type() const = 0;

};


} // namespace resource_manager
} // namespace basic

#endif //INCLUDED_basic_resource_manager_ResourceManager_HH
