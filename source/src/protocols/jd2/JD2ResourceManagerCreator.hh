// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManagerCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_jd2_JD2ResourceManagerCreator_HH
#define INCLUDED_protocols_jd2_JD2ResourceManagerCreator_HH

// Unit Headers
#include <basic/resource_manager/ResourceManagerCreator.hh>

// Package headers
#include <basic/resource_manager/ResourceManager.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>
#include <string>

namespace protocols {
namespace jd2 {

/// @brief The %JD2ResourceManagerCreator is responsible for instantiating the JD2ResourceManager
/// for the ResourceManagerFactory
class JD2ResourceManagerCreator : public basic::resource_manager::ResourceManagerCreator
{
public:
	
	~JD2ResourceManagerCreator() override;

	
	basic::resource_manager::ResourceManager *
	create_resource_manager() const override;

	
	std::string manager_type() const override;


};

} // namespace resource_manager
} // namespace basic

#endif
