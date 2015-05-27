// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/AntibodyInfoRMFallbackConfigurationCreator.hh
/// @brief  
/// @author Michael Pacella mpacella88@gmail.com

#ifndef INCLUDED_protocols_antibody_AntibodyInfofallbackconfiguration_creator_HH
#define INCLUDED_protocols_antibody_AntibodyInfofallbackconfiguration_creator_HH

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

//C++ headers
#include <string>

namespace protocols {
namespace antibody {

class AntibodyInfoRMFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator
{
public:
	virtual
	basic::resource_manager::FallbackConfigurationOP
	create_fallback_configuration() const;

	virtual
	std::string resource_description() const;

};


} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_AntibodyInfoRM_fallback_configuration_creator_HH
