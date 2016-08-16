// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/AntibodyInfoRMFallbackConfiguration.fwd.hh
/// @brief  forward header for OptionsSystemFallback class
/// @author Michael Pacella mpacella88@gmail.com

#ifndef INCLUDED_protocols_antibody_AntibodyInfoRMfallbackconfiguration_FWD_HH
#define INCLUDED_protocols_antibody_AntibodyInfoRMfallbackconfiguration_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {

class AntibodyInfoRMFallbackConfiguration;
typedef utility::pointer::shared_ptr< AntibodyInfoRMFallbackConfiguration > AntibodyInfoRMFallbackConfigurationOP;
typedef utility::pointer::shared_ptr< AntibodyInfoRMFallbackConfiguration const > AntibodyInfoRMFallbackConfigurationCOP;

} // antibody
} // protocols

#endif // INCLUDED_protocols_antibody__fallback_configuration_FWD_HH
