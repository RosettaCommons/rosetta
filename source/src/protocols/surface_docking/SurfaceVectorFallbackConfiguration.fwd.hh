// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/surface_docking/SurfaceVectorFallbackConfiguration.fwd.hh
/// @brief  forward header for OptionsSystemFallback class
/// @author Michael Pacella mpacella88@gmail.com

#ifndef INCLUDED_protocols_surface_docking_surface_vector_fallback_configuration_FWD_HH
#define INCLUDED_protocols_surface_docking_surface_vector_fallback_configuration_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace surface_docking {

class SurfaceVectorFallbackConfiguration;
typedef utility::pointer::shared_ptr< SurfaceVectorFallbackConfiguration > SurfaceVectorFallbackConfigurationOP;
typedef utility::pointer::shared_ptr< SurfaceVectorFallbackConfiguration const > SurfaceVectorFallbackConfigurationCOP;

} // surface_docking
} // protocols

#endif // INCLUDED_protocols_surface_docking_surface_vector_fallback_configuration_FWD_HH
