// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/solid_surface/SurfaceEnergies.fwd.hh
/// @brief  SurfaceEnergies forward declarations header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Mike Pacella (mpacella@gmail.com)


#ifndef INCLUDED_core_scoring_solid_surface_SurfaceEnergies_fwd_hh
#define INCLUDED_core_scoring_solid_surface_SurfaceEnergies_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace solid_surface {

class SurfaceEnergies;
typedef utility::pointer::shared_ptr< SurfaceEnergies       > SurfaceEnergiesOP;
typedef utility::pointer::shared_ptr< SurfaceEnergies const > SurfaceEnergiesCOP;

} // namespace solid_surface
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_solid_surface_SurfaceEnergies_FWD_HH
