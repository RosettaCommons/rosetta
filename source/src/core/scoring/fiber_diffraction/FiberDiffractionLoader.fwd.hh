// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionLoader.fwd.hh
/// @brief  Fiber diffraction data forward creator
/// @author Wojciech Potrzebowski and Ingemar Andre

#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionLoader_fwd_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionLoader_fwd_hh

//Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace fiber_diffraction {

class FiberDiffractionLoader;
typedef utility::pointer::shared_ptr< FiberDiffractionLoader > FiberDiffractionLoaderOP;
typedef utility::pointer::shared_ptr< FiberDiffractionLoader const > FiberDiffractionLoaderCOP;

} // namespace fiber_diffraction
} //namespace scoring
} //namespace core

#endif
