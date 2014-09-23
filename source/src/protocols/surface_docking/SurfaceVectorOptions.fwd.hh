// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/surface_docking/SurfaceVectorOptions.fwd.hh
/// @brief  SurfaceVectorOptions forward declarations header
/// @author Michael Pacella mpacella88@gmail.com

#ifndef INCLUDED_protocols_surface_docking_SurfaceVectorOptions_FWD_HH
#define INCLUDED_protocols_surface_docking_SurfaceVectorOptions_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace surface_docking {

// Forward
class SurfaceVectorOptions;

typedef utility::pointer::shared_ptr< SurfaceVectorOptions > SurfaceVectorOptionsOP;
typedef utility::pointer::shared_ptr< SurfaceVectorOptions const > SurfaceVectorOptionsCOP;

} //namespace surface_docking
} //namespace protocols

#endif //INCLUDED_protocols_SurfaceVectorOptions_FWD_HH
