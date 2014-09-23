// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/hotspot_hashing/PlaceSurfaceProbe.fwd.hh
/// @brief  forward declaration for PlaceSurfaceProbe
/// @author  Alex Ford fordas@uw.edu


#ifndef INCLUDED_protocols_hotspot_hashing_movers_PlaceSurfaceProbe_fwd_hh
#define INCLUDED_protocols_hotspot_hashing_movers_PlaceSurfaceProbe_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

// Forward
class PlaceSurfaceProbe;

// Types
typedef utility::pointer::shared_ptr< PlaceSurfaceProbe >  PlaceSurfaceProbeOP;
typedef utility::pointer::shared_ptr< PlaceSurfaceProbe const >  PlaceSurfaceProbeCOP;

} // namespace movers
} //namespace hotspot_hashing
} // namespace protocols

#endif
