// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/hotspot_hashing/PlaceProbe.fwd.hh
/// @brief  forward declaration for PlaceProbe
/// @author  Sergey Lyskov


#ifndef INCLUDED_protocols_hotspot_hashing_movers_PlaceProbe_fwd_hh
#define INCLUDED_protocols_hotspot_hashing_movers_PlaceProbe_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols
{
namespace hotspot_hashing
{
namespace movers
{

// Forward
class PlaceProbeMover;

// Types
typedef utility::pointer::shared_ptr< PlaceProbeMover >  PlaceProbeMoverOP;
typedef utility::pointer::shared_ptr< PlaceProbeMover const >  PlaceProbeMoverCOP;

} // namespace movers
} //namespace hotspot_hashing
} // namespace protocols

#endif
