// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/hotspot_hashing/LSMGridProbe.fwd.hh
/// @brief  forward declaration for LSMGridProbe
/// @author  Alex Ford fordas@uw.edu


#ifndef INCLUDED_protocols_hotspot_hashing_filters_LSMGridProbe_fwd_hh
#define INCLUDED_protocols_hotspot_hashing_filters_LSMGridProbe_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols
{
namespace hotspot_hashing
{
namespace filters
{

// Forward
class LSMGridProbe;

// Types
typedef utility::pointer::owning_ptr< LSMGridProbe >  LSMGridProbeOP;
typedef utility::pointer::owning_ptr< LSMGridProbe const >  LSMGridProbeCOP;

} // namespace filters
} //namespace hotspot_hashing
} // namespace protocols

#endif
