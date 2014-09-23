// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/docking/DockinLowRes.fwd.hh
/// @brief forward declarations for low-res docking
/// @author Jeff Gray

#ifndef INCLUDED_protocols_docking_DockFilters_fwd_hh
#define INCLUDED_protocols_docking_DockFilters_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace docking {

class DockingLowResFilter;
typedef utility::pointer::shared_ptr< DockingLowResFilter >  DockingLowResFilterOP;
typedef utility::pointer::shared_ptr< DockingLowResFilter const >  DockingLowResFilterCOP;

class DockingHighResFilter;
typedef utility::pointer::shared_ptr< DockingHighResFilter >  DockingHighResFilterOP;
typedef utility::pointer::shared_ptr< DockingHighResFilter const >  DockingHighResFilterCOP;


} //docking
} //protocols

#endif
