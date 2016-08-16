// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/denovo_design/SheetDB.fwd.hh
/// @brief  SheetDB forward header
/// @author Tom Linsky

#ifndef INCLUDED_protocols_denovo_design_components_SheetDB_fwd_hh
#define INCLUDED_protocols_denovo_design_components_SheetDB_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace denovo_design {
namespace components {

// Forward
class SheetDB;

// Types
typedef utility::pointer::shared_ptr< SheetDB > SheetDBOP;
typedef utility::pointer::shared_ptr< SheetDB const > SheetDBCOP;

typedef utility::pointer::weak_ptr< SheetDB > SheetDBAP;
typedef utility::pointer::weak_ptr< SheetDB const > SheetDBCAP;

} // namespace components
} // namespace denovo_design
} // namespace protocols

#endif
