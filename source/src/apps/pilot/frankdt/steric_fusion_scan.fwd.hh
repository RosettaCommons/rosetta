// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/steric_fusion_scan.fwd.hh
/// @brief searches two proteins for the most sterically hindering fusion
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_apps_pilot_frankdt_steric_fusion_scan_fwd_hh
#define INCLUDED_apps_pilot_frankdt_steric_fusion_scan_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace apps {
namespace pilot {
namespace frankdt {

class steric_fusion_scan;

typedef utility::pointer::shared_ptr< steric_fusion_scan > steric_fusion_scanOP;
typedef utility::pointer::shared_ptr< steric_fusion_scan const > steric_fusion_scanCOP;

} //apps
} //pilot
} //frankdt

#endif //INCLUDED_apps_pilot_frankdt_steric_fusion_scan_fwd_hh
