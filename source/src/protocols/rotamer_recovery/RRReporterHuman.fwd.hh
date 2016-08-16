// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rotamer_recovery/RRReporterHuman.fwd.hh
/// @brief Report how well rosetta is a recovering the conformation of
/// an experimentally validated structure in a human readable format.
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_rotamer_recovery_RRReporterHuman_FWD_HH
#define INCLUDED_protocols_rotamer_recovery_RRReporterHuman_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rotamer_recovery {

class PerNativeRRReporterHuman;
typedef utility::pointer::shared_ptr< PerNativeRRReporterHuman > PerNativeRRReporterHumanOP;
typedef utility::pointer::shared_ptr< PerNativeRRReporterHuman const > PerNativeRRReporterHumanCOP;

class RRReporterHuman;
typedef utility::pointer::shared_ptr< RRReporterHuman > RRReporterHumanOP;
typedef utility::pointer::shared_ptr< RRReporterHuman const > RRReporterHumanCOP;


}//rotamer_recovery
}//protocols

#endif //INCLUDED_protocols_rotamer_recovery_RRReporterHuman_FWD_HH
