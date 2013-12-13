// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_recovery/RRComparerElecDensDiff.fwd.hh
/// @brief Measure rotamer recovery by change in sidechain electron density correlation 
/// @author Patrick Conway (ptconway@uw.edu)

#ifndef INCLUDED_protocols_rotamer_recovery_RRComparerElecDensDiff_fwd_hh
#define INCLUDED_protocols_rotamer_recovery_RRComparerElecDensDiff_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace rotamer_recovery{

class RRComparerElecDensDiff;
typedef utility::pointer::owning_ptr< RRComparerElecDensDiff > RRComparerElecDensDiffOP;
typedef utility::pointer::owning_ptr< RRComparerElecDensDiff const > RRComparerElecDensDiffCOP;


}//rotamer_recovery
}//protocols

#endif //INCLUDED_protocols_rotamer_recovery_RRComparerElecDensDiff_FWD_HH
