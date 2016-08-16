// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rotamer_recovery/RRProtocol.fwd.hh
/// @brief  Manage how the rotamer recovery test should be preformed
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_rotamer_recovery_RRProtocol_fwd_hh
#define INCLUDED_protocols_rotamer_recovery_RRProtocol_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rotamer_recovery {

class RRProtocol;
typedef utility::pointer::shared_ptr< RRProtocol > RRProtocolOP;
typedef utility::pointer::shared_ptr< RRProtocol const > RRProtocolCOP;

} // namespace
} // namespace

#endif //include guard
