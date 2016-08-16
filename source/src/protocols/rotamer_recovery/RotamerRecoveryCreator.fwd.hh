// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/RotamerRecoveryCreator.hh
/// @brief  Forward Header for base class for Rotamer Recovery Creators for the RotamerRecovery load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_rotamer_recovery_RotamerRecoveryCreator_fwd_hh
#define INCLUDED_protocols_rotamer_recovery_RotamerRecoveryCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rotamer_recovery {

class RRProtocolCreator;
typedef utility::pointer::shared_ptr< RRProtocolCreator > RRProtocolCreatorOP;
typedef utility::pointer::shared_ptr< RRProtocolCreator const > RRProtocolCreatorCOP;

class RRComparerCreator;
typedef utility::pointer::shared_ptr< RRComparerCreator > RRComparerCreatorOP;
typedef utility::pointer::shared_ptr< RRComparerCreator const > RRComparerCreatorCOP;

class RRReporterCreator;
typedef utility::pointer::shared_ptr< RRReporterCreator > RRReporterCreatorOP;
typedef utility::pointer::shared_ptr< RRReporterCreator const > RRReporterCreatorCOP;

} //namespace
} //namespace

#endif
