// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/InterchainEnergy.cc
/// @brief  Statistically derived rotamer pair potentials
/// @details For docking (or between chains) only those residues at the interface
///      and between the two interfaces need to be evaluated
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_scoring_InterfaceInfo_fwd_hh
#define INCLUDED_protocols_scoring_InterfaceInfo_fwd_hh


// Unit headers

// Package headers

// Project headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++


namespace protocols {
namespace scoring {

class InterfaceInfo;
typedef utility::pointer::shared_ptr< InterfaceInfo > InterfaceInfoOP;
typedef utility::pointer::shared_ptr< InterfaceInfo const > InterfaceInfoCOP;

} // ns scoring
} // ns core

#endif
