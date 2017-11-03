// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/bogus.cc
/// @author Andrew Leaver-Fay
/// @details Fix the static-initialization-order fiasco w/ global variable initialization
/// Previously, these globals were in separate compilation units, and the order of their
/// initialization was not guaranteed. Putting them into a single compilation unit fixes that.

// Unit headers
#include <core/id/bogus.hh>

// Package headers
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>

namespace core {
namespace id {

AtomID const GLOBAL_BOGUS_ATOM_ID( AtomID::BOGUS_ATOM_ID() );

AtomID const GLOBAL_CHAINBREAK_BOGUS_ATOM_ID( AtomID::CHAINBREAK_BOGUS_ATOM_ID() );
StubID const GLOBAL_BOGUS_STUB_ID( StubID::BOGUS_STUB_ID() );
DOF_ID const GLOBAL_BOGUS_DOF_ID( AtomID::BOGUS_ATOM_ID(), PHI );

NamedAtomID const GLOBAL_BOGUS_NAMED_ATOM_ID( NamedAtomID::BOGUS_NAMED_ATOM_ID() );
NamedAtomID const GLOBAL_CHAINBREAK_BOGUS_NAMED_ATOM_ID( NamedAtomID::CHAINBREAK_BOGUS_NAMED_ATOM_ID() );

/// @brief Globals
TorsionID const GLOBAL_BOGUS_TORSION_ID( TorsionID::BOGUS_TORSION_ID() );

void
initialize_core_id_globals()
{
	// Calling this function ensures that all of the global variables in core::id
}

} // namespace id
} // namespace core
