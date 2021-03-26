// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaOutputsBase.cc
/// @brief A pure virtual base class for the outputs of trRosetta.  Derived classes are for particular trRosetta versions (to allow for future versions providing additional outputs).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/trRosetta/trRosettaOutputsBase.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.trRosetta.trRosettaOutputsBase" );


namespace protocols {
namespace trRosetta {

/// @brief Default constructor.
trRosettaOutputsBase::trRosettaOutputsBase() = default;

/// @brief Destructor.
trRosettaOutputsBase::~trRosettaOutputsBase() = default;

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
trRosettaOutputsBaseOP
trRosettaOutputsBase::clone() const {
	return utility::pointer::make_shared< trRosettaOutputsBase >( *this );
}

} //trRosetta
} //protocols
