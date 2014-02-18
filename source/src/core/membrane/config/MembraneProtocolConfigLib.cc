// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MembraneProtocolConfigLib.cc
///
/// @brief Stores membrane specific protocol changes
/// @detail Stores option changes, weights patches, and apply changes for specific protocols
///         that are now membrane compatible
///
/// @author Rebecca Alford

#ifndef INCLUDED_core_membrane_config_MembraneProtocolConfigLib_cc
#define INCLUDED_core_membrane_config_MembraneProtocolConfigLib_cc

// Unit Headers
#include <core/membrane/config/MembraneProtocolConfigLib.hh>

// Package Headers
#include <basic/Tracer.hh>

// Project Headers

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// C++ Headers

static basic::Tracer TR( "core.membrane.config.MembraneProtocolConfigLib" );

namespace core {
namespace membrane {
namespace config {

/// @brief Constructors
MembraneProtocolConfigLib::MembraneProtocolConfigLib() : utility::pointer::ReferenceCount()
{}

/// @brief Destructor
MembraneProtocolConfigLib::~MembraneProtocolConfigLib()
{}

} // config
} // membrane
} // core

#endif INCLUDED_core_membrane_config_MembraneProtocolConfigLib_cc