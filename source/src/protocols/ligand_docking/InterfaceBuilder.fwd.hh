// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/ligand_docking/InterfaceBuilder.fwd.hh
/// @brief  Makes an interface
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_InterfaceBuilder_fwd_hh
#define INCLUDED_protocols_ligand_docking_InterfaceBuilder_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_docking {

class InterfaceBuilder;

typedef utility::pointer::shared_ptr< InterfaceBuilder > InterfaceBuilderOP;
typedef utility::pointer::shared_ptr< InterfaceBuilder const > InterfaceBuilderCOP;
typedef utility::vector1< InterfaceBuilderOP > InterfaceBuilderOPs;
typedef utility::vector1< InterfaceBuilderCOP > InterfaceBuilderCOPs;

} //namespace ligand_docking
} //namespace protocols

#endif
