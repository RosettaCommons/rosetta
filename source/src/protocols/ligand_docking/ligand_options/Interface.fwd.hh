// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/protocols/ligand_docking/ligand_options/Interface.fwd.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_ligand_options_Interface_fwd_hh
#define INCLUDED_protocols_ligand_docking_ligand_options_Interface_fwd_hh

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {
namespace ligand_options {

/// @brief info for each residue- is it part of the interface and if so, what ligands is it near
struct InterfaceInfo;

/// @brief For each residue is it in the interface, a mobile region or a non-mobile region?
class Interface;

} //namespace ligand_options
} //namespace ligand_docking
} //namespace protocols

#endif
