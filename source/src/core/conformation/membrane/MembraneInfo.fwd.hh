// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/conformation/membrane/MembraneInfo.fwd.hh
///
/// @brief		Information about the membrane bilayer and its relationship with the protein(s)
/// @details	MembraneInfo is responsible for describing attributes of the 
/// 			membrane bilayer *and* the position & orientation of the bilayer
///				in 3D space. All of the data members work together to accomplish 
///				this representation. And the players are: 
///					= A pointer to the Pose Conformation, which includes an MEM residue
///					= The lipid type represented in this membrane (Default: DOPC: (18:1/18:1))
///					= Profiles describing chemical characteristics of the bilayer. Including: 
///						= Electrostatics
///						= Hydrogen bonding potential
///						= Polarity
///					= Topology of transmembrane spans
///					= Per-residue lipophilicity
///					= Membrane thickness & steepness - derived from chemical profiles
///
///				This object is a member of Conformation and should be accessed by
///				pose.conofrmation().membrane_info(). DO NOT access the MEM residue
///				outside of the framework!!! 
///
///				Last Updated: 7/23/15
/// @author		Rebecca Faye Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {
		
/// @brief MembraneInfo describes the membrane bilayer and its relationship with the protein
class MembraneInfo;
typedef utility::pointer::shared_ptr< MembraneInfo > MembraneInfoOP;
typedef utility::pointer::shared_ptr< MembraneInfo const > MembraneInfoCOP;
		
} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembraneInfo_fwd_hh

