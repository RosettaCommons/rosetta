// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/icoor_support.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_icoor_support_hh
#define INCLUDED_core_chemical_icoor_support_hh

// Project headers
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/AtomICoor.fwd.hh>

#include <core/kinematics/tree/Atom.fwd.hh>

#include <map>

namespace core {
namespace chemical {

/// @brief Attempt to find new ICOOR definitions for entries in the ResidueType which rely on a now-deleted upper/lower connection
/// (Connection type is specified by the passed type.)
void clean_up_dangling_connect( core::chemical::MutableResidueType & restype, ICoordAtomIDType remove_type );

//The core::kinematics::tree::Atom objects take care of the internal/xyz exchange.
//We just need to chain them together appropriately.
typedef std::map< core::chemical::VD, core::kinematics::tree::AtomOP > VdTreeatomMap;

class RerootRestypeVisitor;
class RerootEdgeSorter;

void
reroot_restype(
	core::chemical::MutableResidueType & restype,
	core::chemical::ResidueGraph const & graph,
	core::chemical::VD root);

void
fill_ideal_xyz_from_icoor(
	core::chemical::MutableResidueType & restype,
	core::chemical::ResidueGraph const & graph);

} // chemical
} // core

#endif
