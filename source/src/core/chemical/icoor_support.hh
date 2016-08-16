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
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>

#include <core/kinematics/tree/Atom.fwd.hh>

#include <map>

namespace core {
namespace chemical {

//The core::kinematics::tree::Atom objects take care of the internal/xyz exchange.
//We just need to chain them together appropriately.
typedef std::map< core::chemical::VD, core::kinematics::tree::AtomOP > VdTreeatomMap;

class RerootRestypeVisitor;
class RerootEdgeSorter;

void
reroot_restype(
	core::chemical::ResidueType & restype,
	core::chemical::ResidueGraph const & graph,
	core::chemical::VD root);

void
fill_ideal_xyz_from_icoor(
	core::chemical::ResidueType & restype,
	core::chemical::ResidueGraph const & graph);

} // chemical
} // core

#endif
