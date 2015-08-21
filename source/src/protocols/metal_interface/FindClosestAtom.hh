// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/metal_interface/FindClosestAtom.hh
/// @brief Finds the closest atom in a given residue to a point (usually a zinc atom).
/// @author Bryan Der

#ifndef INCLUDED_protocols_metal_interface_FindClosestAtom_HH
#define INCLUDED_protocols_metal_interface_FindClosestAtom_HH

#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <string>

namespace protocols {
namespace metal_interface {

//helper function - iterates over all sidechain non-carbon heavy atoms of res to find the one closest to xyz
std::string find_closest_atom( core::conformation::Residue const & res, core::Vector const & xyz );

}//namespace metal_interface
}//namespace protocols

#endif // INCLUDED_protocols_metal_interface_FindClosestAtom_HH
