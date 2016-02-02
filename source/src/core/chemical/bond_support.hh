// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/bond_support.hh
/// @brief support functions for class Bond; functions that
/// should not be included as part of the class.
/// @author Steven Combs


#ifndef INCLUDED_core_chemical_bond_support_hh
#define INCLUDED_core_chemical_bond_support_hh

#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh>

namespace core {
namespace chemical {

/// @brief Determine which bonds are in rings, and set the BondRingness property of each
void find_bonds_in_rings(ResidueType & res, bool const complex_ring_detection = false);
void complex_ring_detection( ResidueType & res);
void quick_ring_detection( ResidueType & res);
    
utility::vector1<VD> get_connecting_atoms(ResidueType const & res, ED const & edge);
utility::vector1<VD> get_connecting_atoms(ResidueGraph const & res, ED const & edge);
ED get_bond(ResidueType const & res, VD const & source, VD const & target);

//this will create a bond length based on gasteiger atom type definitions of bonds
Real create_bond_length(gasteiger::GasteigerAtomTypeData const & atom1,
	gasteiger::GasteigerAtomTypeData const & atom2, BondName bond_type);

/// @brief Find which bonds are rotatatable (chi) bonds
/// Returns a list of four vds representing the chi
utility::vector1<VDs> find_chi_bonds( ResidueType const & restype );

/// @brief Is the given chi a proton chi with the proton attached to an atom attached to an non-sp3 atom?
/// @details The use case is to see if the proton chi should flat or staggered with rotamers
bool is_sp2_proton_chi( core::Size chi, ResidueType const & restype );

}
}


#endif

