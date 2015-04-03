// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Bond.fwd.hh
/// @author Rocco Moretti (rmorettase@gmail.com), Gordon Lemmon, Phil Bradley


#ifndef INCLUDED_core_chemical_Bond_fwd_hh
#define INCLUDED_core_chemical_Bond_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace chemical {

enum BondName
{
	UnknownBond=0,
	SingleBond=1,
	DoubleBond=2,
	TripleBond=3,
	AromaticBond=4,
	OrbitalBond=5
};

enum BondOrder
{
	UnknownBondOrder=0,
	SingleBondOrder=1,
	DoubleBondOrder=2,
	TripleBondOrder=3,
	OrbitalBondOrder,
	PseudoBondOrder=99 // Not a true bond - just there for connectivity representation.
};

/// @brief As with the BCL, bond conjugability is more about the atom types on
/// either end of the bond than about the bond itself. A conjugatable bond is
/// one where it is known that the atoms on *both* sides of the bond can participate
/// in a conjugated system. Double, triple and aromatic bonds are always conjugatable,
/// and single bonds are conjugatable if *both* atoms are in other triple, double,
/// or aromatic systems

enum BondConjugability
{
	UnknownConjugability,
	NotConjugableBond,
	ConjugableBond
};

enum BondRingness
{
	UnknownRingness,
	BondNotInRing,
	BondInRing
};

/// @brief Proper aromaticity implies participation in a 4n+2 electron ring system.
/// For simple single-double alternation, see conjugatable bond.
enum BondAromaticity
{
	UnknownAromaticity,
	NonaromaticBond,
	IsAromaticBond
};

// @brief Bond isometry for double bonds.
enum BondIsometry
{
	UnknownIsometry,
	NoBondIsometry,
	EIsometry,
	ZIsometry
};

class Bond;

//typedef  utility::pointer::owning_ptr< Bond >  BondOP;
//typedef  utility::pointer::owning_ptr< Bond const >  BondCOP;
//typedef  utility::vector1< BondOP >  BondOPs;
//typedef  utility::vector1< BondCOP >  BondCOPs;

} // chemical
} // core


#endif // INCLUDED_core_chemical_Bond_FWD_HH
