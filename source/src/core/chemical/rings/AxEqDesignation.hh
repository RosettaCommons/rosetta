// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/AxEqDesignation.hh
/// @brief   Enumeration definition for AxEqDesignation.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_rings_AxEqDesignation_HH
#define INCLUDED_core_chemical_rings_AxEqDesignation_HH

namespace core {
namespace chemical {
namespace rings {

/// @brief:  The axial/equatorial designation for bonds/substituents on ring systems.
/// @details Axial: "...[B]onds to ring atoms (and molecular entities attached to such bonds) are... axial... [if] the
/// bonds make a relatively large... angle... with the plane containing or passing closest to a majority of the ring
/// atoms.  ...[A]xial bonds are approximately parallel to the C3 axis....\n
/// Equatorial: "...[B]onds to ring atoms (and molecular entities attached to such bonds) are... equatorial... [if] the
/// bonds make a relatively small... angle... with the plane containing or passing closest to a majority of the ring
/// atoms.  ...[E]quatorial bonds [are] approximately parallel to two of the ring bonds.
/// @ref     http://goldbook.iupac.org/A00546.html
/// @note    NEITHER has the value of 0 and so can be used in conditionals as false.
enum AxEqDesignation {
	NEITHER = 0,
	AXIAL,
	EQUATORIAL
};

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_AxEqDesignation_HH
