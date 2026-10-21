// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/modifications/ValenceHandler.hh
/// @brief  A bunch of functions to determine where to place atoms based on hybridization of atoms and number of bonds
///    It should be noted that a vector of coordinates are returned. This is because it creates the angles for all
///    Hydrogens that can be placed. If you are combining a fragment, you really only need the first index of coords
/// @author Rosetta conversion: Steven Combs / Sandeep Kothiwale

#ifndef INCLUDED_core_chemical_modifications_ValenceHandler_hh
#define INCLUDED_core_chemical_modifications_ValenceHandler_hh

#include <core/chemical/MutableResidueType.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>

namespace core {
namespace chemical {
namespace modifications {

//the driver for this file. all the other commands are called by determine coordinates; however, you can use the other functions as much as you want
utility::vector1<numeric::xyzVector<core::Real> > determine_coordinates(core::chemical::MutableResidueType const & res, core::chemical::VD const & atom);
//numeric::xyzVector get_idealized_geometry();

//! point X in B->A->X where A, B and X are on the same line
//! @brief calculates point X in B->A->X where A, B and X are on the same line
//! @param VECTOR_A first point
//! @param VECTOR_B second point
//! @param DISTANCE_XA distance between X and VECTOR_A
//! @return point X in B->A->X where A, B and X are on the same line
numeric::xyzVector<core::Real>  linear_coordinates(numeric::xyzVector<core::Real>  const & vector_a, numeric::xyzVector<core::Real>  const & vector_b, core::Real const distance_xa);


//! @brief calculates coordinates using angle information (point X in B,C->A->X)
//! @param VECTOR_A first point
//! @param VECTOR_B second point
//! @param VECTOR_C third point
//! @param DISTANCE_XA distance between X and VECTOR_A
//! @param ANGLE_XAB angle between X, VECTOR_A, and VECTOR_B
//! @param ANGLE_XAC angle between X, VECTOR_A, and VECTOR_C
//! @param SIDE true if on a side
//! @param VECTOR_SIDE vector of the side
//! @return point X in B,C->A->X
numeric::xyzVector<core::Real>  angle_coordinates(
	const numeric::xyzVector<core::Real>  &VECTOR_A,
	const numeric::xyzVector<core::Real>  &VECTOR_B,
	const numeric::xyzVector<core::Real>  &VECTOR_C,
	const core::Real DISTANCE_XA,
	const core::Real ANGLE_XAB,
	const core::Real ANGLE_XAC,
	const bool SIDE,
	const numeric::xyzVector<core::Real> &VECTOR_SIDE
);

//! @brief calculates the unit vector starting from one linal::Vector3D to another
//! @param ORIGIN vector of origin
//! @param TARGET target vector
//! @return the unit vector between ORIGIN and TARGET
numeric::xyzVector<core::Real>  triganol_coordinates
(
	const numeric::xyzVector<core::Real>  &VECTOR_A,
	const numeric::xyzVector<core::Real>  &VECTOR_B,
	const numeric::xyzVector<core::Real>  &VECTOR_C,
	const double DISTANCE_XA
);


//! point X in B,C,D->A->X
//! @brief calculates point X in B,C,D->A->X
//! @param VECTOR_A first point
//! @param VECTOR_B second point
//! @param VECTOR_C third point
//! @param VECTOR_D fourth point
//! @param DISTANCE_XA distance between X and VECTOR_A
//! @return point X in B,C,D->A->X
numeric::xyzVector<core::Real>  tetrahedral_coordinates
(
	const numeric::xyzVector<core::Real>  &VECTOR_A,
	const numeric::xyzVector<core::Real>  &VECTOR_B,
	const numeric::xyzVector<core::Real>  &VECTOR_C,
	const numeric::xyzVector<core::Real>  &VECTOR_D,
	const double DISTANCE_XA
);


//! @brief calculates projection angle between two linal::Vector3D
//! @param VECTOR_A first vector (point)
//! @param VECTOR_B second vector (point)
//! @return projection angle between two linal::Vector3D
double projection_angle( const numeric::xyzVector<core::Real>  &VECTOR_A, const numeric::xyzVector<core::Real>  &VECTOR_B);

double projection_angle_cosin( const numeric::xyzVector<core::Real>  &VECTOR_A, const numeric::xyzVector<core::Real>  &VECTOR_B);

//! @brief dihedral angle between four points (A->B -x-> C->D)
//! @brief see http://en.wikipedia.org/wiki/Dihedral_angle for reference
//! @param VECTOR_A first vector (point)
//! @param VECTOR_B second vector (point)
//! @param VECTOR_C third vector (point)
//! @param VECTOR_D fourth vector (point)
//! @return dihedral angle between four points
double dihedral_coordinates
(
	const numeric::xyzVector<core::Real> &VECTOR_A,
	const numeric::xyzVector<core::Real> &VECTOR_B,
	const numeric::xyzVector<core::Real> &VECTOR_C,
	const numeric::xyzVector<core::Real> &VECTOR_D
);


}
}
}

#endif /* VALENCEHANDLER_HH_ */
