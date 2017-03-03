// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kinematic_closure/bridgeObjects_nonredundant.hh
/// @brief Non-redundant bridgeObject function.
/// @detail The origional bridgeObject function need some redundant inputs which are just
/// place holders and never used in actual calculation. Here I wrote a wrapper of the original
/// function such that no redundant information are needed. However, this wrapper means even
/// more moving around redundant data. TODO: rewrite the function such that redundant data don't
/// appear at all.
/// @author xingjiepan (xingjiepan@gmail.com)

#ifndef INCLUDED_numeric_kinematic_closure_bridgeObjects_nonredundant_hh
#define INCLUDED_numeric_kinematic_closure_bridgeObjects_nonredundant_hh

#include <numeric/types.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/fixedsizearray1.fwd.hh>

namespace numeric {
namespace kinematic_closure {


/// @brief Nonredundant version of the bridgeObject function
/// @detail stub1 are the coordinates of the first pivot and 2 atoms preceeding it;
/// stub2 are the coordinates of the third pivot and 2 atoms after it;
/// torsions_chain1 are the torsions from pivot1 to pivot2, which has a length of (pivot2 - pivot1 - 2);
/// torsions_chain2 are the torsions from pivot2 to pivot3, which has a length of (pivot3 - pivot2 - 2);
/// angles are the bond angles from pivot1 to pivot3, which has a length of (pivot3 - pivot1 + 1);
/// bonds are the bond lengths from pivot1 to pivot3, which has a length of (pivot3 - pivot1);
/// pivot_torsions are the solutions of pivot tosions whose dimension is nsol * 6 where nsol is the number
/// of solutions.
void bridgeObjects_nonredundant(const utility::vector1<utility::fixedsizearray1<Real,3> >& stub1,
	const utility::vector1<utility::fixedsizearray1<Real,3> >& stub2,
	const utility::vector1<numeric::Real> & torsions_chain1, const utility::vector1<numeric::Real> & torsions_chain2,
	const utility::vector1<numeric::Real> & angles, const utility::vector1<numeric::Real> & bonds,
	utility::vector1<utility::vector1<Real> >& pivot_torsions,
	int &nsol);



} //numeric
} //kinematic_closure


#endif //numeric/kinematic_closure_bridgeObjects_nonredundant_hh

