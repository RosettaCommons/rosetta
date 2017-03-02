// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/kinematic_closure/bridgeObjects_nonredundant.cc
/// @brief Non-redundant bridgeObject function.
/// @author xingjiepan (xingjiepan@gmail.com)

#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/bridgeObjects_nonredundant.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// STD
#include <iostream>


namespace numeric {
namespace kinematic_closure {


void bridgeObjects_nonredundant(const utility::vector1<utility::fixedsizearray1<Real,3> >& stub1,
	const utility::vector1<utility::fixedsizearray1<Real,3> >& stub2,
	const utility::vector1<numeric::Real> & torsions_chain1, const utility::vector1<numeric::Real> & torsions_chain2,
	const utility::vector1<numeric::Real> & angles, const utility::vector1<numeric::Real> & bonds,
	utility::vector1<utility::vector1<Real> >& pivot_torsions,
	int &nsol){
	// Count the number of atoms

	unsigned int len_chain1 = torsions_chain1.size() + 3; // Number of atoms on chain1, including pivot1 and pivot2
	unsigned int len_chain2 = torsions_chain2.size() + 3; // Number of atoms on chain1, including pivot2 and pivot3
	unsigned int len_total = len_chain1 + len_chain2 + 3; // Number of all atoms
	int pivot1 = 3 + 2; // The original bridgeObjects function assumes the first pivot is 5, so I have add 2 dummy atoms
	int pivot2 = 2 + len_chain1 + 2;
	int pivot3 = len_total - 2 + 2;

	assert( angles.size() == len_total - 4);
	assert( bonds.size() == len_total - 5);

	// Create inputs for the bridgeObjects function

	utility::vector1<utility::fixedsizearray1<numeric::Real,3> > atoms(len_total + 2);
	utility::vector1<numeric::Real> dt(len_total + 2);
	utility::vector1<numeric::Real> da(len_total + 2);
	utility::vector1<numeric::Real> db(len_total + 2);
	utility::vector1<int> pivots(3);
	utility::vector1<int> order(3);

	// Assign some dummy postions to the first two atoms
	atoms[1] = stub1[1];
	atoms[1][1] += 1.0;
	atoms[2] = stub1[1];
	atoms[2][2] += 1.0;

	for ( unsigned int i=1; i<=3; ++i ) {
		atoms[pivot1 - 3 + i] = stub1[i];
		atoms[pivot3 - 1 + i] = stub2[i];
	}

	for ( unsigned int i=1; i<=torsions_chain1.size(); ++i ) {
		dt[pivot1 + i] = torsions_chain1[i];
	}

	for ( unsigned int i=1; i<=torsions_chain2.size(); ++i ) {
		dt[pivot2 + i] = torsions_chain2[i];
	}

	for ( unsigned int i=1; i<=angles.size(); ++i ) {
		da[pivot1 - 1 + i] = angles[i];
	}

	for ( unsigned int i=1; i<=bonds.size(); ++i ) {
		db[pivot1 - 1 + i] = bonds[i];
	}

	pivots[1] = pivot1;
	pivots[2] = pivot2;
	pivots[3] = pivot3;
	for ( unsigned int i=1; i<=3; ++i ) {order[i] = i;}

	// Create containers for the outputs from the bridgeObjects function

	utility::vector1<utility::vector1<numeric::Real> > t_ang;
	utility::vector1<utility::vector1<numeric::Real> > b_ang;
	utility::vector1<utility::vector1<numeric::Real> > b_len;

	// Solve the problem by calling the original bridgeObjects function

	bridgeObjects(atoms, dt, da, db, pivots, order, t_ang, b_ang, b_len, nsol);

	// Reorgnize the outputs

	pivot_torsions.clear();
	pivot_torsions.resize(nsol);

	for ( int i=1; i <= nsol; ++i ) {
		for ( int j=1; j <= 3; ++j ) {
			pivot_torsions[i].push_back(t_ang[i][pivots[j]-1]);
			pivot_torsions[i].push_back(t_ang[i][pivots[j]]);
		}
	}

}


} //numeric
} //kinematic_closure


