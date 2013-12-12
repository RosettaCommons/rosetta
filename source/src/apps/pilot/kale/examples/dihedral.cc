// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

// Usually you wouldn't need to use either of these typedefs, because they are 
// already defined in the core library.  Here I'm just trying to avoid using 
// the core library, to simplify the example.

typedef numeric::Real Real;
typedef numeric::xyzVector<Real> Vector;

int main(int argc, char* argv[]) {

	Vector a, b, c, d;
	Real angle;

	// These coordinates were taken from the omega angle of an ideal alpha helix.  
	// The dihedral angle of these four points should be near 180.
	
	a = Vector(-1.165, -3.357, -1.789);
	b = Vector(-0.641, -1.928, -1.789);
	c = Vector(-1.201, -1.087, -2.662);
	d = Vector(-0.790,  0.300, -2.759);

	angle = dihedral(a, b, c, d);

	std::cout << angle <<std::endl;
}


