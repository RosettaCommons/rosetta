// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// /// @file    numeric/geometry/projection_area.hh
// /// @brief   Function to get an approximate projection area.
// /// @author  SM Bargeen Alam Turzo  <turzo.1@osu.edu>

// Unit header
#include <numeric/geometry/projection_area.hh>

// Numeric headers
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/MathMatrix.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#define HEAVY_ATOM_EFFECTIVE_RADIUS 1.91 // effective vdw radius for heavy atoms

namespace numeric {
namespace geometry {
Real projection_area(utility::vector1 < Real > const &xs, utility::vector1 < Real > const &ys, utility::vector1 < Real > const &elements_rad, Real const probe_radius)
{ // Takes vectors of x-coordinates, y-coordinates, and element radius
	//Real probe_radius = probe_radius; // probe radius variable
	Real max_atom_radius = HEAVY_ATOM_EFFECTIVE_RADIUS; // max atom radius variable is set to 1.91 (Recall effective radius for C N O S P from Impact is 1.91 and hydrogen is 1.20)
	Real x_min = floor(*std::min_element( xs.begin(), xs.end() ) - max_atom_radius - probe_radius); // Getting x_min from pdb
	Real x_max = ceil(*std::max_element( xs.begin(), xs.end() ) + max_atom_radius + probe_radius); // Getting x_max from pdb
	utility::vector1 < Real > xs_new; // Creating a new xs list so that x_min can be set to zero.


	for ( Size i = 1; i <= xs.size(); ++i ) {
		xs_new.push_back(xs[i]-x_min);
	}

	Real x_max_new = x_max -x_min; // Since x_min is zero, x_max_new is subtracted from the old x_min.
	Real y_min = floor(*std::min_element( ys.begin(), ys.end() ) - max_atom_radius - probe_radius); // Same thing is done for y-coords of pdb
	Real y_max = ceil(*std::max_element( ys.begin(), ys.end() ) + max_atom_radius + probe_radius);
	utility::vector1 < Real > ys_new;
	for ( Size i = 1; i <= ys.size(); ++i ) {
		ys_new.push_back(ys[i]-y_min);
	}

	Real y_max_new = y_max - y_min;

	MathMatrix< Size > Grid((Size)(x_max_new)+10,(Size)(y_max_new)+10); // Creating a grid from x_max_new and y_max_new plus 10A extra in both dimensions because Steffen likes it.

	Size Grid_rows = Grid.get_number_rows(); // Getting the number of grid rows, required later when counting 1s.

	Real sqrt2inverse = 0.70710678118;
	for ( Size i =1; i <= xs.size(); ++i ) { // Looping through the list x-coordinates
		Real element_vdw_r = elements_rad[i] + probe_radius;
		Real diag_dist = element_vdw_r * sqrt2inverse;
		Grid( (Size)(xs_new[i])+5 , (Size)(ys_new[i])+5 ) = 1;
		Grid( (Size)(xs_new[i]+element_vdw_r)+5 , (Size)(ys_new[i])+5 )= 1;
		Grid( (Size)(xs_new[i]+diag_dist)+5 , (Size)(ys_new[i]+diag_dist)+5 )= 1;
		Grid( (Size)(xs_new[i])+5, (Size)(ys_new[i]+element_vdw_r)+5 )= 1;
		Grid( (Size)(xs_new[i]-diag_dist)+5 , (Size)(ys_new[i]+diag_dist)+5 )= 1;
		Grid( (Size)(xs_new[i]-element_vdw_r)+5 , (Size)(ys_new[i])+5 )= 1;
		Grid( (Size)(xs_new[i]-diag_dist)+5 , (Size)(ys_new[i]-diag_dist)+5 )= 1;
		Grid( (Size)(xs_new[i])+5 , (Size)(ys_new[i]-element_vdw_r)+5 )= 1;
		Grid( (Size)(xs_new[i]+diag_dist)+5 , (Size)(ys_new[i]-diag_dist)+5 )= 1;
	}
	Size area_from_row = 0;
	for ( Size i=0; i<Grid_rows; ++i ) { // Use total number of rows to count 1s
		MathVector< Size > each_row = Grid.get_row(i);
		for ( const Size *it = each_row.begin(); it != each_row.end(); ++it ) {
			// std::cout << *it;
			area_from_row += *it;
		}
		//std::cout<< "   | area from row = "<< area_from_row << "  |" << std::endl;
	}
	return area_from_row; // Return counted 1s which is the area of the shadow.
}
} // namespace geometry
} // namespace numeric
