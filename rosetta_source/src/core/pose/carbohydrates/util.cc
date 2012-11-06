// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/carbohydrates/util.cc
/// @brief   Utility function definitions for carbohydrate-containing poses..
/// @author  labonte

// Unit headers
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/types.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>


// Construct tracer.
static basic::Tracer TR("core.pose.carbohydrates.util");


namespace core {
namespace pose {
namespace carbohydrates {

using namespace std;
using namespace core;


// Calculate and return the phi angle between a saccharide residue of the given pose and the previous residue.
/// @details This special-case function for carbohydrate phis is necessary, because of the following:\n
/// For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...\n
/// However, for aldopyranoses, BB X+1 is defined as: CX-OX-UPPER1-UPPER2.\n
/// CHI 1 is O5-C1-O1-HO1, which for an internal residue with virtual atoms for O1 and HO1, is the same as phi(n),
/// provided the virtual atoms are made to move with any rotation of BB X+1.  The same concept holds for aldo-
/// furanoses; however, ketoses are more complicated.  The cyclic oxygen must be the reference for phi, yet CHI 2 at
/// the anomeric position is defined with C1 as the reference atom, not the cyclic oxygen (O5 for furanoses, O6 for
/// pyranoses).  To complicate matters further, two virtual atoms in a row in a CHI gives NAN, so CHI angles cannot be
/// used after all.  Thus, we have to use vector calculus to calculate phi....
core::Angle
calculate_carbohydrate_phi(Pose const & pose, uint const sequence_position) {
	using namespace numeric;
	using namespace conformation;

	if (sequence_position == 1) {
		bool is_1st_residue_of_branch = false;  // TEMP
		if (is_1st_residue_of_branch) {
			// TODO: When branching is implemented, this must return the 1st torsion angle back to the main chain.
			return 180.0;  // TEMP
		} else {
			TR.Warning << "Phi is undefined for polysaccharide residue 1 unless part of a branch; "
					"returning 0.0." << endl;
			return 0.0;
		}
	}

	// Get two residues of interest.
	ResidueCAP res_n = & pose.residue(sequence_position);
	ResidueCAP res_n_minus_1 = & pose.residue(sequence_position - 1);

	// Get reference atom numbers.
	uint cyclic_O_num;
	if (res_n->carbohydrate_info()->is_aldose()) {
		cyclic_O_num = res_n->carbohydrate_info()->ring_size() - 1;
	} else /*is ketose*/ {
		cyclic_O_num = res_n->carbohydrate_info()->ring_size();
	}
	uint anomeric_C_num = res_n->carbohydrate_info()->anomeric_carbon();
	uint x = res_n_minus_1->carbohydrate_info()->mainchain_glycosidic_bond_acceptor();

	// Set the atom names of the four reference atoms.
	string O_cyclic = "O" + string(1, cyclic_O_num + '0');
	string C_anomeric = "C" + string(1, anomeric_C_num + '0');
	string OX = "O" + string(1, x + '0');
	string CX = "C" + string(1, x + '0');

	// Obtain the position vectors (a, b, c, d) of the four reference atoms.
	Vector a = res_n->xyz(O_cyclic);
	Vector b = res_n->xyz(C_anomeric);
	Vector c = res_n_minus_1->xyz(OX);
	Vector d = res_n_minus_1->xyz(CX);

	return dihedral_degrees(a, b, c, d);
}

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core
