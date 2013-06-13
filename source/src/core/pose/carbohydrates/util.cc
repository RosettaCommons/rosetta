// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/carbohydrates/util.cc
/// @brief   Utility function definitions for carbohydrate-containing poses.
/// @author  Labonte

// Unit headers
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/Pose.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/types.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

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

// Scan through a saccharide residue's connections to find the residue from which it follows or branches.
/// @return  The sequence position of the residue before this one (n-1) or the residue in the parent chain from which
/// the branch occurs or zero if N/A, i.e., if this is the lower terminus.
core::uint
find_seqpos_of_parent_residue(conformation::Residue const & residue) {
	uint seqpos = residue.seqpos();
	Size n_connections = residue.n_residue_connections();

	// Search backwards for speed, since "non-polymeric" connections come after polymeric ones.
	for (uint i = n_connections; i >= 1; --i) {
		uint connection_partner_position = residue.residue_connection_partner(i);
		if (connection_partner_position < seqpos) {  // assumes PDB places branches after main chains
			return connection_partner_position;
		}
	}
	TR.Warning << "This residue is a lower terminus! Returning 0." << endl;
	return 0;
}

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
	using namespace utility;
	using namespace conformation;

	// Get the 1st residue of interest.
	ResidueCAP res_n = & pose.residue(sequence_position);

	if (res_n->is_lower_terminus()) {
		TR.Warning << "Phi is undefined for the first polysaccharide residue of a chain unless part of a branch; "
				"returning 0.0." << endl;
		return 0.0;
	}

	// Get the 2nd residue of interest.
	// (res_n_minus_1 is a misnomer for the lower termini of branches.)
	ResidueCAP res_n_minus_1 = & pose.residue(find_seqpos_of_parent_residue(*res_n));


	// Set the atom names of the four reference atoms.
	string ref1;  // O(cyclic) for cyclic saccharides; C? for linear saccharides
	if (res_n->carbohydrate_info()->is_cyclic()) {
		ref1 = res_n->carbohydrate_info()->cyclic_oxygen_name();
	} else /* is linear */ {
		ref1 = "C?";  // TODO: Figure out how linear polysaccharides are handled by IUPAC.
	}

	string ref2 = res_n->carbohydrate_info()->anomeric_carbon_name();  // always the anomeric carbon

	string ref3 = res_n_minus_1->atom_name(res_n_minus_1->connect_atom(*res_n));  // OX(n-1) for polysaccharides

	string ref4;  // CX(n-1) for polysaccharides; an AA carbon for the 1st saccharide residue of a glycopeptide
	if (res_n_minus_1->is_carbohydrate()) {
		uint x = atoi(&ref3[2]);  // 3rd column (index 2) is the atom number
		ref4 = "C" + string(1, x + '0');  // "CX"
	} else /* is amino acid residue */ {
		// Get list of indices of all atoms connected to reference atom 3.
		vector1<uint> atom_indices = res_n_minus_1->bonded_neighbor(res_n_minus_1->atom_index(ref3));

		// Search for heavy atoms.  Since a residue connection is not an atom, it must be the next reference atom.
		Size n_indices = atom_indices.size();
		for (uint i = 1; i < n_indices; ++i) {
			if (!res_n_minus_1->atom_is_hydrogen(atom_indices[i])) {
				ref4 = res_n_minus_1->atom_name(atom_indices[i]);
				break;
			}
		}
	}


	// Obtain the position vectors (a, b, c, d) of the four reference atoms.
	Vector a = res_n->xyz(ref1);
	Vector b = res_n->xyz(ref2);
	Vector c = res_n_minus_1->xyz(ref3);
	Vector d = res_n_minus_1->xyz(ref4);

	return dihedral_degrees(a, b, c, d);
}

// Return the number of degrees by which the phi angle between a saccharide residue of the given pose and the previous
// residue differs from the BB torsion used by Rosetta.
/// @remarks See the details for calculate_carbohydrate_phi() for an explanation on why this method is necessary.
core::Angle
carbohydrate_phi_offset_from_BB(Pose const & pose, uint const sequence_position)
{
	using namespace id;

	// Get the actual value of phi.
	Angle actual_phi = calculate_carbohydrate_phi(pose, sequence_position);

	// TODO: Refactor completely

	// Get the appropriate BB torsion (found on the previous residue).
	uint x = pose.residue_type(sequence_position - 1).carbohydrate_info()->mainchain_glycosidic_bond_acceptor();
	Angle bb_torsion = pose.torsion(TorsionID(sequence_position - 1, BB, x + 1));

	// Return the difference.
	return actual_phi - bb_torsion;
}

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core
