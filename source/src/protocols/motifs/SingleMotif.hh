// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SingleMotif.hh
/// @brief Header for interaction motif between residues
/// @author havranek, sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_SingleMotif_hh
#define INCLUDED_protocols_motifs_SingleMotif_hh

// Unit Headers
#include <protocols/motifs/SingleMotif.fwd.hh>

// Package Headers
#include <protocols/motifs/Motif.hh>

#include <utility/vector1.hh>


// Project Headers

// Utility Headers

// C++ Headers

namespace protocols {
namespace motifs {

class SingleMotif : public Motif {

public:
	// This is a 'by hand' way to associate a jump with the coordinate system-forming
	// atoms.
	SingleMotif(
		std::string const resname1,
		std::string const res1_atom1,
		std::string const res1_atom2,
		std::string const res1_atom3,
		std::string const resname2,
		std::string const res2_atom1,
		std::string const res2_atom2,
		std::string const res2_atom3,
		core::kinematics::Jump const & orientation
	);

	// This constructor gets the jump for you - PDB numbering + chain id
	SingleMotif(
		core::pose::Pose const & pose,
		Size const pdb_residue_position_1,
		char const pdb_chain_id1,
		std::string const res1_atom1,
		std::string const res1_atom2,
		std::string const res1_atom3,
		Size const pdb_residue_position_2,
		char const pdb_chain_id2,
		std::string const res2_atom1,
		std::string const res2_atom2,
		std::string const res2_atom3
	);

	// This constructor takes a jump, used for ligand motifs
	SingleMotif(
		std::string const resname1,
		std::string const res1_atom1,
		std::string const res1_atom2,
		std::string const res1_atom3,
		std::string const res2_atom1,
		std::string const res2_atom2,
		std::string const res2_atom3,
		core::kinematics::Jump const & orientation
	);

	// This constructor gets the jump for you - Rosetta numbering
	SingleMotif(
		core::pose::Pose const & pose,
		Size const pdb_residue_position_1,
		std::string const res1_atom1,
		std::string const res1_atom2,
		std::string const res1_atom3,
		Size const pdb_residue_position_2,
		std::string const res2_atom1,
		std::string const res2_atom2,
		std::string const res2_atom3
	);

	// This constructor forms a motif using preset atom types from a static map
	SingleMotif(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	);

	// Copy constructor
	SingleMotif(
		SingleMotif const & src
	);

	// Overloaded operator for output
	friend std::ostream & operator <<(
		std::ostream & os,
		SingleMotif const & mot
	);

};

}
}

#endif // INCLUDED_protocols_motifs_SingleMotif
