// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Motif.hh
/// @brief Header for interaction motif between residues
/// @author havranek, sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_Motif_hh
#define INCLUDED_protocols_motifs_Motif_hh

// Unit Headers
#include <protocols/motifs/Motif.fwd.hh>

// Package Headers

// Project Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <iosfwd>
#include <map>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

class Motif : public utility::pointer::ReferenceCount
{

static std::map < std::string, utility::vector1 < std::string > > motifAtomIDs;
static std::map < std::string, utility::vector1 < std::string > > basebaseAtomIDs;

public:
	// This is a 'by hand' way to associate a jump with the coordinate system-forming
	// atoms.
	Motif(
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
	Motif(
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

	// This constructor gets the jump for you - Rosetta numbering
	Motif(
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
	Motif(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	);

	// Ligand motif constructor: This constructor gets the jump for you - Rosetta numbering
	Motif(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		utility::vector1< Size >  const res2_atoms
	);

// Ligand motif search constructor: No residues, just atom names and a jump
	Motif(
		std::string const resname1,
		std::string const res1_atom1,
		std::string const res1_atom2,
		std::string const res1_atom3,
		std::string const res2_atom1,
		std::string const res2_atom2,
		std::string const res2_atom3,
		core::kinematics::Jump const & orientation
	);

	// Copy constructor
	Motif(
		Motif const & src
	);

	// Clone method
	MotifOP
	clone() const;

	// Destructor
	virtual ~Motif();

	// Residue check_res is res1 of motif
	virtual bool
	forward_check(
		core::conformation::Residue const & check_res
	) const;

	// Residue check_res is res2 of motif
	virtual bool
	backward_check(
		core::conformation::Residue const & check_res
	) const;

	// Position in pose is same type as motif
	virtual bool
	apply_check(
		core::pose::Pose const & pose,
		Size const pos
	) const;

	// Make atom integers for this motif, unused as of now, may want in future
	//virtual void
	//generate_atom_ints();

	virtual core::pack::rotamer_set::RotamerSetOP
	build_rotamers(
		core::pose::Pose & pose,
		Size const rotamer_build_position,
		Size const ex_,
		bool res2 = false
	) const;

	// Build a rotamer library using this motif and one given position to serve as the start
	// point for the motif.  The third, optional argument is the position at which to build
	// the rotamers.  If neglected, finds the closest position in the pose and uses that.
	virtual core::pack::rotamer_set::RotamerSetOP
	build_inverted_rotamers(
		core::pose::Pose & pose,
		Size const motif_anchor_position,
		bool & use_forward,
		Size rotamer_build_position = 0
	) const;

	// The place_residue function moves the mobile residue based on motif geometry
	virtual void
	place_residue(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool one_three = true
	) const;

	virtual void
	place_atoms(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		utility::vector1< std::string > const & atoms,
		bool one_three = true
	) const;

	virtual void
	place_atom(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		core::conformation::Atom & atm,
		bool one_three = true,
		std::string const atomtype = "C1'"
	) const;

	virtual void
	place_residue_(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool forward,
		bool one_three = true
	) const;

	virtual void
	place_atom_(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool forward,
		core::conformation::Atom & atm,
		bool one_three = true,
		std::string const atomtype = "C1'"
	) const;

	virtual void
	place_atoms_(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool forward,
		utility::vector1< std::string > const & atoms,
		bool one_three = true
	) const;

	// Ligand motifs place_residue
	virtual void
	place_residue(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		Size const & res2_atom1_index_in,
		Size const & res2_atom2_index_in,
		Size const & res2_atom3_index_in,
		bool one_three = true
	) const;

	// Ligand motifs place_atoms
	virtual void
	place_atoms(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		utility::vector1< Size > const & atoms,
		Size const & res2_atom1_index_in,
		Size const & res2_atom2_index_in,
		Size const & res2_atom3_index_in,
		bool one_three = true
	) const;

	// Ligand motifs place_atom
	virtual void
	place_atom(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		core::conformation::Atom & atm,
		Size const & res2_atom1_index_in,
		Size const & res2_atom2_index_in,
		Size const & res2_atom3_index_in,
		Size const & atomtype,
		bool one_three = true
	) const;

	// Ligand motifs place_residue_
	virtual void
	place_residue_(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool forward,
		Size const & res2_atom1_index_in,
		Size const & res2_atom2_index_in,
		Size const & res2_atom3_index_in,
		bool one_three = true
	) const;

	// Ligand motifs place_atom_
	virtual void
	place_atom_(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool forward,
		core::conformation::Atom & atm,
	  Size const & res2_atom1_index_in,
	  Size const & res2_atom2_index_in,
    Size const & res2_atom3_index_in,
		Size const & atomtype,
		bool one_three = true
	) const;

	// Ligand motifs place_atoms_
	virtual void
	place_atoms_(
		core::conformation::Residue const & fixed,
		core::conformation::Residue & mobile,
		bool forward,
		utility::vector1< Size > const & atoms,
	  Size const & res2_atom1_index_in,
	  Size const & res2_atom2_index_in,
    Size const & res2_atom3_index_in,
		bool one_three = true
	) const;

	// For output
	virtual void
	print(
		std::ostream & os
	) const;

	// For output
	virtual
	std::string print() const;

	// Overloaded operator for output
	friend std::ostream & operator <<(
		std::ostream & os,
		Motif const & mot
	);

	// Stores remark, havranek uses
	void
	store_remark(
		std::string const & remark_in
	);

	// Stores path of motif for future output if motif was read from pdb
	void
	store_path(
		std::string const & path_in
	);

	// Accessors
	std::string const & restype_name1() const { return restype_name1_; }
	std::string const & res1_atom1_name() const { return res1_atom1_name_; }
	std::string const & res1_atom2_name() const { return res1_atom2_name_; }
	std::string const & res1_atom3_name() const { return res1_atom3_name_; }

	/// @brief WARNING res*_atom*_int() values are not consistently initialized in constructors
	int const & res1_atom1_int() const { return res1_atom1_int_; }
	int const & res1_atom2_int() const { return res1_atom2_int_; }
	int const & res1_atom3_int() const { return res1_atom3_int_; }
	Size const & res1_atom1_index() const { return res1_atom1_index_; }
	Size const & res1_atom2_index() const { return res1_atom2_index_; }
	Size const & res1_atom3_index() const { return res1_atom3_index_; }
	std::string const & restype_name2() const { return restype_name2_; }
	std::string const & res2_atom1_name() const { return res2_atom1_name_; }
	std::string const & res2_atom2_name() const { return res2_atom2_name_; }
	std::string const & res2_atom3_name() const { return res2_atom3_name_; }

	/// @brief WARNING res*_atom*_int() values are not consistently initialized in constructors
	int const & res2_atom1_int() const { return res2_atom1_int_; }
	int const & res2_atom2_int() const { return res2_atom2_int_; }
	int const & res2_atom3_int() const { return res2_atom3_int_; }
	Size const & res2_atom1_index() const { return res2_atom1_index_; }
	Size const & res2_atom2_index() const { return res2_atom2_index_; }
	Size const & res2_atom3_index() const { return res2_atom3_index_; }
	core::kinematics::Jump const & forward_jump() const { return forward_jump_; }
	core::kinematics::Jump const & backward_jump() const { return backward_jump_; }
	bool has_remark() const { return has_remark_;}
	bool has_path() const { return has_path_;}
	std::string const & remark() const { return remark_;}
	std::string const & path() const { return path_;}

private:

	// *********************************************************
	// IMPORTANT: Motif has *A LOT* of constructors.
	// If you add a new data members be sure you initialize them
	// in each and *every* constructor
	// *********************************************************

	// names (name3) of the residues in motif and each of the 6 motif atoms
	std::string restype_name1_;
	std::string res1_atom1_name_;
	std::string res1_atom2_name_;
	std::string res1_atom3_name_;
	int res1_atom1_int_;
	int res1_atom2_int_;
	int res1_atom3_int_;
	Size res1_atom1_index_;
	Size res1_atom2_index_;
	Size res1_atom3_index_;
	std::string restype_name2_;
	std::string res2_atom1_name_;
	std::string res2_atom2_name_;
	std::string res2_atom3_name_;
	int res2_atom1_int_;
	int res2_atom2_int_;
	int res2_atom3_int_;
	Size res2_atom1_index_;
	Size res2_atom2_index_;
	Size res2_atom3_index_;
	// Jumps associated with the motif
	core::kinematics::Jump forward_jump_;
	core::kinematics::Jump backward_jump_;
	// Associated information for future output
	std::string remark_;
	std::string path_;
	bool has_remark_;
	bool has_path_;

};

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_Motif
