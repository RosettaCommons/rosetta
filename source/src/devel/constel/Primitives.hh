// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Declarations of primitive functions and classes used by the constel program.
/// @author jk
/// @author Andrea Bazzoli

#ifndef INCLUDED_Primitives_hh
#define INCLUDED_Primitives_hh

#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/chemical/AA.hh>

using core::Size;
using core::pose::Pose;


namespace devel {
namespace constel {


/// @brief Returns the residue number of a residue in a pose.
core::Size get_pose_resnum(int const pdbnum, char const pdbchn, Pose& ps);


/// @brief A class to represent the mutation of a residue.
///
struct ResMut {

	char saa; // start amino acid type
	char eaa; // end amino acid type
	char cid; // identifier of the chain that the residue belongs to
	int pdbn; // pdb number of the residue in the chain
	Size psn; // number of the residue in its Rosetta pose

	ResMut(char s, char e, char c, int d, Size r) : saa(s), eaa(e), cid(c), pdbn(d), psn(r) {}
};

/// @brief Outputs all pair-constellations between a given pair of residues
void pair_constel_set_idx2(Size const i, Size const j, Pose const& pose_init );

/// @brief Outputs to file a constellation obtained from mutating a pair of residues
void out_pair_constel(ResMut const& mut1, ResMut const& mut2, int const cslnum,
	Pose& ps);

/// @brief Outputs all triple-constellations among a given triple of residues
void triple_constel_set_idx3(Size const i, Size const j, Size const k,
	Pose const& pose_init );

/// @brief Outputs to file a constellation obtained from mutating a triple of residues
void out_triple_constel(ResMut const& mut1, ResMut const& mut2,
	ResMut const& mut3, int const cslnum, Pose& ps);


/// @brief A class to hold data structures and functions shared by filters that
///  consider hydrogen bonding.
///
class HBondCommon {

public:

	/// @brief common database for the computation of hydrogen bonds
	static core::scoring::hbonds::HBondDatabaseCOP hb_database;

	static void init();

	/// @brief Tells whether a residue's moiety forms hydrogen bonds.
	static bool is_rmoi_hbonded(Pose const& ps, utility::vector1<Size> const& cnl,
		Size const im, bool const is_donor,
		utility::vector1<std::string> const& hb_atoms);
};


/// @brief A class to hold data structures and functions common to filters that
///  consider how groups of atoms within a constellation are oriented.
///
class OrientCommon {

public:

	/// @brief Tells whether atom 'low' is closer to atom 'tgt' than atom 'hi' is.
	static bool is_closer_to_tgt(numeric::xyzVector<core::Real> const& low,
		numeric::xyzVector<core::Real> const& hi,
		numeric::xyzVector<core::Real> const& tgt);
};


/// @brief A class to hold data structures and functions shared by filters that
///  check whether the amino acids and atoms they require in the constellation
///  are indeed present.
///
///
class PresenceCommon {

public:

	/// @brief Records the presence of given amino acid types in a given constellation.
	static bool are_aa_pres(core::pose::Pose const& ps,
		utility::vector1<Size> const& cnl,
		utility::vector1<core::chemical::AA> const& aa_typs,
		utility::vector1<Size>& aa_idxs);

	/// @brief Records the presence of given atoms in a given residue.
	static bool are_atoms_pres(core::conformation::Residue const& res,
		utility::vector1<std::string> const& anams,
		utility::vector1<Size>& aidxs);
};


} // constel
} // devel

#endif
